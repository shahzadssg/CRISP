"""
CRISP_protocol_faithful.py
==========================

Protocol-faithful implementation of CRISP: Circuit-pRivate Single-Image
Steganography with Permutations.

Carol never sees (row, col).  Wire values flow from
gate to gate as IMAGES, not as scalars.  Each wire is identified by
(image, channel): the image carries the wire's bit in the LSB of a
specific channel at every pixel — including, but not specially, at
the secret pixel.  Routing one wire to another wire's logical role at
a downstream gate is a per-pixel channel-copy operation that does not
depend on (row, col) and processes every pixel uniformly.

Bob is the only party that calls ext().  He calls it exactly once per
primary output, on the image that carries it.

Public objects (Carol sees these):
  - all input covers (with primary inputs and ancillaries embedded)
  - all output covers (fresh, before computation)
  - the circuit specification (gates in random topological order)
  - per-gate (pi_in, pi_out) ∈ S3 × S3

Secret objects (only Alice and Bob):
  - (row, col)
"""

from __future__ import annotations
import itertools
import random
import time
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
ALL_PERMS: List[Tuple[int, int, int]] = list(itertools.permutations([0, 1, 2]))
GATE_TYPES = {"FREDKIN": 3}


# ===========================================================================
# Circuit data structure
# ===========================================================================
class Circuit:
    def __init__(self):
        self.inputs:      List[str]                  = []
        self.ancillaries: List[Tuple[str, int]]      = []
        self.gates:       List[Dict]                 = []
        self.outputs:     List[str]                  = []

    def add_input(self, name):            self.inputs.append(name)
    def add_ancillary(self, name, value): self.ancillaries.append((name, value))
    def add_gate(self, t, ins, outs):     self.gates.append(
        {"type": t, "inputs": ins, "outputs": outs})
    def add_output(self, name):           self.outputs.append(name)


def parse_circuit(spec: str) -> Circuit:
    circuit = Circuit()
    for line in spec.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tokens = line.split()
        if tokens[0] == "INPUT":
            for name in "".join(tokens[1:]).split(","):
                circuit.add_input(name.strip())
        elif tokens[0] == "ANCILLARY":
            for part in "".join(tokens[1:]).split(","):
                wire, value = part.split("=")
                circuit.add_ancillary(wire.strip(), int(value.strip()))
        elif tokens[0] in GATE_TYPES:
            gate_type   = tokens[0]
            expected_io = GATE_TYPES[gate_type]
            arrow_idx   = tokens.index("->")
            ins  = [x.strip() for x in "".join(tokens[1:arrow_idx]).split(",")]
            outs = [x.strip() for x in "".join(
                tokens[arrow_idx + 1:]).split(",")]
            if len(ins) != expected_io or len(outs) != expected_io:
                raise ValueError(
                    f"{gate_type} gate must have {expected_io} inputs/outputs")
            circuit.add_gate(gate_type, ins, outs)
        elif tokens[0] == "OUTPUT":
            for name in "".join(tokens[1:]).split(","):
                circuit.add_output(name.strip())
    return circuit


# ===========================================================================
# Randomised topological sort
# ===========================================================================
def _build_dependency_graph(gates):
    graph        = defaultdict(list)
    in_degree    = defaultdict(int)
    wire_to_gate = defaultdict(list)
    for i, gate in enumerate(gates):
        for out in gate["outputs"]:
            wire_to_gate[out].append(i)
    for i, gate in enumerate(gates):
        for inp in gate["inputs"]:
            for producer in wire_to_gate[inp]:
                graph[producer].append(i)
                in_degree[i] += 1
    sources = [i for i in range(len(gates)) if in_degree[i] == 0]
    return graph, in_degree, sources


def randomized_topological_sort(gates) -> List[int]:
    graph, in_degree, sources = _build_dependency_graph(gates)
    order: List[int] = []
    while sources:
        gate_idx = random.choice(sources)
        order.append(gate_idx)
        sources.remove(gate_idx)
        for dep in graph[gate_idx]:
            in_degree[dep] -= 1
            if in_degree[dep] == 0:
                sources.append(dep)
    if len(order) != len(gates):
        raise ValueError("Circuit contains a cycle")
    return order


# ===========================================================================
# Steganography primitives
# ===========================================================================
def fredkin(c: int, x: int, y: int) -> Tuple[int, int, int]:
    return (c, y, x) if c == 1 else (c, x, y)


def perm_gen() -> Tuple[Tuple[int, int, int], Tuple[int, int, int]]:
    """Sample two independent uniform-random permutations from S3."""
    return ALL_PERMS[random.randrange(6)], ALL_PERMS[random.randrange(6)]


def loc_gen(image_shape: Tuple[int, int, int]) -> Tuple[int, int]:
    """Alice samples (row, col) uniformly. The secret never leaves Alice/Bob."""
    h, w, _ = image_shape
    return random.randrange(h), random.randrange(w)


def emb_bit(image: np.ndarray,
            bit: int,
            row: int, col: int,
            channel: int) -> np.ndarray:
    """
    Alice: write a single bit into the LSB of (row, col, channel).

    All other pixels and channels are untouched; the rest of the image
    is the natural cover.  Returns a fresh stego image (input is not
    modified).
    """
    stego = image.copy()
    stego[row, col, channel] = (
        int(stego[row, col, channel]) & 0xFE) | int(bit)
    return stego


def comp(stego_in:     np.ndarray,
         output_cover: np.ndarray,
         pi_in:        Tuple[int, int, int],
         pi_out:       Tuple[int, int, int]) -> np.ndarray:
    """
    Carol: apply Fredkin pixel-wise across the entire image.

    Reads logical input role t from channel pi_in[t] of *stego_in* at
    every pixel.  Writes logical output role t into channel pi_out[t]
    of a fresh copy of *output_cover* at every pixel.  Carol processes
    every pixel identically and never sees (row, col).
    """
    src    = stego_in.astype(np.uint16)
    result = output_cover.copy().astype(np.uint16)
    mask   = np.uint16(0xFE)

    c_bits = src[:, :, pi_in[0]] & np.uint16(1)
    x_bits = src[:, :, pi_in[1]] & np.uint16(1)
    y_bits = src[:, :, pi_in[2]] & np.uint16(1)

    swap = (c_bits == 1)
    out  = [
        c_bits,
        np.where(swap, y_bits, x_bits).astype(np.uint16),
        np.where(swap, x_bits, y_bits).astype(np.uint16),
    ]
    for t in range(3):
        ch = pi_out[t]
        result[:, :, ch] = (result[:, :, ch] & mask) | out[t]
    return result.astype(np.uint8)


def ext_bit(image: np.ndarray,
            row: int, col: int,
            channel: int) -> int:
    """
    Bob: read the LSB at (row, col, channel).  Called only by Bob, only
    on circuit-output images, and only at primary-output extraction.
    """
    return int(image[row, col, channel]) & 1


# ===========================================================================
# Wire routing  (per-pixel, position-oblivious)
# ===========================================================================
#
# A wire is represented by a "channel image": an (h, w, 3) uint8 array
# whose channel `chan` carries the wire's bit at every pixel — including
# at the secret pixel — in the LSB.  The other two channels are filler
# from the carrier.
#
# When a gate G_i needs three input bits (c, x, y) with channel layout
# pi_in[0] = chan(c), pi_in[1] = chan(x), pi_in[2] = chan(y), Carol
# assembles a stego image whose channel pi_in[t] holds the corresponding
# wire's bit at every pixel.  This is a per-pixel channel-copy that
# does NOT depend on (row, col) — it operates uniformly on every pixel.
# ===========================================================================

def assemble_stego(input_cover:  np.ndarray,
                   wire_sources: List[Tuple[np.ndarray, int]],
                   pi_in:        Tuple[int, int, int]) -> np.ndarray:
    """
    Build the stego image for a gate by copying wire LSBs into the
    correct channels of a fresh input cover.

    Args:
        input_cover:   fresh (h, w, 3) cover for this gate
        wire_sources:  three (source_image, source_channel) pairs, one
                       per logical role (control, x, y)
        pi_in:         logical-role -> physical-channel assignment for
                       this gate

    Returns:
        stego image where for every pixel:
          channel pi_in[t]  has its LSB equal to the source LSB of role t
          all other bits     come from input_cover unchanged
    """
    stego = input_cover.copy().astype(np.uint16)
    mask  = np.uint16(0xFE)
    for t, (src_img, src_ch) in enumerate(wire_sources):
        target_ch     = pi_in[t]
        src_lsb_plane = src_img[:, :, src_ch].astype(np.uint16) & np.uint16(1)
        stego[:, :, target_ch] = (stego[:, :, target_ch] & mask) | src_lsb_plane
    return stego.astype(np.uint8)


# ===========================================================================
# Alice: prepare source images for primary inputs and ancillaries
# ===========================================================================
def alice_prepare_sources(circuit:     Circuit,
                          source_pool: List[np.ndarray],
                          row: int, col: int,
                          input_bits:    Dict[str, int]
                          ) -> Dict[str, Tuple[np.ndarray, int]]:
    """
    Alice creates one stego image per primary-input wire and per
    ancillary, embedding the bit at (row, col) in a fixed channel
    (channel 0 — irrelevant to security since the secret IS the
    position, and the channel is wire-specific).

    Returns a wire-to-(image, channel) map for all source wires.
    """
    sources: Dict[str, Tuple[np.ndarray, int]] = {}
    pool_idx = 0
    for wire in circuit.inputs:
        cover = source_pool[pool_idx]
        pool_idx += 1
        bit = input_bits[wire]
        # Channel 0 is conventional; pi_in at the consuming gate will
        # route it correctly.  No channel secrecy is claimed here.
        stego = emb_bit(cover, bit, row, col, channel=0)
        sources[wire] = (stego, 0)
    for wire, value in circuit.ancillaries:
        cover = source_pool[pool_idx]
        pool_idx += 1
        stego = emb_bit(cover, value, row, col, channel=0)
        sources[wire] = (stego, 0)
    return sources


# ===========================================================================
# Carol: protocol-faithful evaluator
# ===========================================================================
def carol_evaluate(circuit:        Circuit,
                   gate_order:     List[int],
                   wire_sources:   Dict[str, Tuple[np.ndarray, int]],
                   input_covers:   List[np.ndarray],
                   output_covers:  List[np.ndarray],
                   permutations:   List[Tuple[Tuple[int, int, int],
                                              Tuple[int, int, int]]]
                   ) -> Tuple[Dict[str, Tuple[np.ndarray, int]], List[np.ndarray]]:
    """
    Carol evaluates the circuit gate-by-gate.  She never sees (row, col).

    Args:
        wire_sources:  wire -> (image, channel) for primary inputs and
                       ancillaries (provided by Alice with bits embedded)
        input_covers:  one fresh cover per gate (in gate-execution order)
        output_covers: one fresh cover per gate (in gate-execution order)
        permutations:  (pi_in, pi_out) per gate (in gate-execution order)

    Returns:
        wire_to_image:  wire -> (image, channel) including outputs
        gate_outputs:   list of output images, one per gate, in execution order
    """
    wires = dict(wire_sources)
    gate_outputs: List[np.ndarray] = []

    for step, gate_idx in enumerate(gate_order):
        gate = circuit.gates[gate_idx]
        c_wire, x_wire, y_wire    = gate["inputs"]
        co_wire, xo_wire, yo_wire = gate["outputs"]

        sources = [wires[c_wire], wires[x_wire], wires[y_wire]]
        pi_in, pi_out = permutations[step]

        stego   = assemble_stego(input_covers[step], sources, pi_in)
        out_img = comp(stego, output_covers[step], pi_in, pi_out)
        gate_outputs.append(out_img)

        wires[co_wire] = (out_img, pi_out[0])
        wires[xo_wire] = (out_img, pi_out[1])
        wires[yo_wire] = (out_img, pi_out[2])

    return wires, gate_outputs


# ===========================================================================
# Bob: extract primary outputs
# ===========================================================================
def bob_extract(circuit:      Circuit,
                wire_to_image: Dict[str, Tuple[np.ndarray, int]],
                row: int, col: int) -> Dict[str, int]:
    """
    Bob extracts the primary output bits.  Called once per output wire.
    """
    out_bits: Dict[str, int] = {}
    for wire in circuit.outputs:
        img, ch = wire_to_image[wire]
        out_bits[wire] = ext_bit(img, row, col, ch)
    return out_bits


# ===========================================================================
# Top-level orchestration
# ===========================================================================
def run_crisp(circuit_spec:  str,
              input_bits:    Dict[str, int],
              image_shape:   Tuple[int, int, int] = (128, 128, 3),
              seed:          Optional[int] = None
              ) -> Tuple[Dict[str, int], Dict[str, float]]:
    """
    End-to-end execution.  Returns the primary-output bits and timing.

    The function plays all three roles for demonstration, but each
    party's view is honoured strictly:
      - Carol does not receive (row, col)
      - Bob receives (row, col) and the output images Carol produces
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    t_start = time.perf_counter()

    circuit    = parse_circuit(circuit_spec)
    n_gates    = len(circuit.gates)
    n_sources  = len(circuit.inputs) + len(circuit.ancillaries)
    h, w, _    = image_shape

    # ---- Alice ----
    row, col = loc_gen(image_shape)

    source_pool = [
        np.random.randint(0, 256, image_shape, dtype=np.uint8)
        for _ in range(n_sources)
    ]
    wire_sources = alice_prepare_sources(
        circuit, source_pool, row, col, input_bits)

    # ---- Public setup (visible to Carol) ----
    gate_order = randomized_topological_sort(circuit.gates)
    permutations = [perm_gen() for _ in range(n_gates)]
    input_covers = [
        np.random.randint(0, 256, image_shape, dtype=np.uint8)
        for _ in range(n_gates)
    ]
    output_covers = [
        np.random.randint(0, 256, image_shape, dtype=np.uint8)
        for _ in range(n_gates)
    ]

    t_carol_start = time.perf_counter()

    # ---- Carol (no row/col!) ----
    wire_to_image, gate_outputs = carol_evaluate(
        circuit, gate_order, wire_sources,
        input_covers, output_covers, permutations,
    )

    t_carol = time.perf_counter() - t_carol_start

    # ---- Bob ----
    out_bits = bob_extract(circuit, wire_to_image, row, col)

    t_total = time.perf_counter() - t_start

    timing = {
        "carol_seconds": t_carol,
        "carol_per_gate_seconds": t_carol / max(n_gates, 1),
        "total_seconds": t_total,
        "n_gates": n_gates,
    }
    return out_bits, timing


# ===========================================================================
# Tests
# ===========================================================================
def run_unit_tests(verbose: bool = False) -> bool:
    """
    Exhaustive single-gate tests: 8 input triples × 36 permutation
    pairs = 288 cases.  Verifies that ext(comp(emb(...))) recovers the
    Fredkin output for every (input triple, pi_in, pi_out, secret pos).
    """
    h, w   = 16, 16
    cover  = np.random.randint(0, 256, (h, w, 3), dtype=np.uint8)
    out_cv = np.random.randint(0, 256, (h, w, 3), dtype=np.uint8)
    src1   = np.random.randint(0, 256, (h, w, 3), dtype=np.uint8)
    src2   = np.random.randint(0, 256, (h, w, 3), dtype=np.uint8)
    src3   = np.random.randint(0, 256, (h, w, 3), dtype=np.uint8)
    passed = failed = 0

    for c, x, y in itertools.product([0, 1], repeat=3):
        expected = fredkin(c, x, y)
        for pi_in in ALL_PERMS:
            for pi_out in ALL_PERMS:
                row = random.randrange(h)
                col = random.randrange(w)

                # Each input bit lives in a different source image
                # to exercise full assemble_stego behaviour.
                src1_emb = emb_bit(src1, c, row, col, 0)
                src2_emb = emb_bit(src2, x, row, col, 0)
                src3_emb = emb_bit(src3, y, row, col, 0)
                wire_sources = [(src1_emb, 0), (src2_emb, 0), (src3_emb, 0)]

                stego  = assemble_stego(cover, wire_sources, pi_in)
                out_im = comp(stego, out_cv, pi_in, pi_out)
                got    = (
                    ext_bit(out_im, row, col, pi_out[0]),
                    ext_bit(out_im, row, col, pi_out[1]),
                    ext_bit(out_im, row, col, pi_out[2]),
                )
                if got == expected:
                    passed += 1
                else:
                    failed += 1
                    if verbose:
                        print(f"FAIL ({c},{x},{y}) pi_in={pi_in} "
                              f"pi_out={pi_out}: expected {expected} got {got}")
    print(f"Unit tests: {passed}/288 passed, {failed} failed")
    return failed == 0


# ---------------------------------------------------------------------------
# Benchmark circuit (matches the paper's implementation)
# Computes f(A,B,C) = (A AND C) OR (NOT A AND B) OR (NOT B AND NOT C)
# 14 Fredkin gates, 3 inputs, 2 ancillaries.
# ---------------------------------------------------------------------------
BENCHMARK_CIRCUIT = """
INPUT A, B, C
ANCILLARY const0 = 0
ANCILLARY const1 = 1
FREDKIN A, C, const0 -> A1, C1, AC
FREDKIN A, const1, const0 -> A2, notA, A3
FREDKIN notA, B, const0 -> notA1, B1, notAB
FREDKIN B, const1, const0 -> B2, notB, B3
FREDKIN C, const1, const0 -> C2, notC, C3
FREDKIN notB, notC, const0 -> notB1, notC1, notBnotC
FREDKIN AC, const1, const0 -> AC1, notAC, AC2
FREDKIN notAB, const1, const0 -> notAB1, not_notAB, notAB2
FREDKIN notAC, not_notAB, const0 -> notAC1, not_notAB1, not_T1
FREDKIN not_T1, const1, const0 -> not_T1_1, T1, not_T1_2
FREDKIN T1, const1, const0 -> T1_2, not_T1_3, T1_3
FREDKIN notBnotC, const1, const0 -> notBnotC2, not_notBnotC, notBnotC3
FREDKIN not_T1_3, not_notBnotC, const0 -> not_T1_4, not_notBnotC1, not_result
FREDKIN not_result, const1, const0 -> not_result_1, result, not_result_2
OUTPUT result
"""


def run_truth_table() -> bool:
    """Verify f(A,B,C) = (A&C) | (~A&B) | (~B&~C) for all 8 inputs."""
    all_ok = True
    for A, B, C in itertools.product([0, 1], repeat=3):
        expected = int((A & C) | ((1 - A) & B) | ((1 - B) & (1 - C)))
        bits = {"A": A, "B": B, "C": C}
        out, _ = run_crisp(BENCHMARK_CIRCUIT, bits)
        got = out["result"]
        ok  = (got == expected)
        all_ok = all_ok and ok
        print(f"A={A} B={B} C={C}: expected {expected} got {got}  "
              f"{'ok' if ok else 'FAIL'}")
    return all_ok


# ===========================================================================
# Main
# ===========================================================================
if __name__ == "__main__":
    print("=" * 64)
    print("CRISP — protocol-faithful implementation")
    print("Carol never receives (row, col). Bob extracts only.")
    print("=" * 64)

    print("\n[Unit tests]")
    ok1 = run_unit_tests()

    print("\n[Truth table on benchmark circuit]")
    ok2 = run_truth_table()

    print("\n[Timing sample]")
    _, timing = run_crisp(BENCHMARK_CIRCUIT, {"A": 1, "B": 0, "C": 1})
    print(f"  {timing}")

    print(f"\nAll tests passed: {ok1 and ok2}")
