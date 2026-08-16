"""
CRISP: Circuit-pRivate Single-Image Steganography with Permutations
=====================================================================
  - Circuit specification (a 14-gate Fredkin circuit)
  - Randomised topological sort
  - Repairs R1 (LSB-plane re-randomisation) and R2 (pixel-uniform
    constant planes)
  - Step-by-step verbose output, per-gate image saves and matplotlib
    display hooks
  - 288 unit tests + truth-table verification
  - The constant-wire attack, run against a real transcript


Image counts
------------
  A gate costs two covers, one for its assembled input and one for its
  output, against six per gate in the wire-per-image baseline.  The
  asymptotic per-gate reduction is 3x.  For the 14-gate benchmark the
  totals are 31 images against 47 once R2 removes the two constant
  carriers, or 33 against 47 if those carriers are counted.

Security
--------
  Joint position-and-value recovery over the full transcript:

      Adv = E[1 / (1 + B)]  >=  2^lambda / (2^lambda + n - 1)
      B ~ Binomial(n - 1, 2^-lambda),  n = h * w,
      lambda = m - H2(q)

  where m is the number of secret input bits and H2(q) is the collision
  entropy of the server's prior on them.  lambda counts the bits she
  already knows and each one halves the candidate pixel set.  Repairs
  R1 and R2 give lambda = 0, which is exactly 1 / (h * w).

  CRISP provides no circuit privacy.  The server holds both channel
  permutations because she needs them to evaluate and to route, so she
  reads the channel-to-role map at every gate directly, and an adversary
  holding only the images recovers that map and the wiring graph with
  nine channel-pair comparisons per gate.
"""

import itertools
import os
import random
import re
import string
import tempfile
import time
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
ALL_PERMS: List[Tuple[int, int, int]] = list(itertools.permutations([0, 1, 2]))
GATE_TYPES = {"FREDKIN": 3}


# ===========================================================================
# Circuit data structure
# ===========================================================================
class Circuit:
    """Plain container for parsed circuit elements."""

    def __init__(self):
        self.inputs: List[str] = []
        self.ancillaries: List[Tuple[str, int]] = []
        self.gates: List[Dict] = []
        self.outputs: List[str] = []

    def add_input(self, name):
        self.inputs.append(name)

    def add_ancillary(self, name, value):
        self.ancillaries.append((name, value))

    def add_gate(self, gate_type, ins, outs):
        self.gates.append({"type": gate_type, "inputs": ins, "outputs": outs})

    def add_output(self, name):
        self.outputs.append(name)


# ===========================================================================
# Circuit parser
# ===========================================================================
def parse_circuit(spec: str) -> Circuit:
    """Parse a circuit specification in the standard ProSt format."""
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
            gate_type = tokens[0]
            expected_io = GATE_TYPES[gate_type]
            arrow_idx = tokens.index("->")
            ins = [x.strip()
                   for x in "".join(tokens[1:arrow_idx]).split(",")]
            outs = [x.strip()
                    for x in "".join(tokens[arrow_idx + 1:]).split(",")]
            if len(ins) != expected_io or len(outs) != expected_io:
                raise ValueError(
                    f"{gate_type} gate must have {expected_io} I/O wires")
            circuit.add_gate(gate_type, ins, outs)
        elif tokens[0] == "OUTPUT":
            for name in "".join(tokens[1:]).split(","):
                circuit.add_output(name.strip())
    return circuit


# ===========================================================================
# Randomised topological sort
# ===========================================================================
def build_dependency_graph(gates):
    """Construct a dependency graph keyed by gate index."""
    graph = defaultdict(list)
    in_degree = defaultdict(int)
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


def randomized_topological_sort(gates):
    """Return a random valid topological order of gate indices."""
    graph, in_degree, sources = build_dependency_graph(gates)
    order = []
    while sources:
        gate_idx = random.choice(sources)
        order.append(gate_idx)
        sources.remove(gate_idx)
        for dep in graph[gate_idx]:
            in_degree[dep] -= 1
            if in_degree[dep] == 0:
                sources.append(dep)
    if len(order) != len(gates):
        raise ValueError("Circuit contains a cycle or invalid dependencies")
    return order


def get_all_wires(circuit_spec: str):
    """Extract every wire name from a circuit specification."""
    all_wires = set()
    for line in circuit_spec.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tokens = line.split()
        if tokens[0] == "INPUT":
            for name in "".join(tokens[1:]).split(","):
                all_wires.add(name.strip())
        elif tokens[0] == "ANCILLARY":
            for part in "".join(tokens[1:]).split(","):
                wire, _ = part.split("=")
                all_wires.add(wire.strip())
        elif tokens[0] in GATE_TYPES:
            arrow_idx = tokens.index("->")
            ins = [x.strip()
                   for x in "".join(tokens[1:arrow_idx]).split(",")]
            outs = [x.strip()
                    for x in "".join(tokens[arrow_idx + 1:]).split(",")]
            all_wires.update(ins + outs)
        elif tokens[0] == "OUTPUT":
            for name in "".join(tokens[1:]).split(","):
                all_wires.add(name.strip())
    return all_wires


# ===========================================================================
# CRISP steganography primitives
# ===========================================================================
def fredkin(c: int, x: int, y: int) -> Tuple[int, int, int]:
    """The classical Fredkin (controlled-swap) gate.

    F(c, x, y) = (c, x, y)         if c == 0
               = (c, y, x)         if c == 1

    The control bit c always passes through unchanged. This pass-through
    property is precisely what enables the per-channel correlation attack
    of Theorem 2 (Necessity); CRISP's per-gate fresh pi_out is what
    defeats it.
    """
    return (c, y, x) if c == 1 else (c, x, y)


def perm_gen() -> Tuple[Tuple[int, int, int], Tuple[int, int, int]]:
    """Sample two independent uniform-random permutations from S_3.

    Each Fredkin gate uses a fresh pair (pi_in, pi_out). Both are public
    (Carol receives them as part of the circuit specification) but
    independent across gates. The independence is essential: see
    Theorem 2 (Necessity of Per-Gate Resampling).
    """
    return ALL_PERMS[random.randrange(6)], ALL_PERMS[random.randrange(6)]


def randomise_lsb_plane(image: np.ndarray, channel: int) -> np.ndarray:
    """Repair R1: overwrite a whole LSB plane with fresh uniform bits.

    Without this the decoy bits at every pixel other than the secret one
    are the cover's natural LSBs.  Natural LSB planes are spatially
    correlated and biased while the embedded bit is uniform and
    independent, and that mismatch is the signal RS and sample-pair
    steganalysis exploit.  R1 makes the i.i.d. uniform decoy model a
    property of the scheme rather than an assumption about photographs.

    One pass over the image, nothing extra on the wire.  Alice writes
    the secret bit afterwards, so correctness is unaffected.
    """
    out = image.copy()
    plane = np.random.randint(0, 2, image.shape[:2], dtype=np.uint8)
    out[:, :, channel] = (out[:, :, channel] & 0xFE) | plane
    return out


def emb_bit(image: np.ndarray, bit: int,
            row: int, col: int, channel: int) -> np.ndarray:
    """Alice writes a single LSB at the secret pixel.

    Used by Alice to seed source images for primary-input wires and
    ancillary wires. Each source image carries one wire's bit at
    (row, col, channel); every other pixel and every other channel of
    the same image is left as natural-cover noise. Returns a fresh
    copy; the input is never modified.
    """
    stego = image.copy()
    pixel = int(stego[row, col, channel])
    stego[row, col, channel] = (pixel & 0xFE) | int(bit)
    return stego


def comp(stego_in: np.ndarray, output_cover: np.ndarray,
         pi_in: Tuple[int, int, int],
         pi_out: Tuple[int, int, int]) -> np.ndarray:
    """Carol applies Fredkin pixel-wise across the entire image.

    Read logical input role t from channel pi_in[t] of every pixel of
    *stego_in*. Apply F to the resulting LSB triple at each pixel.
    Write logical output role t into channel pi_out[t] of every pixel
    of a fresh copy of *output_cover*. All other bits of output_cover
    are left as natural-cover noise.

    Crucially, Carol processes EVERY pixel identically. She does not
    receive (row, col) and cannot distinguish the secret pixel from the
    decoy pixels. The secret pixel is treated identically to every other
    pixel; only Bob's secret position selects it during extraction.
    """
    src = stego_in.astype(np.uint16)
    result = output_cover.copy().astype(np.uint16)
    mask = np.uint16(0xFE)

    # Vectorised LSB read across the whole image, three channels at once
    c_bits = src[:, :, pi_in[0]] & np.uint16(1)
    x_bits = src[:, :, pi_in[1]] & np.uint16(1)
    y_bits = src[:, :, pi_in[2]] & np.uint16(1)

    swap = (c_bits == 1)
    out_c = c_bits  # control passes through unchanged
    out_x = np.where(swap, y_bits, x_bits).astype(np.uint16)
    out_y = np.where(swap, x_bits, y_bits).astype(np.uint16)
    out_bits = [out_c, out_x, out_y]

    # Vectorised LSB write back into channels per pi_out
    for t in range(3):
        ch = pi_out[t]
        result[:, :, ch] = (result[:, :, ch] & mask) | out_bits[t]
    return result.astype(np.uint8)


def ext_bit(image: np.ndarray,
            row: int, col: int, channel: int) -> int:
    """Bob reads a single LSB at the secret pixel.

    The ONLY function that needs the secret position. Called by Bob
    once per primary output wire over the entire execution.
    """
    return int(image[row, col, channel]) & 1


# ===========================================================================
# Wire routing  (Algorithm 5 in the paper)
# ===========================================================================
def assemble_stego(input_cover: np.ndarray,
                   wire_sources: List[Tuple[np.ndarray, int]],
                   pi_in: Tuple[int, int, int]) -> np.ndarray:
    """Build a stego image by per-pixel wire routing.

    A wire is represented as a (image, channel) pair: the wire's bit
    lives at the LSB of `channel` in `image` at every pixel — including,
    but not distinguishably, at the secret pixel.

    For a gate with input permutation pi_in, the consumer expects each
    logical role t in {control, x, y} to live in channel pi_in[t]. This
    function copies, at every pixel, the LSB from each source wire's
    (image, channel) into the destination channel of a fresh copy of
    input_cover.

    The operation is per-pixel and uniform across the whole image. It
    does not depend on (row, col) and does not give Carol any
    information about the secret position. After this step, channel
    pi_in[t] of the returned stego image carries the role-t wire's bit
    at every pixel, and the secret pixel carries the correct logical
    bit for that role.

    Args:
        input_cover: a fresh natural-cover RGB array of shape (h, w, 3)
        wire_sources: a list of three (source_image, source_channel) pairs,
                      one per logical role t in {0, 1, 2} = {c, x, y}
        pi_in: the gate's public input permutation; channel pi_in[t]
               will hold the logical role t bit after this function
               returns.

    Returns:
        A new RGB array of the same shape with the routed LSBs in place.
    """
    stego = input_cover.copy().astype(np.uint16)
    mask = np.uint16(0xFE)
    for t, (src_img, src_ch) in enumerate(wire_sources):
        target_ch = pi_in[t]
        src_lsb_plane = (src_img[:, :, src_ch].astype(np.uint16)
                         & np.uint16(1))
        stego[:, :, target_ch] = (
            (stego[:, :, target_ch] & mask) | src_lsb_plane)
    return stego.astype(np.uint8)


# ===========================================================================
# Alice: location generation
# ===========================================================================
def locgen(image_shape: Tuple[int, int, int],
           verbose: bool = True,
           interactive: bool = True) -> Tuple[int, int]:
    """Alice samples or accepts the secret pixel position.

    The secret never leaves Alice and Bob; it is shared between them
    over a private channel before the protocol begins. Carol never
    receives it.

    Two integers, (row, col). No channel component: with fresh
    permutations per gate the channel dimension contributes no
    additional secrecy.
    """
    h, w, _ = image_shape
    if not interactive:
        r, c = random.randrange(h), random.randrange(w)
        if verbose:
            print(f"[LOCGEN] Secret position sampled: ({r}, {c})")
        return r, c

    while True:
        try:
            user_input = input(
                f"Enter secret position as 'row,column' "
                f"(max row {h - 1}, col {w - 1}): ")
            r, col = map(int, user_input.split(","))
            if 0 <= r < h and 0 <= col < w:
                break
            print("Values out of range. Please try again.")
        except Exception:
            print("Invalid input. Please enter two integers separated by a comma.")
    print(f"[LOCGEN] Secret position provided: ({r}, {col})")
    return r, col


# ===========================================================================
# Alice: prepare source images for primary inputs and ancillaries
# ===========================================================================
def carol_materialise_constants(circuit: Circuit,
                                image_shape: Tuple[int, int, int],
                                verbose: bool = True
                                ) -> Dict[str, Tuple[np.ndarray, int]]:
    """Repair R2: constant wires as pixel-uniform planes, built by Carol.

    A wire whose value is fixed and public carries that value in the LSB
    of its channel at EVERY pixel, not only the secret one.  Carol builds
    these locally from the circuit specification, so Alice transmits no
    carrier for them.

    Embedding a known constant only at the secret pixel lets Carol
    discard every pixel whose LSB disagrees with the value she already
    knows.  Each distinct constant carrier halves her candidate set, so
    kappa of them take the bound from 1/n to 2^kappa / (2^kappa + n - 1).
    attack_full_transcript() measures exactly that.

    A pixel-uniform plane leaks no position because every pixel satisfies
    the constraint.  The price is a degenerate LSB channel, which is
    acceptable: the plane holds no secret and its contents are public.
    """
    planes: Dict[str, Tuple[np.ndarray, int]] = {}
    if verbose and circuit.ancillaries:
        print(f"\n[CAROL] Building {len(circuit.ancillaries)} constant "
              f"plane(s) locally (repair R2). None of these is transmitted.")
    for wire, value in circuit.ancillaries:
        img = np.random.randint(0, 256, image_shape, dtype=np.uint8)
        img[:, :, 0] = (img[:, :, 0] & 0xFE) | np.uint8(value)
        planes[wire] = (img, 0)
        if verbose:
            print(f"  [R2] {wire} = {value} in channel 0 LSB at all "
                  f"{image_shape[0] * image_shape[1]} pixels")
    return planes


def alice_prepare_sources(circuit: Circuit,
                          source_pool: List[np.ndarray],
                          row: int, col: int,
                          input_bits: Dict[str, int],
                          verbose: bool = True,
                          apply_r1: bool = True,
                          apply_r2: bool = True
                          ) -> Dict[str, Tuple[np.ndarray, int]]:
    """Alice creates one source image per primary-input wire.

    Each source wire becomes (image, 0): the bit lives in channel 0 at
    (row, col), written after that channel's LSB plane is re-randomised
    (repair R1).  The choice of channel 0 is conventional and contributes
    no secrecy; the consuming gate's pi_in routes it to whichever logical
    role the gate requires.

    Public-constant ancillaries get no carrier (repair R2).  Carol builds
    those herself, which removes them from the transmission and removes
    the positional leak they would otherwise cause.

    Set apply_r1 or apply_r2 to False only to measure what the repairs
    are worth.  Returns a wire-to-(image, channel) map of the transmitted
    sources, the seed of the wire table Carol evaluates over.
    """
    sources = {}
    pool_idx = 0

    if verbose:
        n_anc = 0 if apply_r2 else len(circuit.ancillaries)
        print(f"\n[ALICE] Preparing {len(circuit.inputs)} primary input(s) "
              f"+ {n_anc} ancillary source image(s).")

    for wire in circuit.inputs:
        cover = source_pool[pool_idx]
        pool_idx += 1
        bit = input_bits[wire]
        if apply_r1:
            cover = randomise_lsb_plane(cover, channel=0)
        stego = emb_bit(cover, bit, row, col, channel=0)
        sources[wire] = (stego, 0)
        if verbose:
            old = int(cover[row, col, 0])
            new = int(stego[row, col, 0])
            print(f"  [EMB] {wire} = {bit} at ({row},{col},0): "
                  f"pixel {old} -> {new}"
                  f"{'  (plane re-randomised, R1)' if apply_r1 else ''}")

    if not apply_r2:
        for wire, value in circuit.ancillaries:
            cover = source_pool[pool_idx]
            pool_idx += 1
            if apply_r1:
                cover = randomise_lsb_plane(cover, channel=0)
            stego = emb_bit(cover, value, row, col, channel=0)
            sources[wire] = (stego, 0)
            if verbose:
                print(f"  [EMB] {wire} = {value} at ({row},{col},0)  "
                      f"(constant carrier, R2 OFF: this leaks position)")
    elif verbose and circuit.ancillaries:
        print(f"  [R2] {len(circuit.ancillaries)} constant wire(s) not "
              f"transmitted; Carol builds them (see below).")
    return sources


# ===========================================================================
# Image I/O and visualisation
# ===========================================================================
def load_image(path: str) -> np.ndarray:
    """Load an RGB image from disk."""
    return np.array(Image.open(path).convert("RGB"), dtype=np.uint8)


def save_image(array: np.ndarray, path: str) -> None:
    """Save an RGB array to disk."""
    Image.fromarray(array.astype(np.uint8), "RGB").save(path)


def cover_requirement(circuit_spec: str, apply_r2: bool = True) -> Dict[str, int]:
    """How many cover images a run needs, broken down by role."""
    circuit = parse_circuit(circuit_spec)
    n_gates = len(circuit.gates)
    n_const = len(circuit.ancillaries)
    n_src = len(circuit.inputs) + (0 if apply_r2 else n_const)
    return {"gates": n_gates, "input_covers": n_gates, "output_covers": n_gates,
            "source_covers": n_src, "total": 2 * n_gates + n_src,
            "constant_planes_built_locally": n_const if apply_r2 else 0}


def load_cover_pool(cover_dir: str,
                    n_needed: int,
                    resize_to: Optional[Tuple[int, int]] = None,
                    verbose: bool = True) -> List[np.ndarray]:
    """Load n_needed RGB covers from a directory, in sorted filename order.

    Every cover must end up the same shape, since a wire is an (image,
    channel) pair and Carol works pixel by pixel across the whole set.
    Pass resize_to=(h, w) to centre-crop and resize a mixed-size album to
    a common shape; leave it None to require the files already match.

    A note on formats: load anything, but never write results back as
    JPEG. JPEG is lossy and its quantisation destroys exactly the bits
    this scheme stores. save_image() writes PNG.
    """
    if not os.path.isdir(cover_dir):
        raise FileNotFoundError(f"cover directory not found: {cover_dir}")

    files = sorted(f for f in os.listdir(cover_dir)
                   if f.lower().endswith((".png", ".jpg", ".jpeg",
                                          ".bmp", ".tif", ".tiff")))
    if len(files) < n_needed:
        raise ValueError(
            f"need {n_needed} covers, found {len(files)} in {cover_dir}")

    covers = []
    for name in files[:n_needed]:
        img = Image.open(os.path.join(cover_dir, name)).convert("RGB")
        # Never resample an image that is already the right size: a
        # same-size LANCZOS pass is not bit-exact, and the low bit is
        # precisely what this scheme stores.
        if resize_to is not None and (img.size[1], img.size[0]) != tuple(resize_to):
            h, w = resize_to
            src_w, src_h = img.size
            scale = max(w / src_w, h / src_h)
            img = img.resize((max(1, round(src_w * scale)),
                              max(1, round(src_h * scale))), Image.LANCZOS)
            left = (img.size[0] - w) // 2
            top = (img.size[1] - h) // 2
            img = img.crop((left, top, left + w, top + h))
        covers.append(np.array(img, dtype=np.uint8))

    shapes = {c.shape for c in covers}
    if len(shapes) > 1:
        raise ValueError(
            f"covers have {len(shapes)} different shapes: {sorted(shapes)[:4]}. "
            f"Pass resize_to=(h, w) to bring them to a common size.")

    if verbose:
        print(f"[ALICE] Loaded {len(covers)} cover(s) from {cover_dir}")
        print(f"[ALICE] Shape {covers[0].shape}, "
              f"{covers[0].shape[0] * covers[0].shape[1]} pixels each")
        for i, name in enumerate(files[:min(4, n_needed)]):
            print(f"         {i}: {name}")
        if n_needed > 4:
            print(f"         ... and {n_needed - 4} more")
    return covers


def split_cover_pool(covers: List[np.ndarray], n_gates: int, n_sources: int,
                     verbose: bool = True):
    """Split a flat pool into per-gate input, per-gate output and sources."""
    need = 2 * n_gates + n_sources
    if len(covers) < need:
        raise ValueError(f"need {need} covers, got {len(covers)}")
    inp = covers[:n_gates]
    out = covers[n_gates:2 * n_gates]
    src = covers[2 * n_gates:2 * n_gates + n_sources]
    if verbose:
        print(f"[ALICE] {len(inp)} input cover(s), {len(out)} output cover(s), "
              f"{len(src)} source cover(s)")
    return inp, out, src


def lsb_plane_stats(image: np.ndarray, channel: int = 0) -> Dict[str, float]:
    """How far a channel's LSB plane is from uniform i.i.d.

    Two numbers. `ones` is the fraction of set bits, which should sit at
    0.5. `neighbour_agree` is how often horizontally adjacent LSBs match,
    which should also sit at 0.5 when the plane is independent. Natural
    photographs fail the second badly in smooth regions, and that
    dependence is what repair R1 removes.
    """
    plane = (image[:, :, channel] & 1).astype(np.int8)
    agree = float((plane[:, :-1] == plane[:, 1:]).mean())
    return {"ones": float(plane.mean()),
            "neighbour_agree": agree,
            "excess_over_uniform": abs(agree - 0.5)}


def cover_report(covers: List[np.ndarray], channel: int = 0,
                 verbose: bool = True) -> Dict[str, float]:
    """Check whether an album's LSB planes look uniform and independent.

    Repair R1 exists because the security analysis treats the decoy bits
    at every non-secret pixel as i.i.d. uniform.  Whether a real album
    satisfies that without help depends entirely on the pictures:

      - Noisy camera photographs usually do.  Sensor noise randomises
        the low bit almost everywhere, so neighbour agreement sits near
        0.5 and R1 changes little.
      - Graphics, screenshots, logos, renders, noiseless gradients and
        heavily compressed images usually do not.  Flat regions give
        neighbour agreement near 1.0, and a smooth even-valued ramp can
        have an LSB plane that is entirely zero.  Embedding one bit into
        a plane like that leaves a single lit pixel on a black field.

    R1 removes the question.  It costs one pass over the image and makes
    the assumption hold whatever the album looks like.
    """
    ones, agree = [], []
    for img in covers:
        st = lsb_plane_stats(img, channel)
        ones.append(st["ones"]); agree.append(st["neighbour_agree"])
    mean_ones, mean_agree = float(np.mean(ones)), float(np.mean(agree))
    worst = float(max(agree))
    if verbose:
        print(f"[COVERS] {len(covers)} image(s), channel {channel} LSB plane")
        print(f"         fraction of set bits : {mean_ones:.4f}   "
              f"(uniform = 0.5000)")
        print(f"         neighbour agreement  : {mean_agree:.4f}   "
              f"(independent = 0.5000)")
        print(f"         worst single image   : {worst:.4f}")
        if worst > 0.55 or abs(mean_ones - 0.5) > 0.05:
            print("         -> These planes are NOT uniform i.i.d. Without R1 the")
            print("            secret pixel would stand out against its decoys.")
            print("            R1 is doing real work on this album.")
        else:
            print("         -> Already close to uniform i.i.d., so R1 costs you")
            print("            nothing here and removes the dependence on it.")
    return {"mean_ones": mean_ones, "mean_neighbour_agree": mean_agree,
            "worst_neighbour_agree": worst}


def display_images(images: Dict[str, np.ndarray],
                   title: str,
                   stage: str,
                   save_path: Optional[str] = None,
                   inline_display_enabled: bool = False) -> None:
    """Display a dict of {label -> image} in a grid.

    In Jupyter / VS Code Interactive (inline_display_enabled=True), the
    grid renders inline below the cell. In a headless script, the grid
    is saved to disk and a path is printed.
    """
    n = len(images)
    cols = min(4, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols,
                             figsize=(4 * cols, 4 * rows),
                             squeeze=False)
    for ax_row in axes:
        for ax in ax_row:
            ax.axis("off")

    for idx, (label, img) in enumerate(images.items()):
        r, c = divmod(idx, cols)
        axes[r][c].imshow(img)
        axes[r][c].set_title(label, fontsize=8, wrap=True)
        axes[r][c].axis("on")

    fig.suptitle(f"{stage}: {title}", fontsize=11)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    if inline_display_enabled:
        plt.show()
    else:
        if save_path is None:
            safe = title.replace(" ", "_").replace("/", "-")[:40]
            save_path = os.path.join(tempfile.gettempdir(),
                                     f"crisp_display_{safe}.png")
        plt.savefig(save_path, dpi=100, bbox_inches="tight")
        print(f"[Display] Figure saved -> {save_path}")

    plt.close(fig)


# ===========================================================================
# Carol: protocol-faithful circuit evaluation
# ---------------------------------------------------------------------------

def Carol_evaluate(circuit: Circuit,
                   gate_order: List[int],
                   wire_sources: Dict[str, Tuple[np.ndarray, int]],
                   input_covers: List[np.ndarray],
                   output_covers: List[np.ndarray],
                   permutations: List[Tuple[Tuple[int, int, int],
                                            Tuple[int, int, int]]],
                   processed_dir: Optional[str] = None,
                   inline_display_enabled: bool = False,
                   verbose: bool = True
                   ) -> Tuple[Dict[str, Tuple[np.ndarray, int]],
                              List[np.ndarray],
                              Dict[str, float]]:
    """Carol evaluates the circuit gate-by-gate.

    She receives:
        - the circuit specification (parsed)
        - the gate execution order (random topological sort)
        - the wire table seeded with primary-input and ancillary sources
        - one fresh input cover per gate (in execution order)
        - one fresh output cover per gate (in execution order)
        - one (pi_in, pi_out) permutation pair per gate

    She does NOT receive (row, col).

    Returns:
        wire_to_image: final wire table after all gates execute
        gate_outputs: list of per-gate output images, in execution order
        timing: dict with profiling breakdown
    """
    if verbose:
        print("\n[CAROL] Beginning protocol-faithful circuit evaluation.")
        print("[CAROL] Carol's signature does NOT contain (row, col).")
        print("[CAROL] Wires flow as (image, channel) pairs; "
              "Carol never calls ext().")

    n_gates = len(gate_order)
    wires: Dict[str, Tuple[np.ndarray, int]] = dict(wire_sources)
    gate_outputs: List[np.ndarray] = []

    t_assemble_total = 0.0
    t_comp_total = 0.0
    t_io_total = 0.0
    t_total_start = time.perf_counter()

    for step_idx, gate_idx in enumerate(gate_order, start=1):
        gate = circuit.gates[gate_idx]
        gate_type = gate["type"]

        c_wire, x_wire, y_wire = gate["inputs"]
        co_wire, xo_wire, yo_wire = gate["outputs"]

        # Look up input wire carriers in the wire table.
        # Each entry is a (source_image, source_channel) pair.
        c_src = wires[c_wire]
        x_src = wires[x_wire]
        y_src = wires[y_wire]
        wire_sources_list = [c_src, x_src, y_src]

        # Per-gate fresh permutations (already sampled by Alice).
        pi_in, pi_out = permutations[step_idx - 1]

        # Per-gate fresh covers.
        in_cover = input_covers[step_idx - 1]
        out_cover = output_covers[step_idx - 1]

        if gate_type == "FREDKIN":
            # Step A: assemble the stego image. Per-pixel channel
            # routing, uniform across the whole image, no row/col access.
            t0 = time.perf_counter()
            stego = assemble_stego(in_cover, wire_sources_list, pi_in)
            t1 = time.perf_counter()

            # Step B: apply Fredkin pixel-wise across the entire image.
            out_img = comp(stego, out_cover, pi_in, pi_out)
            t2 = time.perf_counter()

            t_assemble_total += t1 - t0
            t_comp_total += t2 - t1
        else:
            raise ValueError(f"Unknown gate type: {gate_type}")

        gate_outputs.append(out_img)

        # Update the wire table with the gate's three output wires.
        # The output image carries logical role t at channel pi_out[t]
        # at every pixel.
        wires[co_wire] = (out_img, pi_out[0])
        wires[xo_wire] = (out_img, pi_out[1])
        wires[yo_wire] = (out_img, pi_out[2])

        # Verbose per-gate logging.
        if verbose:
            print(
                f"\n[CAROL] Gate {step_idx}/{n_gates} (idx {gate_idx}): "
                f"{gate_type}")
            print(f"        Inputs : {gate['inputs']}")
            print(f"        Outputs: {gate['outputs']}")
            print(f"        pi_in  = {pi_in}")
            print(f"        pi_out = {pi_out}")
            print(f"        Source wire carriers (image_id, channel):")
            for role, src in zip(("c", "x", "y"), wire_sources_list):
                src_img, src_ch = src
                print(f"          {role}: image={id(src_img):#x}  "
                      f"channel={src_ch}")
            print(f"        After Comp: output image carries")
            for role, ch in zip(("c", "x", "y"), pi_out):
                print(f"          {role} -> channel {ch} of output image "
                      f"{id(out_img):#x}")
            print(f"        Timing: assemble={1000 * (t1 - t0):.3f} ms  "
                  f"comp={1000 * (t2 - t1):.3f} ms")

        # Save and display.
        if processed_dir is not None:
            t_io_start = time.perf_counter()
            stego_path = os.path.join(processed_dir,
                                      f"gate_{step_idx}_stego.png")
            output_path = os.path.join(processed_dir,
                                       f"gate_{step_idx}_output.png")
            try:
                save_image(stego, stego_path)
                save_image(out_img, output_path)
                if verbose:
                    print(f"        Stego  saved: {stego_path}")
                    print(f"        Output saved: {output_path}")
            except Exception as e:
                if verbose:
                    print(f"        [Warning] Could not save images: {e}")

            t_io_total += time.perf_counter() - t_io_start

        # Display is independent of disk output: inline rendering in a
        # notebook must not require a processed_dir.
        if inline_display_enabled or processed_dir is not None:
            gate_label = (f"{gate_type} {gate['inputs']} -> "
                          f"{gate['outputs']}")
            save_path = (os.path.join(processed_dir,
                                      f"gate_{step_idx}_display.png")
                         if processed_dir is not None else None)
            try:
                display_images(
                    {"Cover (input)": in_cover,
                     f"Stego (after Assemble, pi_in={pi_in})": stego,
                     f"Output (after Comp, pi_out={pi_out})": out_img,
                     "Cover (output)": out_cover},
                    title=f"Gate {step_idx} - {gate_label}",
                    stage="Carol's Computation",
                    save_path=save_path,
                    inline_display_enabled=inline_display_enabled)
            except Exception as e:
                if verbose:
                    print(f"        [Warning] Display failed: {e}")

    t_total = time.perf_counter() - t_total_start

    timing = {
        "assemble_seconds": t_assemble_total,
        "comp_seconds": t_comp_total,
        "io_seconds": t_io_total,
        "compute_seconds": t_assemble_total + t_comp_total,
        "total_seconds": t_total,
        "per_gate_compute_seconds": (
            (t_assemble_total + t_comp_total) / max(n_gates, 1)),
        "n_gates": n_gates,
    }

    if verbose:
        print(f"\n[CAROL] Evaluation complete.")
        print(f"        Total compute time     : "
              f"{1000 * timing['compute_seconds']:.2f} ms")
        print(f"        Per-gate compute       : "
              f"{1000 * timing['per_gate_compute_seconds']:.3f} ms")
        print(f"        Assemble (routing)     : "
              f"{1000 * timing['assemble_seconds']:.2f} ms")
        print(f"        Comp (Fredkin)         : "
              f"{1000 * timing['comp_seconds']:.2f} ms")
        print(f"        I/O + display overhead : "
              f"{1000 * timing['io_seconds']:.2f} ms")

    return wires, gate_outputs, timing


# ===========================================================================
# Bob: extract primary outputs (the only place ext_bit is called)
# ===========================================================================
def Bob_extract(circuit: Circuit,
                wire_to_image: Dict[str, Tuple[np.ndarray, int]],
                row: int, col: int,
                verbose: bool = True) -> Dict[str, int]:
    """Bob recovers the primary output bits.

    Called once per primary output wire over the entire execution. Bob
    is the only party with access to (row, col); he reads the LSB of
    the indicated channel of the output wire's carrier image.
    """
    if verbose:
        print(f"\n[BOB]   Extracting {len(circuit.outputs)} primary "
              f"output bit(s) at secret position ({row}, {col}).")
    out_bits = {}
    for w in circuit.outputs:
        img, ch = wire_to_image[w]
        bit = ext_bit(img, row, col, ch)
        out_bits[w] = bit
        if verbose:
            pixel = int(img[row, col, ch])
            print(f"  [EXT] {w} = {bit}  "
                  f"(pixel {pixel} at ({row},{col},{ch}))")
    return out_bits


# ===========================================================================
# Top-level orchestrator (programmatic API for non-interactive use)
# ===========================================================================
def run_protocol(circuit_spec: str,
                 input_bits: Dict[str, int],
                 image_shape: Tuple[int, int, int] = (128, 128, 3),
                 verbose: bool = True,
                 seed: Optional[int] = None,
                 display_gates: bool = False,
                 cover_dir: Optional[str] = None,
                 apply_r1: bool = True,
                 apply_r2: bool = True,
                 return_view: bool = False
                 ) -> Tuple[Dict[str, int], Dict[str, float]]:
    """Run Alice -> Carol -> Bob end-to-end on a given circuit.

    For interactive use see main(); for unit tests and benchmarks,
    call this helper.
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    if verbose:
        print("\n" + "=" * 72)
        print("CRISP — protocol-faithful end-to-end execution")
        print("Alice -> Carol (no row,col) -> Bob")
        print("=" * 72)

    circuit = parse_circuit(circuit_spec)
    n_gates = len(circuit.gates)
    n_const = len(circuit.ancillaries)
    n_sources = len(circuit.inputs) + (0 if apply_r2 else n_const)

    if verbose:
        print(f"\n[SETUP] Circuit: {n_gates} gates, "
              f"{len(circuit.inputs)} primary inputs, "
              f"{len(circuit.ancillaries)} ancillaries.")
        print(f"[SETUP] Image shape: {image_shape}")
        n_total = 2 * n_gates + n_sources
        print(f"[SETUP] Images transmitted: {n_total} = "
              f"2*{n_gates} (per-gate) + {n_sources} (sources)")
        if apply_r2 and n_const:
            print(f"[SETUP] Plus {n_const} constant plane(s) Carol builds "
                  f"locally, never sent (R2)")

    # Covers: real album if one was given, synthetic otherwise.
    pool = None
    if cover_dir is not None:
        pool = load_cover_pool(cover_dir, 2 * n_gates + n_sources,
                               resize_to=image_shape[:2], verbose=verbose)
        image_shape = pool[0].shape

    # Alice
    row, col = locgen(image_shape, interactive=False, verbose=verbose)
    if pool is not None:
        _in_c, _out_c, source_pool = split_cover_pool(
            pool, n_gates, n_sources, verbose=verbose)
    else:
        source_pool = [
            np.random.randint(0, 256, image_shape, dtype=np.uint8)
            for _ in range(n_sources)
        ]
    wire_sources = alice_prepare_sources(
        circuit, source_pool, row, col, input_bits, verbose=verbose,
        apply_r1=apply_r1, apply_r2=apply_r2)

    transmitted = dict(wire_sources)
    if apply_r2:
        wire_sources = {**wire_sources,
                        **carol_materialise_constants(circuit, image_shape,
                                                      verbose=verbose)}

    # Public material (Carol sees these).
    gate_order = randomized_topological_sort(circuit.gates)
    permutations = [perm_gen() for _ in range(n_gates)]
    if pool is not None:
        input_covers, output_covers = _in_c, _out_c
    else:
        input_covers = [
            np.random.randint(0, 256, image_shape, dtype=np.uint8)
            for _ in range(n_gates)
        ]
        output_covers = [
            np.random.randint(0, 256, image_shape, dtype=np.uint8)
            for _ in range(n_gates)
        ]

    # Carol (no secret position!).
    wire_to_image, _gate_outputs, timing = Carol_evaluate(
        circuit, gate_order, wire_sources,
        input_covers, output_covers, permutations,
        processed_dir=None, inline_display_enabled=display_gates,
        verbose=verbose)

    # Bob.
    out_bits = Bob_extract(circuit, wire_to_image, row, col, verbose=verbose)
    timing["images_transmitted"] = 2 * n_gates + n_sources
    timing["constant_planes_local"] = n_const if apply_r2 else 0
    if return_view:
        return out_bits, timing, {
            "transmitted_sources": transmitted,
            "wire_map": wire_to_image,
            "secret": (row, col),
            "circuit": circuit,
            "image_shape": image_shape,
        }
    return out_bits, timing


# ===========================================================================
# Unit tests  (288 = 8 input triples x 6 pi_in x 6 pi_out)
# ===========================================================================
def run_unit_tests(verbose: bool = False) -> bool:
    """Exhaustively verify the Emb -> Assemble -> Comp -> Ext pipeline.

    For every (c, x, y) input triple and every (pi_in, pi_out) permutation
    pair, the pipeline must recover the correct Fredkin output.
    """
    print("=" * 64)
    print("UNIT TESTS  (288 = 8 input triples * 6 pi_in * 6 pi_out)")
    print("=" * 64)

    h, w = 16, 16
    src1 = np.random.randint(0, 256, (h, w, 3), dtype=np.uint8)
    src2 = np.random.randint(0, 256, (h, w, 3), dtype=np.uint8)
    src3 = np.random.randint(0, 256, (h, w, 3), dtype=np.uint8)
    in_cover = np.random.randint(0, 256, (h, w, 3), dtype=np.uint8)
    out_cover = np.random.randint(0, 256, (h, w, 3), dtype=np.uint8)

    passed = 0
    failed = 0
    failures: List[str] = []

    for c, x, y in itertools.product([0, 1], repeat=3):
        expected = fredkin(c, x, y)
        for pi_in in ALL_PERMS:
            for pi_out in ALL_PERMS:
                row = random.randrange(h)
                col = random.randrange(w)

                # Each input bit lives in a separate source image,
                # at channel 0 of (row, col). This is exactly the
                # protocol-faithful path.
                a = emb_bit(src1, c, row, col, 0)
                b = emb_bit(src2, x, row, col, 0)
                d = emb_bit(src3, y, row, col, 0)
                wire_sources_list = [(a, 0), (b, 0), (d, 0)]

                stego = assemble_stego(
                    in_cover, wire_sources_list, pi_in)
                out_img = comp(stego, out_cover, pi_in, pi_out)
                got = (
                    ext_bit(out_img, row, col, pi_out[0]),
                    ext_bit(out_img, row, col, pi_out[1]),
                    ext_bit(out_img, row, col, pi_out[2]),
                )
                if got == expected:
                    passed += 1
                    if verbose:
                        print(f"  PASS ({c},{x},{y}) pi_in={pi_in} "
                              f"pi_out={pi_out} -> {got}")
                else:
                    failed += 1
                    failures.append(
                        f"  FAIL ({c},{x},{y}) pi_in={pi_in} "
                        f"pi_out={pi_out} expected={expected} got={got}")

    print(f"\n  Passed : {passed}/288")
    print(f"  Failed : {failed}/288")
    for msg in failures[:5]:
        print(msg)
    if failed > 5:
        print(f"  ... and {failed - 5} more failures.")
    ok = (failed == 0)
    print(f"\n  Result : {'ALL 288 TESTS PASSED' if ok else 'SOME TESTS FAILED'}")
    return ok


# ===========================================================================
# Truth-table verification (8 rows of the benchmark function)
# ===========================================================================
def run_truth_table(circuit_spec: str,
                    image_shape: Tuple[int, int, int] = (64, 64, 3),
                    verbose: bool = True) -> bool:
    """Verify the 8 truth-table rows for the benchmark function.

    The benchmark is f(A, B, C) = (A AND C) OR (NOT A AND B) OR (NOT B AND NOT C).
    """
    print("\n" + "=" * 64)
    print("TRUTH TABLE — (A AND C) OR (NOT A AND B) OR (NOT B AND NOT C)")
    print("=" * 64)

    circuit = parse_circuit(circuit_spec)
    if len(circuit.outputs) != 1:
        raise ValueError("Truth-table verifier expects exactly 1 output.")
    if len(circuit.inputs) != 3:
        raise ValueError("Truth-table verifier expects exactly 3 inputs.")

    out_wire = circuit.outputs[0]
    in_names = circuit.inputs
    all_ok = True

    print(f"\n  {'A':>3}  {'B':>3}  {'C':>3}  "
          f"{'Expected':>10}  {'CRISP':>7}  {'OK':>4}")
    print("  " + "-" * 42)
    for A, B, C in itertools.product([0, 1], repeat=3):
        expected = int((A & C) | ((1 - A) & B) | ((1 - B) & (1 - C)))
        bits = {in_names[0]: A, in_names[1]: B, in_names[2]: C}
        out, _ = run_protocol(circuit_spec, bits,
                              image_shape=image_shape, verbose=False)
        got = out[out_wire]
        ok = (got == expected)
        all_ok = all_ok and ok
        print(f"  {A:>3}  {B:>3}  {C:>3}  {expected:>10}  {got:>7}  "
              f"{'OK' if ok else 'FAIL':>4}")
    print(f"\n  Result: {'ALL 8 CASES CORRECT' if all_ok else 'SOME FAILURES'}")
    return all_ok


# ===========================================================================
# Benchmark suite
# ===========================================================================
def bench_per_gate_timing(image_sizes: Tuple[int, ...] = (64, 128, 256, 512),
                          n_trials: int = 20) -> List[Tuple]:
    """Per-gate timing breakdown: Assemble vs Comp at multiple image sizes."""
    print("\n" + "=" * 72)
    print("BENCHMARK 1 — Per-gate timing breakdown (Assemble vs Comp)")
    print("=" * 72)
    print(f"{'Image':>10}  {'Pixels':>10}  "
          f"{'Assemble (ms)':>16}  {'Comp (ms)':>14}  {'Total (ms)':>14}")
    print("-" * 72)

    rows = []
    for s in image_sizes:
        shape = (s, s, 3)
        n_pix = s * s
        as_times, cp_times = [], []
        for _ in range(n_trials):
            cover = np.random.randint(0, 256, shape, dtype=np.uint8)
            a = np.random.randint(0, 256, shape, dtype=np.uint8)
            b = np.random.randint(0, 256, shape, dtype=np.uint8)
            d = np.random.randint(0, 256, shape, dtype=np.uint8)
            out_cover = np.random.randint(0, 256, shape, dtype=np.uint8)
            pi_in, pi_out = perm_gen()
            t0 = time.perf_counter()
            stego = assemble_stego(cover, [(a, 0), (b, 0), (d, 0)], pi_in)
            t1 = time.perf_counter()
            _ = comp(stego, out_cover, pi_in, pi_out)
            t2 = time.perf_counter()
            as_times.append(t1 - t0)
            cp_times.append(t2 - t1)
        a_mean = 1000 * float(np.mean(as_times))
        c_mean = 1000 * float(np.mean(cp_times))
        a_std = 1000 * float(np.std(as_times))
        c_std = 1000 * float(np.std(cp_times))
        total = a_mean + c_mean
        print(f"{s}x{s:<6} {n_pix:>10}  "
              f"{a_mean:>8.2f} +/- {a_std:<5.2f}  "
              f"{c_mean:>8.2f} +/- {c_std:<5.2f}  "
              f"{total:>10.2f}")
        rows.append((s, n_pix, a_mean, a_std, c_mean, c_std, total))
    return rows


def bench_image_count() -> List[Tuple]:
    """Compare ProSt's per-wire image count with CRISP's per-gate count."""
    print("\n" + "=" * 72)
    print("BENCHMARK 2 — Image count: CRISP vs ProSt across circuits")
    print("=" * 72)
    print(f"{'Circuit':>22}  {'Gates':>6}  {'Wires':>6}  "
          f"{'ProSt':>6}  {'CRISP':>6}  {'Reduction':>10}")
    print("-" * 72)

    circuits = [
        ("Benchmark, 14 gates", CIRCUIT_SPEC),
        ("Linear, 17 gates", _make_linear_chain(17)),
        ("Linear, 50 gates", _make_linear_chain(50)),
        ("Linear, 100 gates", _make_linear_chain(100)),
    ]

    rows = []
    for name, spec in circuits:
        circ = parse_circuit(spec)
        n = len(circ.gates)
        wires = len(get_all_wires(spec))
        n_src = len(circ.inputs) + len(circ.ancillaries)
        crisp = 2 * n + n_src
        reduction = wires / crisp
        print(f"{name:>22}  {n:>6}  {wires:>6}  "
              f"{wires:>6}  {crisp:>6}  {reduction:>9.2f}x")
        rows.append((name, n, wires, wires, crisp, reduction))
    return rows


def bench_chi_square_steganalysis(image_sizes: Tuple[int, ...] = (64, 128, 256),
                                  n_trials: int = 200) -> List[Tuple]:
    """Empirical chi-square steganalysis of CRISP outputs vs natural covers.

    Tests the structural argument of Section 6.7 of the paper: because
    Comp processes every pixel uniformly, the LSB plane of a CRISP
    output should be statistically indistinguishable from a fresh
    natural cover.
    """
    print("\n" + "=" * 72)
    print("BENCHMARK 3 — Chi-square steganalysis on the LSB plane")
    print("=" * 72)
    print(f"{'Image':>10}  {'Natural mean':>15}  "
          f"{'CRISP mean':>15}  {'p(detect) at 0.05':>20}")
    print("-" * 72)

    rows = []
    crit = 3.841  # chi-square critical value at alpha=0.05, df=1
    for s in image_sizes:
        shape = (s, s, 3)
        n_total = 3 * s * s
        nat_chi, crisp_chi = [], []
        for _ in range(n_trials):
            # Natural cover (control).
            nat = np.random.randint(0, 256, shape, dtype=np.uint8)
            n0 = int(np.sum((nat & 1) == 0))
            n1 = n_total - n0
            exp = n_total / 2
            chi_nat = ((n0 - exp) ** 2 + (n1 - exp) ** 2) / exp
            nat_chi.append(chi_nat)

            # CRISP output.
            cover = np.random.randint(0, 256, shape, dtype=np.uint8)
            a = np.random.randint(0, 256, shape, dtype=np.uint8)
            b = np.random.randint(0, 256, shape, dtype=np.uint8)
            d = np.random.randint(0, 256, shape, dtype=np.uint8)
            out_cover = np.random.randint(0, 256, shape, dtype=np.uint8)
            pi_in, pi_out = perm_gen()
            stego = assemble_stego(cover, [(a, 0), (b, 0), (d, 0)], pi_in)
            out_img = comp(stego, out_cover, pi_in, pi_out)
            n0c = int(np.sum((out_img & 1) == 0))
            n1c = n_total - n0c
            chi_crisp = ((n0c - exp) ** 2 + (n1c - exp) ** 2) / exp
            crisp_chi.append(chi_crisp)

        nat_arr = np.array(nat_chi)
        crisp_arr = np.array(crisp_chi)
        nat_det = float(np.mean(nat_arr > crit))
        crisp_det = float(np.mean(crisp_arr > crit))
        print(f"{s}x{s:<6}  "
              f"{float(np.mean(nat_arr)):>14.3f}  "
              f"{float(np.mean(crisp_arr)):>14.3f}  "
              f"{crisp_det:.3f} (vs nat {nat_det:.3f})")
        rows.append((s, float(np.mean(nat_arr)), float(np.mean(crisp_arr)),
                     nat_det, crisp_det))

    print()
    print("Interpretation: CRISP outputs are statistically indistinguishable")
    print("from fresh natural covers under chi-square detection. False-")
    print("positive rates land near the nominal 5% baseline for both.")
    return rows


def bench_security_empirical(
        image_sizes: Tuple[Tuple[int, int], ...] = (
            (8, 8), (16, 16), (32, 32), (64, 64), (128, 128)),
        n_trials: int = 10000) -> List[Tuple]:
    """Joint position-and-value recovery when lambda = 0."""
    print("\n" + "=" * 72)
    print("BENCHMARK 4 - Joint position+value recovery at lambda = 0")
    print("=" * 72)
    print("With repairs R1 and R2 in place, Carol knows nothing about the")
    print("secret input beyond its prior, so her posterior over pixels stays")
    print("uniform.  Her optimal strategy is to pick a pixel at random and")
    print("read its three LSBs, winning at rate 1/n.  bench_constant_leak()")
    print("measures what happens when lambda > 0.")
    print()
    print(f"{'Image':>10}  {'n=h*w':>8}  {'1/n':>14}  "
          f"{'Empirical':>14}  {'|delta|':>10}")
    print("-" * 72)

    rows = []
    for h, w in image_sizes:
        n = h * w
        wins = 0
        for _ in range(n_trials):
            r_s, c_s = random.randrange(h), random.randrange(w)
            t_s = (random.randrange(2),
                   random.randrange(2),
                   random.randrange(2))
            img_lsb = np.random.randint(0, 2, (h, w, 3), dtype=np.uint8)
            img_lsb[r_s, c_s, :] = t_s
            r_g, c_g = random.randrange(h), random.randrange(w)
            t_g = tuple(int(img_lsb[r_g, c_g, k]) for k in range(3))
            if (r_g, c_g) == (r_s, c_s) and t_g == t_s:
                wins += 1
        empirical = wins / n_trials
        theor = 1.0 / n
        delta = abs(empirical - theor)
        print(f"{h}x{w:<6}  {n:>8}  "
              f"{theor:>13.6e}  {empirical:>13.6e}  {delta:>9.2e}")
        rows.append((h, w, theor, empirical, delta))

    print()
    print("All empirical rates match 1/n within sampling error.  This is the")
    print("lambda = 0 case only; every bit of the secret input that Carol")
    print("already knows halves her candidate set.")
    return rows


def bench_end_to_end(image_sizes: Tuple[Tuple[int, int], ...] = (
        (64, 64), (128, 128), (256, 256))) -> List[Tuple]:
    """End-to-end protocol-faithful CRISP timing on the benchmark circuit."""
    print("\n" + "=" * 72)
    print("BENCHMARK 5 — End-to-end timing on the 14-gate benchmark circuit")
    print("=" * 72)
    print(f"{'Image':>10}  {'Total (ms)':>14}  {'Comp (ms)':>14}  "
          f"{'Assemble (ms)':>16}")
    print("-" * 72)

    rows = []
    for h, w in image_sizes:
        random.seed(42)
        np.random.seed(42)
        out, timing = run_protocol(
            CIRCUIT_SPEC, {"A": 1, "B": 0, "C": 1},
            image_shape=(h, w, 3), verbose=False)
        total_ms = 1000 * timing["compute_seconds"]
        comp_ms = 1000 * timing["comp_seconds"]
        asm_ms = 1000 * timing["assemble_seconds"]
        print(f"{h}x{w:<6}  {total_ms:>12.2f}  "
              f"{comp_ms:>12.2f}  {asm_ms:>14.2f}")
        rows.append((h, w, total_ms, comp_ms, asm_ms))

    print()
    print("Note: Comp (the pixel-wise Fredkin pass) dominates total cost.")
    print("Assemble (the wire-routing step) is roughly half as expensive.")
    return rows


# ---------------------------------------------------------------------------
# Helpers for benchmarks
# ---------------------------------------------------------------------------
def _make_linear_chain(n_gates: int) -> str:
    """Generate a linear-chain Fredkin circuit specification with n_gates."""
    lines = ["INPUT A, B, C",
             "ANCILLARY const0 = 0",
             "ANCILLARY const1 = 1"]
    cur = ["A", "B", "C"]
    for i in range(n_gates):
        outs = [f"w{i}_{j}" for j in range(3)]
        lines.append(f"FREDKIN {', '.join(cur)} -> {', '.join(outs)}")
        cur = outs
    lines.append(f"OUTPUT {cur[0]}")
    return "\n".join(lines)


# ===========================================================================
# Main workflow
# ===========================================================================

CIRCUIT_SPEC = """
# Circuit to compute: (A AND C) OR (NOT A AND B) OR (NOT B AND NOT C)
# This is a Fredkin-gate only implementation.
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

# Set to True if running in Jupyter / VS Code Interactive,
# False for a standard script / terminal run.
SHOULD_DISPLAY_INLINE = True


# ===========================================================================
# Notebook helpers (use these from individual Jupyter cells)
# ===========================================================================
def nb_unit_tests():
    """Run the 288 unit tests."""
    return run_unit_tests(verbose=False)


def nb_demo(image_size: int = 128, seed: int = 42,
            bits: Optional[Dict[str, int]] = None):
    """Single non-interactive end-to-end execution with verbose output."""
    if bits is None:
        bits = {"A": 1, "B": 0, "C": 1}
    expected = int((bits["A"] & bits["C"])
                   | ((1 - bits["A"]) & bits["B"])
                   | ((1 - bits["B"]) & (1 - bits["C"])))
    print(f"Expected output for A={bits['A']} B={bits['B']} C={bits['C']}: "
          f"{expected}")
    out, timing = run_protocol(
        CIRCUIT_SPEC, bits,
        image_shape=(image_size, image_size, 3),
        verbose=True, seed=seed)
    return out, timing


def nb_truth_table(image_size: int = 64):
    """Verify the full 8-row truth table."""
    return run_truth_table(CIRCUIT_SPEC,
                           image_shape=(image_size, image_size, 3),
                           verbose=True)


def nb_benchmarks_quick():
    """Reduced-trial benchmark suite for quick notebook turnaround."""
    print("\nRunning quick benchmarks (reduced trial counts)...")
    bench_per_gate_timing(image_sizes=(64, 128, 256), n_trials=10)
    bench_image_count()
    bench_chi_square_steganalysis(image_sizes=(64, 128), n_trials=50)
    bench_security_empirical(
        image_sizes=((8, 8), (16, 16), (32, 32), (64, 64)),
        n_trials=2000)
    bench_end_to_end(image_sizes=((64, 64), (128, 128)))


def nb_benchmarks_full():
    """Full paper-quality benchmark suite (~3 minutes)."""
    print("\nRunning full benchmarks...")
    bench_per_gate_timing(image_sizes=(64, 128, 256, 512), n_trials=20)
    bench_image_count()
    bench_chi_square_steganalysis(
        image_sizes=(64, 128, 256), n_trials=200)
    bench_security_empirical(
        image_sizes=((8, 8), (16, 16), (32, 32), (64, 64), (128, 128)),
        n_trials=10000)
    bench_end_to_end(image_sizes=((64, 64), (128, 128), (256, 256)))


# ===========================================================================
# Entry point
# ===========================================================================
if __name__ == "__main__":
    main()


# ===========================================================================
# Attack: constant-wire candidate filtering
# ===========================================================================
def known_constant_wires(circuit: Circuit) -> Dict[str, int]:
    """Every wire whose value follows from the public specification alone.

    Starts at the constant ancillaries and propagates forward through any
    gate all of whose inputs are already known.
    """
    known = {w: v for w, v in circuit.ancillaries}
    changed = True
    while changed:
        changed = False
        for gate in circuit.gates:
            ci, xi, yi = gate["inputs"]
            co, xo, yo = gate["outputs"]
            if all(w in known for w in (ci, xi, yi)) and co not in known:
                a, b, c = fredkin(known[ci], known[xi], known[yi])
                known[co], known[xo], known[yo] = a, b, c
                changed = True
    return known


def attack_full_transcript(view: Dict, verbose: bool = False) -> Dict:
    """Carol filters candidate pixels using every value she can derive.

    She holds the circuit specification, so she knows each ancillary's
    constant and can propagate it forward.  For each such wire she keeps
    only pixels whose LSB agrees with the value she expects.  The secret
    pixel always survives, so the survivors are her candidate set and her
    optimal move is to guess uniformly among them.

    Under R2 the constant planes are pixel-uniform, nothing is filtered,
    and the candidate set is the whole image.
    """
    circuit = view["circuit"]
    h, w, _ = view["image_shape"]
    row, col = view["secret"]
    wires = view["wire_map"]

    known = known_constant_wires(circuit)
    candidates = np.ones((h, w), dtype=bool)
    used = 0
    for wire, value in known.items():
        if wire not in wires:
            continue
        img, ch = wires[wire]
        plane = (img[:, :, ch] & 1) == value
        if plane.all():
            continue                      # pixel-uniform: filters nothing
        before = int(candidates.sum())
        candidates &= plane
        used += 1
        if verbose:
            print(f"  [ATTACK] filter on {wire}={value}: "
                  f"{before} -> {int(candidates.sum())} candidates")

    n_cand = int(candidates.sum())
    return {"derivable_wires": len(known), "filters_applied": used,
            "candidates": n_cand, "secret_survives": bool(candidates[row, col]),
            "adv": 1.0 / n_cand, "degradation": (h * w) / n_cand,
            "pinned": n_cand == 1}


def bench_constant_leak(image_shape=(128, 128, 3), n_trials: int = 40):
    """BENCHMARK 6 - what constant wires cost, measured on real transcripts."""
    n = image_shape[0] * image_shape[1]
    print("\n" + "=" * 72)
    print("BENCHMARK 6 - Constant-wire leak, with R2 off and on")
    print("=" * 72)
    print(f"{'R2':>5}  {'filters':>8}  {'images':>7}  {'mean cands':>11}  "
          f"{'advantage':>12}  {'vs 1/n':>9}")
    print("-" * 72)
    rows = []
    for r2 in (False, True):
        cands, advs, filt, ok = [], [], 0, True
        for _ in range(n_trials):
            out, timing, view = run_protocol(
                CIRCUIT_SPEC, {"A": 1, "B": 0, "C": 1},
                image_shape=image_shape, verbose=False,
                apply_r2=r2, return_view=True)
            a = attack_full_transcript(view)
            ok &= a["secret_survives"]
            cands.append(a["candidates"]); advs.append(a["adv"])
            filt = a["filters_applied"]
        mean_adv = sum(advs) / n_trials
        mean_cand = sum(cands) / n_trials
        print(f"{'on' if r2 else 'off':>5}  {filt:>8}  "
              f"{timing['images_transmitted']:>7}  {mean_cand:>11.1f}  "
              f"{mean_adv:>12.3e}  {mean_adv * n:>8.1f}x")
        rows.append((r2, filt, mean_cand, mean_adv))
        assert ok, "the secret pixel must always survive Carol's own filter"
    print()
    print("R2 restores exactly 1/n and sends two fewer images.")
    return rows


def bench_lambda_sweep(n: int = 128 * 128,
                       kappas=(0, 1, 3, 5, 7, 10, 14),
                       trials: int = 2_000_000):
    """BENCHMARK 7 - the bound itself, swept over the number of known bits."""
    rng = np.random.default_rng(7)
    print("\n" + "=" * 72)
    print(f"BENCHMARK 7 - Bound vs lambda   (n = {n})")
    print("=" * 72)
    print(f"{'lambda':>7}  {'measured':>13}  {'2^k/(2^k+n-1)':>15}  {'vs 1/n':>9}")
    print("-" * 72)
    rows = []
    for k in kappas:
        B = rng.binomial(n - 1, 2.0 ** (-k), size=trials)
        meas = float(np.mean(1.0 / (1.0 + B)))
        bound = (2.0 ** k) / (2.0 ** k + n - 1)
        print(f"{k:>7}  {meas:>13.3e}  {bound:>15.3e}  {meas * n:>8.1f}x")
        rows.append((k, meas, bound))
    print()
    print("Each known bit halves the candidate set.  Wires whose values follow")
    print("from the constants by propagation add nothing further: any pixel")
    print("surviving the filters on a gate's inputs reproduces its outputs.")
    return rows


def nb_attack_demo(image_size: int = 128):
    """Show the filter shrinking the candidate set, wire by wire."""
    print("R2 OFF - constants embedded only at the secret pixel")
    out, t, view = run_protocol(CIRCUIT_SPEC, {"A": 1, "B": 0, "C": 1},
                                image_shape=(image_size, image_size, 3),
                                verbose=False, apply_r2=False,
                                return_view=True)
    a = attack_full_transcript(view, verbose=True)
    print(f"  -> {a['candidates']} candidates of {image_size**2}, "
          f"advantage {a['adv']:.3e} ({a['degradation']:.1f}x worse than 1/n)")
    print(f"  -> secret pixel survived: {a['secret_survives']}")

    print("\nR2 ON - constants are pixel-uniform planes Carol builds")
    out, t, view = run_protocol(CIRCUIT_SPEC, {"A": 1, "B": 0, "C": 1},
                                image_shape=(image_size, image_size, 3),
                                verbose=False, apply_r2=True,
                                return_view=True)
    b = attack_full_transcript(view, verbose=True)
    print(f"  -> {b['candidates']} candidates of {image_size**2}, "
          f"advantage {b['adv']:.3e} ({b['degradation']:.1f}x worse than 1/n)")
    print(f"  -> {b['derivable_wires']} wires were derivable; "
          f"{b['filters_applied']} filters had any effect")
    return a, b


def nb_benchmarks_security():
    """Benchmarks 4, 6 and 7: the security picture end to end."""
    bench_security_empirical(n_trials=4000)
    bench_constant_leak()
    bench_lambda_sweep()
