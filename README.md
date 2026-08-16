# CRISP: Channel-Randomised Single-Image Steganography with Permutations

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

CRISP is a homomorphic steganography scheme. Hidden bits live at a secret pixel of a cover image, and an honest-but-curious cloud server computes Boolean circuits over those bits without learning which pixel is the secret one. It extends ProSt with a single-image, dual-permutation construction that cuts image overhead from one cover per wire to two covers per gate.

The security analysis has two parts, and the second is a negative result. Positional hiding holds, with an exact bound over the full multi-image transcript. Circuit privacy does not hold at all, and we prove it.

Reference implementation for the paper appearing at **ICICS 2026** (Fukui, Japan, Springer LNCS). Extended version with full proofs: [ePrint 2026/424](https://eprint.iacr.org/2026/424).

---

## What this repository contains

```
.
├── code/
│   ├── CRISP.ipynb                  ← annotated walkthrough, executed, with figures
│   ├── CRISP_notebook_module.py     ← the verbose implementation the notebook uses
│   ├── CRISP_protocol_faithful.py   ← minimal core, no logging or plotting
│   └── const_attack_sim.py          ← standalone Monte Carlo over the lambda bound
```

---

## The idea, briefly

ProSt embeds one secret bit per cover image. A circuit with 47 wires needs 47 images.

CRISP puts the three logical inputs of a Fredkin gate (control, x, y) into the three RGB channels of one cover image at a secret pixel `(row, col)`, and the three outputs into a fresh cover at the same pixel. Two permutations `pi_in` and `pi_out`, sampled per gate from S₃, map channels to logical roles. Both are public, because the cloud needs `pi_in` to read the gate's inputs and `pi_out` to route its outputs to the next gate.

Image counts: 31 for the 14-gate benchmark against ProSt's 47, since repair R2 removes the two constant carriers from the transmission. Counting those carriers, as the paper's conservative figures do, gives 33 against 47 (1.42×) and 44 against 61 (1.39×) for a 17-gate variant.

One caveat. The 14-gate figure is the Fredkin synthesis published with ProSt, not a lower bound; a tighter synthesis narrows the gap, since ProSt's count scales with wires and CRISP's with gates. Read the result as a crossover, not a fixed speedup.

---

## Threat model (read this before opening issues)

A common misreading of homomorphic steganography is that it should hide the existence of stego content. That is incompatible with outsourced computation: if the cloud does not know it is running CRISP, it cannot run CRISP. The model is Kerckhoffs-style.

| Carol **sees** | Carol **does not see** |
|---|---|
| The full circuit specification | The secret pixel position `(row, col)` |
| Every per-gate `(pi_in, pi_out)` | The secret bit values, to the extent they are not implied by her prior |
| Every cover image and stego image | |
| The channel-to-role map at every gate | |
| The wiring graph and the gate count | |
| That LSB steganography is in use | |

The right-hand column has one entry that matters. Everything about the circuit is visible.

### What CRISP guarantees

**Positional hiding over the full transcript** (Theorem 3). Carol's optimal joint recovery of the secret pixel and the secret input satisfies

```
Adv_joint = E[1 / (1 + B)]  >=  2^lambda / (2^lambda + n - 1),
    B ~ Binomial(n - 1, 2^-lambda),  n = h * w,  lambda = m - H2(q)
```

where `m` is the number of secret input bits and `H2(q)` is the collision entropy of Carol's prior on them. `lambda` is the number of bits of the secret input Carol already knows, and each one halves the candidate pixel set. With `lambda = 0` this is exactly `1/(h*w)`.

**Correctness** (Theorem 1) and **distribution preservation** (Lemma 4) hold unconditionally.

### What CRISP does not guarantee

**Circuit privacy. Not in any sense, scoped or otherwise.**

Two proofs, in the paper as Propositions 2 and 3. Carol holds `pi_in` and `pi_out` because she needs them to evaluate and to route, so she inverts them in constant time and reads the channel-to-role map at every gate. Resampling `pi_out` per gate changes nothing; a permutation the adversary is handed is not a secret at any resampling rate. An adversary holding only the images and no specification recovers the same map plus the wiring graph, using nine channel-pair XOR comparisons per gate at `O(h*w)` bit operations each, with error at most `(3/4)^(h*w)`. The Fredkin control passes through unchanged, so the true channel pair sums to zero and every other pair disagrees with probability at least 1/4.

If you came here looking for function hiding, this scheme does not provide it. See the discussion of universal circuits and oblivious routing in the paper (Remark 1) for what it would cost.

---

## Two repairs the security analysis requires

The bound above does not hold for the scheme as originally described. Both repairs are cheap and neither affects correctness.

**R1: plane re-randomisation.** Before writing the secret bit at `(row, col)`, overwrite the entire LSB plane of channel 0 of each secret-carrying source image with fresh uniform bits. Without R1 the decoy bits are the cover's natural LSBs, which are spatially correlated and biased, while the embedded bit is uniform and independent. That mismatch is exactly the signal RS and sample-pair steganalysis exploit, and it means the i.i.d. assumption in the analysis is an assumption about photographs rather than a property of the scheme. R1 makes it true by construction. One pass over the image, no extra communication.

**R2: pixel-uniform constant planes.** A wire whose value is fixed and public carries that value in the LSB of its channel at *every* pixel, not only at the secret one. Carol materialises the plane locally from the circuit specification, so no carrier is transmitted for it.

R2 is not optional. Without it, a constant embedded only at the secret pixel lets Carol discard every pixel whose LSB disagrees with the value she already knows, and each constant carrier adds one bit to `lambda`.

Measured end to end against this implementation, 128×128 covers, 60 executions per row. `attack_full_transcript()` plays Carol: it constant-propagates the public specification, then filters candidate pixels on every wire whose value it derives.

| circuit | R2 | κ | images sent | mean candidates | joint advantage | vs 1/n |
|---|---|---|---|---|---|---|
| benchmark, 14 gates | off | 2 | 33 | 4096.4 | 2.44e-4 | 4.0× |
| benchmark, 14 gates | **on** | 0 | **31** | 16384.0 | 6.10e-5 | **1.0×** |

The benchmark carries two constants, so κ = 2 and the candidate set shrinks fourfold. R2 restores exactly `1/n` and drops two images from the transmission while it is at it. A circuit padded with more constant wires pays proportionally: the sweep below covers κ up to 14.

| κ | measured | 2^κ/(2^κ+n−1) | vs 1/n |
|---|---|---|---|
| 0 | 6.10e-5 | 6.10e-5 | 1.0× |
| 3 | 4.88e-4 | 4.88e-4 | 8.0× |
| 7 | 7.81e-3 | 7.75e-3 | **128×** |
| 14 | 6.32e-1 | 5.00e-1 | **10356×** |

At κ = 14 Carol wins outright most of the time.

One subtlety: wires whose values follow from the constants by propagation add nothing further. Any pixel that already survives the filters on a gate's constant inputs reproduces that gate's outputs automatically, so those constraints are implied rather than independent. Only the count of **distinct constant carriers** matters, which `known_constant_wires()` and `attack_full_transcript()` demonstrate together.

Reproduce with `demo_repairs()` or the standalone `const_attack_sim.py`.

The cost of R2 is that a constant plane has a degenerate LSB channel, so it is not undetectable. That is acceptable because the plane carries no secret and its contents are already public, but it does mean the transmitted bundle is easy to spot.

Both repairs are on by default. `run_crisp(..., apply_r1=False, apply_r2=False)` disables them so the leak can be measured; that is the only reason the flags exist. All 288 unit tests and the 8-row truth table pass with the repairs enabled.

`CRISP.ipynb` is generated from the same code and runs end to end with every cell executed.

---

## Quickstart

### Install

```bash
git clone https://github.com/<your-username>/CRISP.git
cd CRISP/code
pip install numpy Pillow matplotlib
```

Python 3.9+. No GPU, no compilation.

### Run the walkthrough

```bash
jupyter notebook CRISP.ipynb
```

Ships already executed, so the outputs and figures render without running anything. It walks through the circuit model and primitives, prints a full verbose run showing every gate with its permutations, wire carriers and per-step timing, displays the cover and stego images and their LSB planes, renders each gate's input and output inline, traces the constant-wire attack filter by filter with R2 off and on, and finishes with seven benchmark suites: per-gate timing, image count, chi-square steganalysis, positional security at λ = 0, the constant-wire leak, the λ sweep, and end-to-end timing.

Set `SHOW_FIGURES = False` in the first cell for headless runs, which saves figures to `/tmp` instead.

### Run the minimal core

```bash
python CRISP_protocol_faithful.py
```

### Use your own cover images

Set `COVER_DIR` in the notebook (section 5b) or pass `cover_dir=` to `run_protocol`:

```python
out, timing = run_protocol(CIRCUIT_SPEC, {"A": 1, "B": 0, "C": 1},
                           image_shape=(128, 128, 3),
                           cover_dir="C:/cover")
```

Files load in sorted filename order. `load_cover_pool(..., resize_to=(h, w))` centre-crops and resizes a mixed-size album to a common shape, which every cover must share since Carol works pixel by pixel across the whole set.

**How many images.** Two covers per gate, one for the assembled input and one for the output, plus one per secret-valued source wire. R2 means constant wires cost nothing. For the 14-gate benchmark that is 14 + 14 + 3 = **31**. Call `cover_requirement(spec)` for any circuit.

**Never save results as JPEG.** JPEG quantisation destroys exactly the bits the scheme stores. `save_image()` writes PNG.

### Are your covers suitable?

`cover_report(covers)` measures whether the LSB planes are uniform and independent, which is what the analysis assumes of the decoy bits. Whether an album satisfies that on its own depends entirely on the pictures:

| cover type | set bits | neighbour agreement | after R1 |
|---|---|---|---|
| flat / screenshot-like | 0.5000 | 0.9921 | 0.5050 |
| noiseless gradient | **0.0000** | **1.0000** | 0.5031 |
| heavily compressed JPEG | 0.5000 | 0.9921 | 0.4995 |
| noisy camera photograph | 0.5029 | 0.5007 | 0.4976 |

Neighbour agreement is how often horizontally adjacent LSBs match, and sits at 0.5 for an independent plane. Noisy photographs are already fine, because sensor noise randomises the low bit almost everywhere. Graphics, screenshots, renders and flat regions are not: the noiseless gradient has an LSB plane that is **entirely zero**, so a single embedded bit would be one lit pixel on a black field.

R1 removes the question for one pass over the image, which is why it is on by default rather than conditional on the album.

### Notes for a real album

- **Already 128x128?** Then `resize_to` is skipped. `load_cover_pool` short-circuits when an image already matches the target shape, so no resampling touches your pixels.
- **Windows paths work.** Temp files use `tempfile.gettempdir()`, so `SHOW_FIGURES = False` writes to `%TEMP%` rather than a POSIX path.
- **You need 31 images minimum** for the 14-gate benchmark. More in the folder is fine; the first 31 in sorted filename order are used.
- **Covers must not be reused across executions.** Two runs on one cover differ only where the LSB planes disagree, which hands Carol the secret position.

---

## Reference implementation

| File | Lines | Purpose |
|---|---|---|
| `CRISP_protocol_faithful.py` | ~500 | Minimal core. Use this for code review and citation. |
| `CRISP.ipynb` | 54 cells | Annotated walkthrough with per-gate `[STAGE]` logging, inline figures and seven benchmark suites. Ships executed with outputs. |
| `CRISP_notebook_module.py` | ~1400 | The verbose implementation: same primitives plus logging, image I/O and display hooks. |
| `const_attack_sim.py` | ~60 | Standalone Monte Carlo over the `lambda` bound. |

Attack and demo entry points in the core module: `known_constant_wires()`, `attack_constant_wires()`, `attack_full_transcript()`, `demo_repairs()`.

### The primitives

| Function | Who calls it | What it does |
|---|---|---|
| `locgen(image_shape)` | Alice | Sample `(row, col)` uniformly. |
| `perm_gen()` | Alice | Sample `(pi_in, pi_out)` for one gate. |
| `emb_bit(image, bit, row, col, ch)` | Alice | Write one LSB at the secret pixel. |
| `assemble_stego(input_cover, wire_sources, pi_in)` | Carol | Per-pixel channel routing. |
| `comp(stego, output_cover, pi_in, pi_out)` | Carol | Pixel-wise Fredkin pass. |
| `ext_bit(image, row, col, ch)` | Bob | Read one LSB at the secret pixel. |

### Verifying that Carol is position-blind

The signatures are the proof:

```python
import inspect
from CRISP_protocol_faithful import Carol_evaluate, Bob_extract

print(inspect.signature(Carol_evaluate))
# (circuit, gate_order, wire_sources, input_covers, output_covers,
#  permutations, processed_dir=None, inline_display_enabled=False, verbose=True)

print(inspect.signature(Bob_extract))
# (circuit, wire_to_image, row, col, verbose=True)
```

`Carol_evaluate` takes neither `row` nor `col` and never calls `ext_bit`. The only path from the secret position to any output runs through `Bob_extract`, called once per primary output wire.

---

## Reproducing the paper's results

| Paper claim | Where to verify it |
|---|---|
| 288 unit tests pass | `nb_unit_tests()` |
| Truth table verifies on the protocol-faithful pipeline | `nb_truth_table()` |
| Image counts: 33 vs 47, and 44 vs 61 | `bench_image_count()` |
| Joint recovery matches `1/(h*w)` when `lambda = 0` | `bench_security_empirical()` |
| Constant wires degrade the bound to `2^kappa/(2^kappa + n - 1)` | `const_attack_sim.py` |
| Chi-square does not separate CRISP outputs from synthetic covers | `bench_chi_square_steganalysis()` |
| Per-gate `comp` dominates cost, `assemble` is roughly half | `bench_per_gate_timing()` |
| End-to-end timing on the 14-gate benchmark | `bench_end_to_end()` |

Note that the ProSt comparison column in the paper's timing table is a cost model charging `3 × comp` per gate, not a measurement. A like-for-like timing against the released ProSt implementation has not been run.

---

## Limitations and known caveats

Flagged because they are real and the paper says so.

1. **No circuit privacy.** Proven, not conjectured. See above and Propositions 2 and 3.

2. **The advantage is polynomially small, not negligible.** It decays as `1/n` in the pixel count, so the bound is statistically small but not cryptographically negligible in the standard sense. Scaling `n` costs bandwidth linearly.

3. **Communication is the real cost.** The 17-gate benchmark at 128×128 moves about 2.06 MB against ProSt's 2.86 MB. That is 2 MB to evaluate a three-input Boolean function. TFHE moves tens of kilobytes and a TEE moves about one. The image-count reduction improves a constant on a cost that is already three orders of magnitude above the alternatives. If you are choosing on bandwidth, choose FHE.

4. **Synthetic covers, not natural images.** The chi-square benchmark uses random covers, which is the easy case for any LSB test. Passing it is necessary, not sufficient. On natural photographs the question is whether R1's re-randomised plane looks like a natural LSB plane, and it does not. Sample-pair analysis and learned detectors (SRM, Yedroudj-Net) on BOSSbase or ALASKA are untested, and we expect them to separate CRISP images from natural covers. This does not affect Theorem 3, which holds under known presence, but CRISP should not be called undetectable.

6. **Honest-but-curious only.** A malicious Carol can corrupt outputs. The paper sketches a Hamming-weight MAC integrity layer but does not prove it. Open.

7. **Fresh keys and fresh covers per execution.** Reusing `(row, col)` enables a differential attack across runs, and reusing a cover lets Carol difference two runs at every pixel. Reuse should be rejected at the API boundary.

8. **Linear-chain implementation.** `Carol_evaluate` handles linear chains. General DAGs need a wire-tracking layer that adds bookkeeping without changing the security analysis. Fan-out itself is free, since a wire is an (image, channel) pair that `assemble_stego` reads without consuming.

9. **CPU-only.** `comp` is embarrassingly parallel and a GPU port should give a large speedup. Not implemented.

---

## Building on prior work

CRISP extends **ProSt**: S. Ahmad and S. Rass, "ProSt: Provably Secure Homomorphic Steganography," *IEEE Access*, vol. 14, pp. 14393–14409, 2026. [doi:10.1109/ACCESS.2026.3656995](https://doi.org/10.1109/ACCESS.2026.3656995)

CRISP inherits ProSt's circuit-specification format and randomised topological sort, and replaces the per-wire image model with a per-gate one. The Fredkin construction, the secret-position threat model, and the Bayesian style of the security proofs are adaptations of the ProSt framework.

ProSt also ships a pass that randomises wire names and pads the circuit with extra ancillary wires and gates. CRISP does not use it and claims nothing for it. Proposition 3 recovers the wiring graph from the images alone without reading a single wire name, so relabelling changes nothing an adversary can act on, and the extra constant wires it introduces cost positional security rather than buying anything.

Repair R2 applies to ProSt as well, since it embeds ancillary constants at the secret pixel too. Corollary 1 is not a CRISP-specific patch.

---

## Citation

```bibtex
@inproceedings{crisp2026,
  author    = {Ahmad, Shahzad and Rass, Stefan},
  title     = {{CRISP}: Channel-Randomised Single-Image Steganography with Permutations},
  booktitle = {Information and Communications Security (ICICS 2026)},
  series    = {Lecture Notes in Computer Science},
  publisher = {Springer},
  year      = {2026}
}
```

## License

LIT Secure and Correct Systems Lab, Johannes Kepler University Linz, Austria
