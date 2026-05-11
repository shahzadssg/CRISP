# CRISP: Circuit-pRivate Single-Image Steganography with Permutations

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)


CRISP is a provably secure homomorphic steganography scheme. Hidden bits live at a secret pixel of a cover image; an honest-but-curious cloud server can compute Boolean circuits over those bits without learning which pixel is the secret one. CRISP advances ProSt on three fronts: a precisely scoped threat model, a single-image dual-permutation construction that drops image overhead, and a self-contained Bayesian security analysis with empirical confirmation.

---

## What this repository contains

```
.
├── code/
    ├── CRISP.ipynb                  ← Jupyter notebook walkthrough
    └── CRISP_protocol_faithful.py   ← minimal core (≈500 lines)
    
```


---

## The idea, briefly

ProSt embeds one secret bit per cover image. A circuit with 47 wires therefore needs 47 images, and the cloud trivially identifies the control wire of every Fredkin gate (since each wire's role is fixed at a specific image).

CRISP changes both:

1. **One image per gate, not per wire.** Three logical inputs (control, x, y) live in the three RGB channels of a single cover image at a secret pixel `(row, col)`. Three logical outputs live in a fresh output cover at the same pixel.

2. **Two independent permutations per gate.** `pi_in` and `pi_out` are sampled freshly per gate from S₃. Both are public (the cloud needs them to compute), but their independence per gate is what defeats the per-channel correlation attack the Fredkin gate's control-bit pass-through would otherwise enable.

The result: 33 images for the 14-gate benchmark vs ProSt's 47 (1.42× reduction); 44 images for the 17-gate obfuscated form vs ProSt's 61 (1.39×); larger circuits amortise to ~1.49×. Joint position+value recovery probability is exactly **1/(h·w)**, independent of the secret-bit prior.

---

## Threat model (read this before opening issues)

A common misreading of homomorphic steganography is that it should hide the existence of stego content. That is incompatible with outsourced computation: if the cloud doesn't know it's running CRISP, it can't run CRISP. The threat model is therefore Kerckhoffs-style:

| Adversary (cloud Carol) **sees** | Adversary **does not see** |
|---|---|
| The full circuit specification | The secret pixel position `(row, col)` |
| Every per-gate `(pi_in, pi_out)` | The secret bit values |
| Every cover image and stego image | The wire-name to logical-purpose mapping (obfuscated) |
| The fact that LSB steganography is in use | |

CRISP guarantees three things within this model:

1. **Positional hiding.** Carol's posterior over the secret pixel stays uniform (Theorem 3).
2. **Value hiding through positional uncertainty.** Joint position+value recovery rate is exactly 1/(h·w) (Theorem 4).
3. **Per-channel circuit privacy.** Independent per-gate `pi_out` defeats the LSB-correlation attack on the control wire (Theorem 2 / Necessity).

CRISP does **not** hide the existence of stego content, the gate count, or the circuit topology beyond wire-role assignments. The last is handled by the obfuscation pass inherited from ProSt.

---

## Quickstart

### Install

```bash
git clone https://github.com/<your-username>/CRISP.git
cd CRISP
pip install numpy Pillow matplotlib
```

Python 3.9+ is required. No GPU, no compilation.

### Run the verbose walkthrough (script mode)

```bash
python CRISP.ipynb
```

This runs the 288 unit tests, obfuscates the benchmark circuit, prompts you for the secret position and the input bits, then walks through Alice's setup, Carol's gate-by-gate evaluation, and Bob's extraction with detailed `[STAGE]` log lines and per-gate matplotlib displays.

### Run from a Jupyter notebook

The script's `argparse` would otherwise choke on Jupyter's kernel arguments, so the file ships with notebook-friendly helpers:

```python
from CRISP_illustrative import (
    nb_unit_tests, nb_demo, nb_truth_table,
    nb_obfuscated_demo, nb_benchmarks_quick, nb_benchmarks_full,
)

nb_unit_tests()              # 288 unit tests
nb_demo(image_size=128)      # single end-to-end run, verbose
nb_truth_table()             # all 8 truth-table rows
nb_obfuscated_demo()         # obfuscation pipeline + end-to-end
nb_benchmarks_quick()        # ~30 seconds
nb_benchmarks_full()         # ~3 minutes, paper-quality numbers
```

Set `SHOULD_DISPLAY_INLINE = True` at the top of the file to render gate-by-gate matplotlib figures inline; set it `False` for headless runs (figures get saved to `/tmp/`).

### Use your own cover images

If `C:/cover/` exists and contains at least `2 * n_gates + n_sources` PNG/JPG files, the script loads them in alphabetical order. Otherwise it generates synthetic 128×128 covers. To use a different directory, edit the `cover_dir` line in `main()`.

---

## Reference implementation

Two files, same scheme, different verbosity:

| File | Lines | Purpose |
|---|---|---|
| `CRISP_protocol_faithful.py` | ~500 | Minimal core. Use this for code review and citation. |
| `CRISP.ipynb` | ~1500 | Verbose. Same primitives, extensive `[STAGE]` logging, per-gate display, five benchmark suites. Use this for understanding the protocol. |

Both implement the same algorithms and produce identical outputs.

### The five primitives

| Function | Who calls it | What it does |
|---|---|---|
| `locgen(image_shape)` | Alice | Sample `(row, col)` uniformly. |
| `perm_gen()` | Alice | Sample `(pi_in, pi_out)` for one gate. |
| `emb_bit(image, bit, row, col, ch)` | Alice | Write one LSB at the secret pixel. |
| `assemble_stego(input_cover, wire_sources, pi_in)` | Carol | Per-pixel channel routing (Algorithm 5). |
| `comp(stego, output_cover, pi_in, pi_out)` | Carol | Pixel-wise Fredkin pass. |
| `ext_bit(image, row, col, ch)` | Bob | Read one LSB at the secret pixel. |

### Verifying that Carol is position-blind

The function signature itself is the proof:

```python
import inspect
from CRISP_illustrative import Carol_evaluate, Bob_extract

print(inspect.signature(Carol_evaluate))
# (circuit, gate_order, wire_sources, input_covers, output_covers,
#  permutations, processed_dir=None, inline_display_enabled=False, verbose=True)

print(inspect.signature(Bob_extract))
# (circuit, wire_to_image, row, col, verbose=True)
```

`Carol_evaluate` does not take `row` or `col` as arguments and does not call `ext_bit`. The only path from the secret position to any output flows through `Bob_extract`, which is called exactly once per primary output wire.

---

## Reproducing the paper's results

| Paper claim | Where to verify it |
|---|---|
| 288 unit tests pass | `nb_unit_tests()` or `python CRISP_illustrative.py` |
| Truth table verifies on the protocol-faithful pipeline | `nb_truth_table()` |
| Image-count reduction: 1.39× obfuscated, 1.42× un-obfuscated, ~1.49× at 100 gates | `bench_image_count()` (Benchmark 2) |
| Joint position+value rate matches 1/(h·w) within sampling error | `bench_security_empirical()` (Benchmark 4) |
| CRISP outputs are statistically indistinguishable from natural covers under chi-square detection | `bench_chi_square_steganalysis()` (Benchmark 3) |
| Per-gate Comp dominates total cost; Assemble is roughly half | `bench_per_gate_timing()` (Benchmark 1) |
| End-to-end timing on 14-gate benchmark | `bench_end_to_end()` (Benchmark 5) |

Quick benchmarks: `nb_benchmarks_quick()` (~30 seconds). Paper-quality: `nb_benchmarks_full()` (~3 minutes on a laptop).

---

## Limitations and known caveats

These are flagged because they are real and the paper says so:

1. **Synthetic covers, not natural images.** The chi-square steganalysis benchmark uses random covers. Passing it is a *necessary* condition for steganalysis resistance, not a sufficient one. Detection by trained CNN steganalysers (SRM, Yedroudj-Net) on natural-image corpora (BOSSbase, ALASKA) has not been evaluated.

2. **Honest-but-curious adversary only.** A malicious Carol can corrupt outputs. The paper sketches a Hamming-weight MAC integrity layer (Lemma 5 preserves the per-pixel weight multiset) but does not prove it; this is open work.

3. **Per-channel circuit privacy, not simulation-based.** Theorem 2 rules out the LSB-correlation attack that fixed-`pi_out` would enable; it does not give garbled-circuits-style simulation security. Lifting to UC-style privacy is open.

4. **Linear-chain implementation.** The current `Carol_evaluate` handles linear chains; general DAGs require a wire-tracking layer that does not change the security analysis but does add bookkeeping.

5. **CPU-only.** `comp` is embarrassingly parallel; a GPU port is straightforward and should give 10-100× speedup. Not yet implemented.

---


## Building on prior work

CRISP is a direct extension of:

- **ProSt** (the predecessor): Provably Secure Homomorphic Steganography. CRISP inherits ProSt's circuit-specification format, obfuscation pipeline, and randomised topological sort verbatim, and replaces ProSt's per-wire image model with a per-gate one.

The Fredkin gate construction, secret-position threat model, and Bayesian-style security proofs are direct adaptations from the ProSt framework.

