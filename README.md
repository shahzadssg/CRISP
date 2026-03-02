
# CRISP: Circuit-pRivate Single-Image Steganography with Permutations

## Overview

CRISP (**C**ircuit-p**R**ivate Single-**I**mage **S**teganography with **P**ermutations) is a Python implementation of a homomorphic steganography scheme designed for secure circuit evaluation in cloud environments. It builds upon the concepts introduced in ProSt, optimizing image usage and enhancing circuit privacy through the use of fully symmetric input and output permutations per gate.

The core idea is to embed binary circuit inputs into a single pixel of an image (the *cover image*) by modifying the least significant bits (LSBs) of its RGB channels. A cloud server Carol, without knowing the secret pixel location, can then perform homomorphic operations across the entire image. The result of the computation is encoded back into the LSBs of a *new, independent output cover image*. The authorized recipient (Bob) can then extract the circuit output from the secret pixel of this processed image.

## Key Features & Differences from ProSt

CRISP shares several core mechanisms with its predecessor, ProSt:

*   **Circuit Specification:** Uses the same 17-gate Fredkin circuit structure.
*   **Obfuscation Pipeline:** Incorporates wire renaming and the addition of dummy ancillaries and gates to obscure the circuit structure.
*   **Randomized Topological Sort:** Gates are executed in a topologically valid but randomized order to further hide computation flow.
*   **Interactive Input:** Prompts the user for circuit input bits and the secret pixel position.
*   **Output Messages & Profiling:** Provides detailed console output and performance metrics.

The **key innovation** in CRISP lies in its image utilization and permutation strategy:

*   **Optimized Image Usage:** ProSt requires one image per wire, leading to 61 images for a 17-gate obfuscated circuit. **CRISP requires only one image per gate**, reducing the image count to 34 for the same circuit (a 1.8x reduction). This significantly lowers storage and transmission overhead.
*   **Enhanced Circuit Privacy with Symmetric Permutations:**
    *   Each Fredkin gate embeds its three logical inputs (control `c`, input `x`, input `y`) into the RGB channels of a *single image* at a secret pixel.
    *   This embedding uses a secret input permutation $\pi_{in} \in S_3$.
    *   The cloud (Carol) applies the Fredkin gate pixel-wise across the entire image using the public $(\pi_{in}, \pi_{out})$ for that specific gate.
    *   Bob reads the three output bits from the same secret pixel using the corresponding $\pi_{out}$.
    *   Both $\pi_{in}$ and $\pi_{out}$ are sampled **independently and uniformly at random from $S_3$** for *every single gate*.
    *   Neither permutation is anchored to any user-chosen channel. This "fully symmetric design" ensures that the channel-to-role assignment is entirely fresh at both the input and output sides of every gate, maximizing circuit privacy.

## Security Considerations (Adapted from ProSt Theorem 1)

The adversary's advantage against CRISP is $\Theta(1/n)$, where $n = h \times w$ represents the total number of possible pixel positions (height × width). In contrast, ProSt used $n = h \times w \times 3$ for pixel-channel positions.

While CRISP's position space is 3 times smaller, it compensates by:
1.  Using 1.8 times fewer images overall.
2.  Adding significant circuit privacy through the per-gate, randomly chosen, fully symmetric input and output permutations.

Both schemes offer cryptographically negligible adversary advantages.

## Technical Details

### Fredkin Gate
The core logic gate used is the Fredkin gate: $F(c, x, y) \rightarrow (c, y \text{ if } c=1 \text{ else } x, x \text{ if } c=1 \text{ else } y)$. It passes the control bit $c$ through unchanged and swaps $x$ and $y$ if $c$ is 1.

### Steganography Primitives
*   `emb(image, c_bit, x_bit, y_bit, row, col, pi_in)`: Alice's function to embed three logical bits into the LSBs of a single pixel's RGB channels according to $\pi_{in}$. Returns a new stego image.
*   `comp(stego, output_cover, pi_in, pi_out)`: Carol's (Carol's) function to apply the Fredkin operation pixel-wise across the entire `stego` image. Inputs are read using $\pi_{in}$ and outputs are written into an independent `output_cover` image using $\pi_{out}$.
*   `ext(image, row, col, pi_out)`: Bob's function to extract the three logical output bits from the secret pixel of the processed image according to $\pi_{out}$.

### Circuit Evaluation Flow
1.  **Alice's Setup:**
    *   Defines the circuit (`CIRCUIT_SPEC`).
    *   Obfuscates the circuit (wire renaming, dummy components).
    *   Loads/generates `n_gates` input cover images and `n_gates` independent output cover images.
    *   Prompts for a secret pixel position `(row, col)`.
    *   Prompts for initial input bits (A, B, C) and embeds them conceptually (log prints mirror ProSt).
2.  **Carol's Computation (`Carol_compute_circuit`):**
    *   Parses the obfuscated circuit.
    *   Performs a randomized topological sort of the gates.
    *   Iterates through each gate in the randomized order:
        *   Samples new, independent $\pi_{in}$ and $\pi_{out}$ for the current gate.
        *   Uses `emb` to "virtually" embed input bits for the current gate into its assigned `input_cover` image (the actual embedding is part of the homomorphic operation for Carol).
        *   Uses `comp` to apply the Fredkin operation across the entire `input_cover` image, producing an `output_image` based on a fresh `output_cover`.
        *   Uses `ext` to extract the *logical* output bits from the `output_image` at the secret pixel.
        *   Stores these logical output bits for subsequent gates.
        *   Saves intermediate stego and output images.
        *   Prints a Base64 representation of the processed images for console output (or saves them to disk if Base64 is not desired for every step).
    *   Returns the final computed wire values.
3.  **Bob's Extraction:**
    *   Receives the final output value from Carol (conceptually extracted from the final image).

## Setup & Usage

### Prerequisites
*   Python 3.x
*   `numpy`
*   `Pillow` (PIL fork)
*   `matplotlib`

You can install these using pip:
```bash
pip install numpy Pillow matplotlib
```

### Directory Structure
The script expects a `cover` directory at the same level as the script for image inputs, or it will generate synthetic images.

.
├── crisp.py
└── C:/cover/  # (Or any other path, default is C:/cover)
    ├── image01.png
    ├── image02.png
    └── ...

If `C:/cover` does not exist or does not contain enough images, the script will generate synthetic 128x128 RGB images. You need at least `2 * num_gates` images if you provide your own (one set for input covers, one for output covers).

### Running the Script

1.  **Save the code:** Save the provided Python code as `crisp.py`.
2.  **Execute from terminal:**
    ```bash
    python crisp.py
    ```
3.  **Interactive Prompts:**
    *   The script will first run unit tests.
    *   It will then prompt you to enter a secret pixel position as `'row,column'` (e.g., `64,64` for a 128x128 image).
    *   Finally, it will ask you to enter binary inputs (0 or 1) for `A`, `B`, and `C` for the example circuit.
4.  **Output:**
    *   The console will display various logging messages, profiling information, and the final circuit output.
    *   During Carol's computation, intermediate images (stego and output for each gate) will be saved to `C:/cover/processed/` (or `/tmp/crisp_display_*.png` if `C:/cover` doesn't exist).
    *   Crucially, the console will also print **Base64 encoded strings** of the displayed images, making the image content part of the textual output stream. You can decode these Base64 strings externally to view the images.

### Circuit Specification

The example circuit `(A AND C) OR (NOT A AND B) OR (NOT B AND NOT C)` is hardcoded in `CIRCUIT_SPEC`. You can modify this string to experiment with different Fredkin-gate based circuits.

## Contributions & Further Research

This implementation is a research prototype. Feel free to explore:
*   **Performance Optimization:** Investigate more efficient NumPy operations or alternative image processing libraries.
*   **Advanced Obfuscation:** Explore different obfuscation techniques for the circuit structure.
*   **Alternative Gates:** Adapt the steganography primitives for other universal gate sets.
*   **Security Analysis:** Conduct more rigorous analysis of the scheme's security bounds and attack vectors, especially regarding the new permutation strategy.
*   **Practical Deployment:** Consider how this could integrate into a real-world cloud service, addressing issues like image size, bandwidth, and real-time performance.

---
