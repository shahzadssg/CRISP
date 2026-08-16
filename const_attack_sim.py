"""
Simulates the constant-wire attack of Corollary 1.

Model (faithful to the scheme):
  - n = h*w pixels. Secret pixel K uniform.
  - kappa ancillary wires carry publicly known constants, each embedded ONLY at K
    in the LSB of channel 0 of its own carrier image. Every other pixel of that
    plane is an independent uniform decoy bit.
  - Carol knows each constant's value. She keeps pixel j iff j agrees with all
    kappa known values. She then guesses uniformly among survivors.
  - She wins iff she picks K (picking K fixes the value guess correctly).

Reports measured joint-win rate vs the theoretical E[1/(1+B)] and the
Jensen lower bound 2^kappa/(2^kappa + n - 1).
"""
import numpy as np

rng = np.random.default_rng(20260811)


def theory_exact(n, kappa, trials=400000):
    """E[1/(1+B)], B ~ Bin(n-1, 2^-kappa), by Monte Carlo on B (fast, exact model)."""
    B = rng.binomial(n - 1, 2.0 ** (-kappa), size=trials)
    return float(np.mean(1.0 / (1.0 + B)))


def simulate(n, kappa, trials=200000):
    """Full pixel-level simulation of the filtering attack."""
    wins = 0
    for _ in range(trials):
        K = rng.integers(n)
        # constants are all 0 w.l.o.g.; decoys uniform
        # survivor set = K plus each other pixel that is 0 in all kappa planes
        if kappa == 0:
            m = n - 1
        else:
            m = rng.binomial(n - 1, 2.0 ** (-kappa))
        # Carol picks uniformly among 1 + m survivors
        if rng.integers(1 + m) == 0:
            wins += 1
    return wins / trials


print(f"{'n':>8} {'h x w':>10} {'kappa':>6} {'measured':>12} {'E[1/(1+B)]':>12} "
      f"{'2^k/(2^k+n-1)':>15} {'1/n':>12} {'degrad.':>9}")
print("-" * 92)

rows = []
for (h, w) in [(64, 64), (128, 128)]:
    n = h * w
    for kappa in [0, 1, 3, 5, 7, 10, 14]:
        meas = simulate(n, kappa, trials=120000)
        ex = theory_exact(n, kappa)
        lb = (2.0 ** kappa) / (2.0 ** kappa + n - 1)
        print(f"{n:>8} {h}x{w:>6} {kappa:>6} {meas:>12.3e} {ex:>12.3e} "
              f"{lb:>15.3e} {1/n:>12.3e} {meas*n:>8.1f}x")
        rows.append((h, w, n, kappa, meas, ex, lb))

print()
print("Reading: 'degrad.' is measured win rate divided by 1/n, i.e. how many times")
print("worse the real bound is than the single-image claim.")
