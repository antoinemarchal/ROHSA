"""Compare fits with the true noiseless cube. Usage: python compare.py <data dir> <name>...
       python compare.py <data dir> --diff <name A> <name B>   (rms and max abs difference of the two fitted models)

Prints, per run, the rms difference between the fitted model and the truth over
all channels, for the blank pixels and for the other pixels.
"""
import sys
import numpy as np

d = sys.argv[1]
truth = np.load(f"{d}/truth.npy")
blank = np.load(f"{d}/blank.npy")
nv, ny, nx = truth.shape
v = np.arange(nv)[:, None] + 1.0


def model(name):
    p = np.loadtxt(f"{d}/runs/{name}.out", comments="#").reshape(ny, nx, -1, 5)
    m = np.zeros_like(truth)
    for g in range(p.shape[2]):
        a, mu, s = (p[:, :, g, k + 2][None] for k in range(3))
        m += a * np.exp(-(v.reshape(nv, 1, 1) - mu) ** 2 / (2 * s ** 2))
    return m


if sys.argv[2] == "--diff":
    r = model(sys.argv[3]) - model(sys.argv[4])
    print(f"{sys.argv[3]} - {sys.argv[4]}: rms {np.sqrt((r ** 2).mean()):.2e}  max abs {np.abs(r).max():.2e}")
    sys.exit()

for name in sys.argv[2:]:
    r = model(name) - truth
    f = lambda mask: np.sqrt((r[:, mask] ** 2).mean())
    print(f"{name:24s} blank {f(blank):.3f}  other {f(~blank):.3f}")
