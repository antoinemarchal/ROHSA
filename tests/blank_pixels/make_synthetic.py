"""Synthetic cube for the blank pixels test (numpy only).

A 32 x 32 x 64 cube with two smooth Gaussian components plus noise, written in
ROHSA's text format together with noise maps:

  cube.dat        intact cube
  cube_blank.dat  the same cube with a 6 x 6 patch and a few single pixels set to 0 in every channel
  rms.dat         noise map, 0.05 everywhere
  rms_blank.dat   noise map with 0 at the blank pixels
  truth.npy, blank.npy  the noiseless cube and the mask of blank pixels

Usage: python make_synthetic.py <out dir> [NYxNX]

With NYxNX (e.g. 30x28) every product is cropped to the first NY x NX pixels of the 32 x 32 map, to get a
map size that is not a power of two (padded grid in ROHSA). The truth is the crop of the true cube.
"""
import os, sys
import numpy as np

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(1)
nv, ny, nx = 64, 32, 32
y, x = np.mgrid[0:ny, 0:nx] / (ny - 1.0)
v = np.arange(nv, dtype=float)[:, None, None] + 1.0   # ROHSA channels are 1-based

truth = np.zeros((nv, ny, nx))
for amp, mu, sig in [(2.0 + 0.5 * x, 18.0 + 4.0 * y, 3.0 + 1.0 * x),
                     (1.5 + 0.4 * y, 42.0 - 3.0 * x, 5.0 + 1.0 * y)]:
    truth += amp * np.exp(-(v - mu) ** 2 / (2.0 * sig ** 2))
sigma = 0.05
cube = truth + sigma * rng.standard_normal(truth.shape)

blank = np.zeros((ny, nx), dtype=bool)
blank[10:16, 12:18] = True                      # 6 x 6 patch
for j, i in [(3, 3), (20, 25), (27, 8), (5, 28), (0, 31)]:
    blank[j, i] = True                          # single pixels, one at the map edge
cube_blank = cube.copy()
cube_blank[:, blank] = 0.0

if len(sys.argv) > 2:
    cy, cx = (int(t) for t in sys.argv[2].lower().split("x"))
    truth, cube, cube_blank, blank = (a[..., :cy, :cx] for a in (truth, cube, cube_blank, blank))
    ny, nx = cy, cx


def write_cube(name, c):
    with open(os.path.join(out, name), "w") as f:
        f.write(f"{nv} {ny} {nx}\n")
        for k in range(nv):
            for j in range(ny):
                for i in range(nx):
                    f.write(f"{k} {j} {i} {c[k, j, i]:.16e}\n")


def write_map(name, m):
    with open(os.path.join(out, name), "w") as f:
        f.write(f"{ny} {nx}\n")
        for j in range(ny):
            for i in range(nx):
                f.write(f"{j} {i} {m[j, i]:.16e}\n")


write_cube("cube.dat", cube)
write_cube("cube_blank.dat", cube_blank)
noise = np.full((ny, nx), sigma)
write_map("rms.dat", noise)
nb = noise.copy()
nb[blank] = 0.0
write_map("rms_blank.dat", nb)
np.save(os.path.join(out, "truth.npy"), truth)
np.save(os.path.join(out, "blank.npy"), blank)
