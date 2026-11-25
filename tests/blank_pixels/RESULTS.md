# Blank pixels test results

Synthetic cube 32 x 32 x 64, two smooth Gaussian components, noise 0.05
(`make_synthetic.py`). Blank pixels: a 6 x 6 patch plus 5 single pixels (one at
the map edge), 41 in total, set to 0 in every channel with noise 0 in the noise
map. rms of (fitted model - true noiseless cube) over all channels, gfortran
16.2, `make` defaults (-O2), master e36bc49.

| run | blank pixels | other pixels |
|---|---|---|
| (a) master, intact cube | 0.015 | 0.015 |
| (b) master, blanked cube, noise 0 at blank pixels | 0.951 | 0.015 |
| (c) this branch, blanked cube | 0.011 | 0.015 |

Run: `python3 make_synthetic.py D; ./run.sh BIN D NAME cube.dat rms.dat; python3 compare.py D NAME`
(`./run.sh BIN D NAME cube_blank.dat rms_blank.dat` for the blanked cube; use `none`
instead of the noise file to run with `noise = .false.`).

Other checks:

- (d) `noise = .false.`, intact and blanked cube: master and branch outputs are byte-identical (`cmp`).
- (e) noise map without zeros (`rms.dat`), intact and blanked cube: byte-identical.
- (f) `-fcheck=all -fbacktrace` builds of master and branch run (b)/(c) without bounds errors
  (rms 0.951 / 0.010 at the blank pixels, 0.015 elsewhere). The branch also runs the
  intact cube, and a 30 x 28 crop of the blanked cube (padded grid), without errors.

## Further checks

`./extra_checks.sh MASTER_BIN BRANCH_BIN WORKDIR` reproduces sections 1 to 3 (optional master variants, see its header).
Same metric as above: rms of (fitted model - true cube), blank pixels / other pixels. Binaries built with `make` defaults.

### 1. `lym = .true.` (n_gauss = 2, lambda_lym_sig = 10, 32 x 32)

| run | blank pixels | other pixels |
|---|---|---|
| master, intact cube, `rms.dat` | 0.015 | 0.016 |
| branch, intact cube, `rms.dat` | 0.015 | 0.016 |
| master, blanked cube | 0.951 | 0.016 |
| branch, blanked cube | 0.142 | 0.016 |

- Intact cube with a noise map without zeros: master and branch outputs are byte-identical.
- The branch is better than master on the blanked cube but worse than without lym (0.011), because the lym term pulls
  sigma2 / sigma1 of a blank pixel toward c_lym = 1. With lambda_lym_sig = 0 the numbers equal the non-lym case
  (blank 0.951 master, 0.011 branch; intact cube identical).
- `-fcheck=all -fbacktrace` builds of master and branch run the blanked lym case without errors.

### 2. Padded grids (`make_synthetic.py D 30x28`, `... 20x32`; truth = crop of the true cube)

(a) Intact cube, `rms.dat` (no zeros). Padded pixels have noise 0 on the padded grid.

| 30 x 28 | blank/all pixels |
|---|---|
| master, 3 runs, default heap | byte-identical to each other, 0.015 |
| master, `-finit-real=snan` and `=zero` builds | byte-identical to master |
| master, `MALLOC_PERTURB_=85` | 0.015 (differs from master by rms 1.9e-06) |
| master, `MALLOC_PERTURB_=170` | 0.573 (differs from master by rms 0.58, max abs 1.7) |
| master + `std_cube = 0`, `MALLOC_PERTURB_=170` | still differs (rms 0.58) |
| master + `std_cube = 0` and `cube = 0` after the allocates, `MALLOC_PERTURB_=170` | byte-identical to its default run |
| branch, default / `MALLOC_PERTURB_=85` / `=170` | byte-identical to each other, 0.015 |
| branch vs master (default heap) | rms 7.6e-06, max abs 1.7e-04 |

| 20 x 32 | |
|---|---|
| master, 3 runs, `-finit-real` builds | as for 30 x 28 |
| master, `MALLOC_PERTURB_=170` | 0.463 vs truth (differs from master by rms 0.46, max abs 1.9) |
| branch, default / `MALLOC_PERTURB_=85` / `=170` | byte-identical to each other, 0.015 |
| branch vs master (default heap) | rms 1.3e-07, max abs 2.6e-06 |

- Master is deterministic with the default allocator (fresh pages are zero) and changes with other heap contents.
  `-finit-real` does not touch allocatable arrays, so it shows nothing.
- Part of the garbage is `std_cube`, but zeroing it alone does not make master stable: the padded part of `cube`
  (also never set on master, `reshape_up` writes only the data window) matters as well. The branch ignores both,
  since the padded pixels have noise 0 and are skipped.
- On `noise = .false.` the padded part of `cube` is still uninitialised (not changed by this branch).

(b) Blanked crop, blank pixels 64 (30 x 28) and 39 (20 x 32):

| crop | master blank / other | branch blank / other |
|---|---|---|
| 30 x 28 | 0.948 / 0.015 | 0.012 / 0.015 |
| 20 x 32 | 0.947 / 0.015 | 0.012 / 0.015 |

### 3. `descent = .false.` (blanked cube, only the bound change acts)

| run | blank pixels | other pixels |
|---|---|---|
| master, intact cube | 0.015 | 0.015 |
| master, blanked cube | 0.951 | 0.015 |
| branch, blanked cube | 0.013 | 0.015 |

- Intact cube with `rms.dat`: master and branch byte-identical.

### 4. On top of PR #9 (local ref pr9-tmp, contains #8)

Both commits cherry-pick onto pr9-tmp without conflicts. Built with `make OPENMP=1`, blanked cube, `OMP_NUM_THREADS` 1 and 4:

| run | blank pixels | other pixels |
|---|---|---|
| pr9-tmp alone, 1 thread = 4 threads (byte-identical) | 0.951 | 0.015 |
| pr9-tmp + fix, 1 thread = 4 threads (byte-identical) | 0.010 | 0.015 |

- Intact cube, `rms.dat`: pr9-tmp with and without the fix byte-identical.
- pr9-tmp + fix differs from the master-based branch (0.011) by rms 6.8e-04, max abs 1.9e-02, from the numerical changes in #8 and #9.
