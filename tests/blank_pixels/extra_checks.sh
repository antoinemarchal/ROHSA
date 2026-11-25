#!/bin/bash
# extra_checks.sh <master binary> <branch binary> <work dir>
# Reproduces the sections of RESULTS.md after the first table (lym, padded grids, descent off). Needs numpy.
# Optional: MASTER_SNAN, MASTER_ZERO = master builds with -finit-real=snan / zero, MASTER_STDZERO = master with only
# `std_cube = 0._xp` added after its allocate, MASTER_ALLZERO = master with `std_cube = 0._xp` and `cube = 0._xp` added
# (default: skipped). The -finit-real builds do not touch allocatable arrays, so they do not test this.
set -e
M=$1; B=$2; W=$3
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$W"
py="python3 $here/make_synthetic.py"
$py "$W/full"
$py "$W/c30x28" 30x28
$py "$W/c20x32" 20x32
same() { cmp -s "$1" "$2" && echo identical || echo DIFFERENT; }

echo "== 1. lym = .true., lambda_lym_sig = 10 (32 x 32)"
d=$W/full
export LYM=.true. LLYM=10.0d0
$here/run.sh $M $d lym_m_int cube.dat rms.dat;        $here/run.sh $B $d lym_b_int cube.dat rms.dat
$here/run.sh $M $d lym_m_blank cube_blank.dat rms_blank.dat; $here/run.sh $B $d lym_b_blank cube_blank.dat rms_blank.dat
echo "intact cube, noise map without zeros, master vs branch: $(same $d/runs/lym_m_int.out $d/runs/lym_b_int.out)"
python3 $here/compare.py $d lym_m_int lym_b_int lym_m_blank lym_b_blank
unset LYM LLYM

echo "== 2. padded grids"
for c in c30x28 c20x32; do
  d=$W/$c; echo "-- $c"
  $here/run.sh $M $d m1 cube.dat rms.dat; $here/run.sh $M $d m2 cube.dat rms.dat; $here/run.sh $M $d m3 cube.dat rms.dat
  MALLOC_PERTURB_=85 $here/run.sh $M $d m_pert85 cube.dat rms.dat
  MALLOC_PERTURB_=170 $here/run.sh $M $d m_pert170 cube.dat rms.dat
  [ -n "$MASTER_SNAN" ] && $here/run.sh $MASTER_SNAN $d m_snan cube.dat rms.dat
  [ -n "$MASTER_ZERO" ] && $here/run.sh $MASTER_ZERO $d m_zero cube.dat rms.dat
  [ -n "$MASTER_STDZERO" ] && { $here/run.sh $MASTER_STDZERO $d m_stdzero cube.dat rms.dat; MALLOC_PERTURB_=170 $here/run.sh $MASTER_STDZERO $d m_stdzero170 cube.dat rms.dat; }
  [ -n "$MASTER_ALLZERO" ] && { $here/run.sh $MASTER_ALLZERO $d m_allzero cube.dat rms.dat; MALLOC_PERTURB_=170 $here/run.sh $MASTER_ALLZERO $d m_allzero170 cube.dat rms.dat; }
  $here/run.sh $B $d b cube.dat rms.dat
  MALLOC_PERTURB_=85 $here/run.sh $B $d b_pert85 cube.dat rms.dat; MALLOC_PERTURB_=170 $here/run.sh $B $d b_pert170 cube.dat rms.dat
  $here/run.sh $M $d mb cube_blank.dat rms_blank.dat; $here/run.sh $B $d bb cube_blank.dat rms_blank.dat
  echo "master run 1 vs run 2: $(same $d/runs/m1.out $d/runs/m2.out); vs run 3: $(same $d/runs/m1.out $d/runs/m3.out)"
  echo "branch vs b_pert85: $(same $d/runs/b.out $d/runs/b_pert85.out); vs b_pert170: $(same $d/runs/b.out $d/runs/b_pert170.out)"
  [ -n "$MASTER_STDZERO" ] && echo "master+std_cube=0, default vs MALLOC_PERTURB_=170: $(same $d/runs/m_stdzero.out $d/runs/m_stdzero170.out)"
  [ -n "$MASTER_ALLZERO" ] && echo "master+std_cube=0+cube=0, default vs MALLOC_PERTURB_=170: $(same $d/runs/m_allzero.out $d/runs/m_allzero170.out)"
  for n in m_pert85 m_pert170 m_snan m_zero; do [ -f $d/runs/$n.out ] && echo "master run 1 vs $n: $(same $d/runs/m1.out $d/runs/$n.out)"; done
  python3 $here/compare.py $d m1 m_pert85 m_pert170 b b_pert170 mb bb
  [ -n "$MASTER_ALLZERO" ] && echo "master+std_cube=0+cube=0, default vs MALLOC_PERTURB_=170: $(same $d/runs/m_allzero.out $d/runs/m_allzero170.out)"
  for n in m_pert85 m_pert170 m_snan m_zero m_stdzero m_stdzero170 m_allzero m_allzero170 b; do [ -f $d/runs/$n.out ] && python3 $here/compare.py $d --diff m1 $n; done
  python3 $here/compare.py $d --diff m1 m2
done

echo "== 3. descent = .false. (32 x 32, blanked cube)"
d=$W/full
export DESCENT=.false.
$here/run.sh $M $d nd_m_int cube.dat rms.dat; $here/run.sh $M $d nd_m_blank cube_blank.dat rms_blank.dat; $here/run.sh $B $d nd_b_blank cube_blank.dat rms_blank.dat
$here/run.sh $B $d nd_b_int cube.dat rms.dat
echo "intact cube, master vs branch: $(same $d/runs/nd_m_int.out $d/runs/nd_b_int.out)"
python3 $here/compare.py $d nd_m_int nd_m_blank nd_b_blank
