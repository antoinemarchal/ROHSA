#!/bin/bash
# run.sh <ROHSA binary> <data dir> <name> <cube> <rms file | none>
# Optional environment: DESCENT=.false. (default .true.), LYM=.true. (default .false.), LLYM (lambda_lym_sig, default 0).
# Runs one fit in <data dir>/runs/<name>.out. "none" runs with noise = .false. (no noise map).
set -e
bin=$1; data=$2; name=$3; cube=$4; rms=$5
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$data/runs"
if [ "$rms" = none ]; then noise=.false.; rmsfile=''; else noise=.true.; rmsfile="$data/$rms"; fi
sed -e "s#@CUBE@#$data/$cube#" -e "s#@RMS@#$rmsfile#" -e "s#@OUT@#$data/runs/$name#" -e "s#@NOISE@#$noise#" \
    -e "s#@DESCENT@#${DESCENT:-.true.}#" -e "s#@LYM@#${LYM:-.false.}#" -e "s#@LLYM@#${LLYM:-0.0d0}#" \
    "$here/params_template.txt" > "$data/runs/$name.params"
"$bin" "$data/runs/$name.params" > "$data/runs/$name.log"
