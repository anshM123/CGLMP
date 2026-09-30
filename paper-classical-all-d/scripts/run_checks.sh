#!/usr/bin/env bash
# Re-runs every numerical cross-check of the paper and writes the outputs to ../logs/.
# Usage (from this folder):  bash run_checks.sh        [set PYTHON=... to choose the interpreter]
set -e
PY="${PYTHON:-python}"
export PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"
mkdir -p ../logs
env_line() { "$PY" -c "import sys, numpy, scipy, mpmath; print('# python', sys.version.split()[0], '| numpy', numpy.__version__, '| scipy', scipy.__version__, '| mpmath', mpmath.__version__)"; }
run() {  # run <logname> <script> <args...>
  local log="../logs/$1.log"; shift
  { echo "# command: python $*"; env_line; } > "$log"
  local t0=$(date +%s)
  "$PY" "$@" >> "$log" 2>&1
  echo "# wall time: $(( $(date +%s) - t0 )) s" >> "$log"
  echo "done: $log"
}
run verify_writeup   verify_writeup.py 9
run verify_junction  verify_junction.py 6 7 8 12 16 24 40 64 101
run verify_chain     verify_chain.py 24,6,6,500 40,10,10,500 48,12,12,300 80,20,20,100 30,4,11,500 50,7,19,300
run cot_rearr        cot_rearr.py 12,3,3 16,4,4 20,5,5 24,6,6 16,3,6 20,4,7
run swap_test        swap_test.py 9,3,3 10,3,3 12,3,4 12,4,4
run classical_sa     classical_sa.py 12
run check_reduction  check_reduction.py
run check_constants  check_constants.py
