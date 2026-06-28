#!/usr/bin/env bash
# Run every bin/*.py selftest. Runs all (does not stop at first), reports a
# summary, exits non-zero if any failed. Used by CI and locally.
#   PATH=/data/alvin/envs/mrd/bin:$PATH bash bin/run_selftests.sh
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
SCRIPTS=(panel_select interrogate mrd_integrate validate build_background
         pyclone_prep normal_evidence sample_id)
fail=0
for s in "${SCRIPTS[@]}"; do
  printf '=== %s ===\n' "$s"
  if python3 "bin/$s.py" selftest; then echo "PASS $s"; else echo "FAIL $s"; fail=1; fi
done
echo "----"
if [[ $fail -eq 0 ]]; then echo "ALL ${#SCRIPTS[@]} SELFTESTS PASS"; else echo "SELFTEST FAILURES"; exit 1; fi
