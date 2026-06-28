#!/usr/bin/env bash
# PostToolUse hook: when an edited file is a bin/<engine>.py that ships a
# `selftest` subcommand, run that selftest so a regression surfaces in-loop
# (not at CI). Exits 2 on failure so the result is fed back to Claude; exits 0
# (silently) for any edit that is not a tracked engine script.
input=$(cat)
fp=$(printf '%s' "$input" | python3 -c 'import sys,json; print(json.load(sys.stdin).get("tool_input",{}).get("file_path",""))' 2>/dev/null)
[[ -n $fp ]] || exit 0
base=$(basename "$fp")
case "$base" in
  panel_select.py|interrogate.py|mrd_integrate.py|validate.py|build_background.py|pyclone_prep.py|normal_evidence.py|sample_id.py) ;;
  *) exit 0 ;;
esac
repo=$(cd "$(dirname "$fp")/.." 2>/dev/null && pwd)
[[ -f $repo/bin/$base ]] || exit 0
# prefer the project env python (has numpy/scipy/pysam); fall back to system
py=python3
[[ -x /data/alvin/envs/mrd/bin/python3 ]] && py=/data/alvin/envs/mrd/bin/python3
if out=$("$py" "$repo/bin/$base" selftest 2>&1); then
  echo "selftest OK: $base"
  exit 0
fi
{ echo "SELFTEST FAILED for $base after edit:"; echo "$out"; } >&2
exit 2
