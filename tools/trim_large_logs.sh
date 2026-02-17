#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob globstar

# Target safely below GitHub's 100MB blob limit
MAX_BYTES="${MAX_BYTES:-95000000}"

# Keep 20% from start + 20% from end (by lines)
KEEP_FRACTION="${KEEP_FRACTION:-0.20}"

PATTERNS=(
  "case_*/slurm/slurm.*.out"
  "case_*/postProcessing/**/foamRun.stdout.log"
)

trim_one() {
  local f="$1"
  [[ -f "$f" ]] || return 0

  local sz
  sz="$(stat -c%s "$f")"
  (( sz > MAX_BYTES )) || return 0

  echo "Trimming: $f (current $(numfmt --to=iec --suffix=B "$sz"))"

  # Keep a local backup of the full file (not tracked unless you add it)
  if [[ ! -f "${f}.full" ]]; then
    cp -p "$f" "${f}.full"
  fi

  local total keep
  total="$(wc -l < "$f" | tr -d ' ')"
  keep="$(awk -v t="$total" -v frac="$KEEP_FRACTION" 'BEGIN{ k=int(t*frac); if(k<1) k=1; print k }')"

  local tmp
  tmp="$(mktemp "${f}.trim.XXXXXX")"
  {
    head -n "$keep" "$f" || true
    echo "----- TRIMMED (${KEEP_FRACTION} head + ${KEEP_FRACTION} tail) original_lines=${total} kept_lines_per_side=${keep} -----"
    tail -n "$keep" "$f" || true
  } > "$tmp"

  local tsize
  tsize="$(stat -c%s "$tmp")"

  if (( tsize < sz )); then
    mv "$tmp" "$f"
    echo " -> new size $(numfmt --to=iec --suffix=B "$tsize")"
  else
    rm -f "$tmp"
    echo " -> no improvement; left unchanged"
  fi
}

for pat in "${PATTERNS[@]}"; do
  for f in $pat; do
    trim_one "$f"
  done
done
