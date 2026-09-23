#!/usr/bin/env bash
# run_queue.sh -- run Phase A fits as independent `matlab -batch` processes, N at a time.
#
#   validation/matlab/run_queue.sh <jobfile> <max_parallel>
#
# <jobfile>: one job per line:  <btc_id> <objective> <t_end_mode> <seed> [instrumented]
# Lines starting with # are skipped. A job whose output folder already contains
# summary.json is skipped (so the queue can be re-run after interruption).
# Logs: private/validation/matlab_out/queue_logs/<btc>__<objective>__<tend>__seed<seed>.log
set -u
cd "$(dirname "$0")/../.."
jobfile=$1; maxp=${2:-4}
logdir=private/validation/matlab_out/queue_logs; mkdir -p "$logdir"

run_job() {
  local btc=$1 obj=$2 tend=$3 seed=$4 instr=${5:-false}
  local name="${btc}__${obj}__${tend}__seed${seed}"
  if [ -f "private/validation/matlab_out/${name}/summary.json" ]; then
    echo "$(date +%T) skip (done): $name"; return 0
  fi
  echo "$(date +%T) start: $name (instrumented=$instr)"
  matlab -batch "addpath('validation/matlab'); run_phase_a('${btc}', 'objective', '${obj}', 't_end_mode', '${tend}', 'seed', ${seed}, 'instrumented', ${instr});" \
    > "$logdir/${name}.log" 2>&1
  echo "$(date +%T) end ($?): $name"
}

while read -r btc obj tend seed instr; do
  btc=${btc//$'\r'/}; obj=${obj//$'\r'/}; tend=${tend//$'\r'/}; seed=${seed//$'\r'/}; instr=${instr//$'\r'/}
  [[ -z "$btc" || "$btc" == \#* ]] && continue
  while [ "$(jobs -rp | wc -l)" -ge "$maxp" ]; do sleep 30; done
  run_job "$btc" "$obj" "$tend" "$seed" "${instr:-false}" &
  sleep 5
done < "$jobfile"
wait
echo "$(date +%T) queue finished"
