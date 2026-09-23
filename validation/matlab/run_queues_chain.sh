#!/usr/bin/env bash
# Run several job files one after another with the same parallelism.
#   validation/matlab/run_queues_chain.sh <max_parallel> <jobfile1> [jobfile2 ...]
cd "$(dirname "$0")/../.."
maxp=$1; shift
for jf in "$@"; do
  echo "$(date +%T) === queue $jf"
  bash validation/matlab/run_queue.sh "$jf" "$maxp"
done
echo "$(date +%T) === all queues finished"
