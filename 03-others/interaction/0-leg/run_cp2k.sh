#!/usr/bin/env bash
set -euo pipefail

echo ""
ulimit -s unlimited

for root_folder in "s+0.00" "s+1.00" "s+2.00"; do
    # Use seq, then normalize to two decimals to match folder names like d+2.50
    for dist in $(seq 2.50 0.25 4.50); do
        dist_fmt=$(printf "%.2f" "$dist")
        dir="${root_folder}/d+${dist_fmt}"

        # Skip if the directory doesn't exist
        if [[ ! -d "$dir" ]]; then
            echo "Skipping missing directory: $dir"
            continue
        fi

        # If cp2k_sub.out exists, don't submit
        if [[ -e "${dir}/cp2k_sub.out" ]]; then
            echo "Already done: ${dir} (cp2k_sub.out found) — skipping sbatch"
            continue
        fi

        # Otherwise, submit
        pushd "$dir" >/dev/null
            cp ../../run.slurm .
            sbatch run.slurm
        popd >/dev/null
    done
done

echo "===== Sequence complete. ====="
echo ""

