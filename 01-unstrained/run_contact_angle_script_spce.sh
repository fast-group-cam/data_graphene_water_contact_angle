#!/bin/bash

echo "----- Run started at $(date '+%Y-%m-%d %H:%M:%S') -----"
echo ""
source ../00-common/software/venv/bin/activate

for folder in ./spce/*/run_prod; do
    for i in {0..9}; do
        if [ -f "${folder}/nvt_prod_${i}.lammpstrj" ]; then
            if [ ! -d "${folder}/angle_${i}" ]; then
                echo "[$(date '+%H:%M:%S')] Running ${folder}/angle_${i}..."
                contact_angle ${folder}/nvt_prod_${i}.lammpstrj -o ${folder}/angle_${i} --no-graphics > /dev/null
                echo "[$(date '+%H:%M:%S')] ...done"
                echo ""
            else
                echo "[$(date '+%H:%M:%S')] Skipped ${folder}/angle_${i} as it already exists"
                echo ""
            fi

            if [ ! -d "${folder}/angle-gibbs_${i}" ]; then
                echo "[$(date '+%H:%M:%S')] Running ${folder}/angle-gibbs_${i}..."
                python ../00-common/software/scripts/gibbs-contact-angle.py ${folder}/nvt_prod_${i}.lammpstrj -o ${folder}/angle-gibbs_${i} > /dev/null
                echo "[$(date '+%H:%M:%S')] ...done"
                echo ""
            else
                echo "[$(date '+%H:%M:%S')] Skipped ${folder}/angle-gibbs_${i} as it already exists"
                echo ""
            fi
        fi
    done
done

echo "----- Completed contact_angle scripts, consolidating results... -----"
echo ""

for folder in ./spce/*/run_prod; do
    python ../00-common/software/scripts/collate-blocks.py ${folder}/angle ${folder}/angle-gibbs -o ${folder}/contact-angle
    echo "[$(date '+%H:%M:%S')] Collated ${folder}/contact-angle"
done

deactivate
echo "----- Run completed at $(date '+%Y-%m-%d %H:%M:%S') -----"
echo ""
