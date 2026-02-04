#!/bin/bash

echo "----- Run started at $(date '+%Y-%m-%d %H:%M:%S') -----"
echo ""
source ../00-common/software/venv/bin/activate

for folder in ./mace/*/run_prod; do
    if [ ! -d "${folder}/contact-angle" ]; then
        for i in {0..9}; do
            a="${folder}/nvt_prod_$(( 3*i )).lammpstrj"
            b="${folder}/nvt_prod_$(( 3*i + 1 )).lammpstrj"
            c="${folder}/nvt_prod_$(( 3*i + 2 )).lammpstrj"
            if [ -f ${a} ] && [ -f ${b} ] && [ -f ${c} ] && [ ! -d "${folder}/angle_${i}" ]; then
                echo "[$(date '+%H:%M:%S')] Running ${folder}/angle_${i}..."
                contact_angle ${a} ${b} ${c} -o ${folder}/angle_${i} > /dev/null
                echo "[$(date '+%H:%M:%S')] ...done"
                echo ""
            fi
            if [ -f ${a} ] && [ -f ${b} ] && [ -f ${c} ] && [[ ${folder} == *"fixed-"* ]] && [ ! -d "${folder}/angle-gibbs_${i}" ]; then
                echo "[$(date '+%H:%M:%S')] Running ${folder}/angle-gibbs_${i}..."
                python ../00-common/software/scripts/gibbs-contact-angle.py ${a} ${b} ${c} -o ${folder}/angle-gibbs_${i} > /dev/null
                echo "[$(date '+%H:%M:%S')] ...done"
                echo ""
            fi
        done
    else
        echo "[$(date '+%H:%M:%S')] Skipped ${folder} as contact-angle folder already exists"
        echo ""
    fi
done

echo "----- Completed contact_angle scripts, consolidating results... -----"
echo ""

for folder in ./mace/*/run_prod; do
    if [ ! -d "${folder}/contact-angle" ]; then
        if [[ ${folder} == *"fixed-"* ]]; then
            python ../00-common/software/scripts/collate-blocks.py ${folder}/angle ${folder}/angle-gibbs -o ${folder}/contact-angle
        else
            python ../00-common/software/scripts/collate-blocks.py ${folder}/angle -o ${folder}/contact-angle
        fi
        echo "[$(date '+%H:%M:%S')] Collated ${folder}/contact-angle"
    fi
done

deactivate
echo "----- Run completed at $(date '+%Y-%m-%d %H:%M:%S') -----"
echo ""
