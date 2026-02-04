#!/bin/bash

echo "----- Run started at $(date '+%Y-%m-%d %H:%M:%S') -----"
echo ""
source ../00-common/software/venv/bin/activate

for FOLDER in s*0**/run_prod; do
    if [ ! -d "${FOLDER}/contact-angle" ]; then
        for i in {0..9}; do
            a="${FOLDER}/nvt_prod_$(( 3*i )).lammpstrj"
            b="${FOLDER}/nvt_prod_$(( 3*i + 1 )).lammpstrj"
            c="${FOLDER}/nvt_prod_$(( 3*i + 2 )).lammpstrj"
            if [ -f ${a} ] && [ -f ${b} ] && [ -f ${c} ] && [ ! -d "${FOLDER}/angle_${i}" ]; then
                echo "[$(date '+%H:%M:%S')] Running ${FOLDER}/angle_${i}..."
                contact_angle ${a} ${b} ${c} -o ${FOLDER}/angle_${i} > /dev/null
                echo "[$(date '+%H:%M:%S')] ...done"
                echo ""
            fi
        done
    else
        echo "[$(date '+%H:%M:%S')] Skipped ${FOLDER} as contact-angle folder already exists"
        echo ""
    fi

    if [ ! -d "${FOLDER}/contact-angle" ]; then
        python ../00-common/software/scripts/collate-blocks.py ${FOLDER}/angle -o ${FOLDER}/contact-angle
        echo "[$(date '+%H:%M:%S')] Collated ${FOLDER}/contact-angle"
        echo ""
    fi

    targets=""
    for i in {0..99}; do
        if [ -f "${FOLDER}/nvt_prod_${i}.lammpstrj" ]; then
            targets="${targets} ${FOLDER}/nvt_prod_${i}.lammpstrj"
        fi
    done
    if [ ! -d "${FOLDER}/contact-angle-aniso" ]; then
        echo "[$(date '+%H:%M:%S')] Running ${FOLDER}/contact-angle-aniso..."
        contact_angle_anisotropic ${targets} -o ${FOLDER}/contact-angle-aniso > /dev/null
        echo "[$(date '+%H:%M:%S')] ...done"
    else
        echo "[$(date '+%H:%M:%S')] Skipped ${FOLDER} as contact-angle-aniso already exists"
    fi
    echo ""

done

deactivate
echo "----- Run completed at $(date '+%Y-%m-%d %H:%M:%S') -----"
echo ""
