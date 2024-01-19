# !/bin/bash

# set data directories
dataDir=/media/yuhui/xxx # replace with your own data directory
cd ${dataDir}

for patDir in xxx; do # replace with your own patient name
{

	echo "*****************************************************************************************************"
	echo ++ ${patDir}
	echo "*****************************************************************************************************"

	cd ${dataDir}/${patDir}
	for runDir in vaper.sft; do 
	{	
		if [ -d ${dataDir}/${patDir}/${runDir} ]; then
			cd ${dataDir}/${patDir}/${runDir}

			run_dsets=($(ls -f rbold*.nii.gz))
			run_num=${#run_dsets[@]}

			echo "******************** outcount censor ***********************"
			deepp=0.0
			deeppdiff=1.0

			while (( $(echo $deepp'=='0.0 | bc -l) )) || (( $(echo $deeppdiff'>'0.015 | bc -l) )) || (( $(echo $deeppdiff'<'-0.015 | bc -l) )); do 
				if (( $(echo $deepp'=='0.0 | bc -l) )); then
					outc=0.006
				else
					if (( $(echo $deeppdiff'>'0.015 | bc -l) )); then
						outc=`bc -l <<< "${outc}+0.0001"`
					elif (( $(echo $deeppdiff'<'-0.015 | bc -l) )); then
						outc=`bc -l <<< "${outc}-0.0001"`
					fi
				fi


				for run in `seq 1 ${run_num}`; do
				{
				    if [ ! -f outcount.bold_mdant.r$run.1D ]; then 
					    3dToutcount -mask brain_mask.nii.gz -fraction -polort 4 -legendre \
					                bold_mdant${run}.nii.gz > outcount.bold_mdant.r$run.1D
					fi

					if [ ! -f outcount.sub_d_bold.r$run.1D ]; then
					    3dToutcount -mask brain_mask.nii.gz -fraction -polort 4 -legendre \
					                sub_d_bold${run}.nii.gz > outcount.sub_d_bold.r$run.1D
					fi

				    1deval -a outcount.bold_mdant.r$run.1D -expr "1-step(a-${outc})" > rm.out.cen.bold_mdant.r$run.1D
				    1deval -a outcount.sub_d_bold.r$run.1D -expr "1-step(a-${outc})" > rm.out.cen.sub_d_bold.r$run.1D
				}&
				done
				wait
				# catenate outlier censor files into a single time series
				cat rm.out.cen.bold_mdant.r*.1D > outcount_bold_mdant_censor.1D
				cat rm.out.cen.sub_d_bold.r*.1D > outcount_sub_d_bold_censor.1D
				rm rm.out.cen*

				1deval -overwrite -a motion_bold_censor.1D -b motion_dant_censor.1D \
					-expr "a*b" > censor_motion.1D

				1deval -overwrite -a censor_motion.1D -b outcount_bold_mdant_censor.1D -c outcount_sub_d_bold_censor.1D \
				       -expr "a*step((1-b)+(1-c))" > outcount_deepbreath_lowmotion.1D

				1deval -overwrite -a censor_motion.1D -b outcount_bold_mdant_censor.1D -c outcount_sub_d_bold_censor.1D \
				       -expr "a*(1-step((1-b)+(1-c)))" > outcount_normbreath_lowmotion.1D

				ktrs_outcount_deep=`1d_tool.py -infile outcount_deepbreath_lowmotion.1D -show_trs_uncensored encoded`
				ktrs_outcount_norm=`1d_tool.py -infile outcount_normbreath_lowmotion.1D -show_trs_uncensored encoded`

				deepn=`1dsum outcount_deepbreath_lowmotion.1D`
				normn=`1dsum outcount_normbreath_lowmotion.1D`
				deepp=`bc -l <<< "${deepn}/(${normn}+${deepn})"`
				deeppdiff=`bc -l <<< "${deepp}-0.15"`

				echo "++ outc = ${outc}, deepn = ${deepn}, normn = ${normn}, deepp = ${deepp}, deeppdiff = ${deeppdiff}"
			done

			echo "++ outcount_deepbreath_lowmotion: $ktrs_outcount_deep ..."
			echo "++ outcount_normbreath_lowmotion: $ktrs_outcount_norm ..."

			subjList='rbold bold_mdant'

			for subj in $subjList; do
			{
				if [ ! -f all_runs.${subj}.nii.gz ]; then
					3dTcat -overwrite -prefix all_runs.${subj}.nii.gz ${subj}*nii.gz -overwrite
				fi

				# 3dTcat -overwrite -prefix deepbreath_biopac.${subj}.nii.gz all_runs.${subj}.nii.gz"[$ktrs_biopac_deep]" &
				# 3dTcat -overwrite -prefix normbreath_biopac.${subj}.nii.gz all_runs.${subj}.nii.gz"[$ktrs_biopac_norm]" &

				3dTcat -overwrite -prefix deepbreath_outcount.${subj}.nii.gz all_runs.${subj}.nii.gz"[$ktrs_outcount_deep]" &
				3dTcat -overwrite -prefix normbreath_outcount.${subj}.nii.gz all_runs.${subj}.nii.gz"[$ktrs_outcount_norm]" &
				wait
				sleep 2

				for way in outcount; do
				{
					3dTstat -overwrite -mean -prefix mean.normbreath_${way}.${subj}.nii.gz normbreath_${way}.${subj}.nii.gz
					
					3dcalc -a deepbreath_${way}.${subj}.nii.gz -b mean.normbreath_${way}.${subj}.nii.gz -c brain_mask.nii.gz \
						-expr "abs(a-b)*step(c)" -prefix vessels_${way}.${subj}.nii.gz -overwrite
					3dTstat -overwrite -mean -prefix mean.vessels_${way}.${subj}.nii.gz vessels_${way}.${subj}.nii.gz

					3dcalc -a deepbreath_${way}.${subj}.nii.gz -b mean.normbreath_${way}.${subj}.nii.gz \
						-expr "abs(a-b)" -prefix vessels_${way}.${subj}.nomask.nii.gz -overwrite
					3dTstat -overwrite -mean -prefix mean.vessels_${way}.${subj}.nomask.nii.gz vessels_${way}.${subj}.nomask.nii.gz
				}&
				done
				wait

				rm *normbreath_*.${subj}.nii.gz deepbreath_*.${subj}.nii.gz 
				rm vessels_*.${subj}*.nii.gz mean.normbreath_${way}.${subj}.nii.gz

				rm all_runs.${subj}.nii.gz
			}&
			done
			wait

			ratio_vessel_goal=0.1

			3dresample -master mean.vessels_outcount.bold_mdant.nii.gz -rmode NN \
				-overwrite -prefix rm.LayerMask_dnsamp.nii.gz -input ../mt.sft/SUMA/LayerMask4smooth.nii.gz

			
			3dcalc -a rm.LayerMask_dnsamp.nii.gz -expr "step(a)" -prefix LayerMask_dnsamp.nii.gz -overwrite

			ratio_vessel=0.0

			while (( $(echo $ratio_vessel'=='0.0 | bc -l) )) || (( $(echo $vesseldiff'>'0.005 | bc -l) )) || (( $(echo $vesseldiff'<'-0.005 | bc -l) )); do 
				if (( $(echo $ratio_vessel'=='0.0 | bc -l) )); then
					vesselthr=55.0
				else
					if (( $(echo $vesseldiff'>'0.005 | bc -l) )); then
						vesselthr=`bc -l <<< "${vesselthr}+0.5"`
					elif (( $(echo $vesseldiff'<'-0.005 | bc -l) )); then
						vesselthr=`bc -l <<< "${vesselthr}-0.5"`
					fi
				fi

				3dcalc -a mean.vessels_outcount.bold_mdant.nii.gz -expr "step(a-${vesselthr})*a" \
					-prefix meanthr.vessels_outcount.bold_mdant.nii.gz -overwrite

				3dcalc -a meanthr.vessels_outcount.bold_mdant.nii.gz -b LayerMask_dnsamp.nii.gz \
					-expr "step(a)*step(b)" -prefix mask_vessels.nii.gz -overwrite

				nvxl_vessel=`3dROIstats -nomeanout -nzvoxels -quiet -mask mask_vessels.nii.gz mask_vessels.nii.gz`

				nvxl_brain=`3dROIstats -nomeanout -nzvoxels -quiet -mask LayerMask_dnsamp.nii.gz LayerMask_dnsamp.nii.gz`

				ratio_vessel=`bc -l <<< "${nvxl_vessel}/${nvxl_brain}"`

				vesseldiff=`bc -l <<< "${ratio_vessel}-${ratio_vessel_goal}"`

				echo "++ vesselthr = ${vesselthr}, ratio_vessel = ${ratio_vessel} ..."

			done

			rm mask_vessels.nii.gz


		fi
			
	}&
	done
	wait


}&
done
wait

