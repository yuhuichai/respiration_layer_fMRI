#!/bin/bash
# replace xxx with your own data and subject directory
# set data directories
top_dir=/media/yuhui/xxx

cd ${top_dir}
for patDir in xxx/mt.sft; do
{
	
	echo "*********************************************************************************************"
	echo " start working in ${patDir} "
	echo "*********************************************************************************************"

	sumaDir=${top_dir}/${patDir}/SUMA
	cd ${sumaDir}
	3dcalc -a aparc+aseg_REN_gm.nii.gz -expr "step(a-45)*1000" \
		-prefix gm_boosted.nii.gz -overwrite
	3dcalc -a aparc+aseg_REN_wmat.nii.gz -expr "amongst(a,1,21,38,41,42,43,44,45)*1000" \
		-prefix wm_boosted.nii.gz -overwrite
	3dcalc -a aparc+aseg_REN_all.nii.gz -expr "step(a)*1000" \
		-prefix all_boosted.nii.gz -overwrite

	cd ${top_dir}/${patDir}
	template=mean.sub_d_dant.beta100.nii.gz
	delta_x=$(3dinfo -adi $template)
	delta_y=$(3dinfo -adj $template)
	delta_z=$(3dinfo -adk $template)
	sdelta_x=$(echo "(($delta_x / 4))"|bc -l)
	sdelta_y=$(echo "(($delta_x / 4))"|bc -l)
	sdelta_z=$(echo "(($delta_z / 4))"|bc -l) # here I only upscale in 2 dimensions.
	3dresample -dxyz $sdelta_x $sdelta_y $sdelta_z -rmode NN -overwrite -prefix scaledXYZ_$template -input $template
	# 3dresample -dxyz $sdelta_x $sdelta_y $delta_z -rmode NN -overwrite -prefix scaledXY_$template -input $template

	cd ${sumaDir}

	3dresample -master ../scaledXYZ_mean.sub_d_dant.beta100.nii.gz -rmode NN \
		-overwrite -prefix gm_upsamp.nii.gz -input gm_boosted.nii.gz

	3dresample -master ../scaledXYZ_mean.sub_d_dant.beta100.nii.gz -rmode NN \
		-overwrite -prefix wm_upsamp.nii.gz -input wm_boosted.nii.gz


	3dresample -master ../scaledXYZ_mean.sub_d_dant.beta100.nii.gz -rmode NN \
		-overwrite -prefix all_upsamp.nii.gz -input all_boosted.nii.gz


	3dcalc -a wm_upsamp.nii.gz -expr "step(a-500)" \
		-prefix wm_upsamp_thr.nii.gz -overwrite
	3dcalc -a all_upsamp.nii.gz -expr "step(a-500)" \
		-prefix all_upsamp_thr.nii.gz -overwrite
	3dcalc -a gm_upsamp.nii.gz -expr "step(a-200)" \
		-prefix gm_upsamp_thr.nii.gz -overwrite


	rm gm_upsamp.nii.gz wm_upsamp.nii.gz subc_upsamp.nii.gz


	3dmask_tool -dilate_inputs 1 -prefix gm_upsamp_thr_d1.nii.gz -overwrite -inputs gm_upsamp_thr.nii.gz

	# extract white matter edge
	3dmask_tool -dilate_inputs -1 -prefix wm_upsamp_thr_e1.nii.gz -overwrite -inputs wm_upsamp_thr.nii.gz
	3dcalc -a wm_upsamp_thr.nii.gz -b wm_upsamp_thr_e1.nii.gz  -c gm_upsamp_thr_d1.nii.gz -expr '(step(a)-step(b))*step(c)' \
		-prefix wm_upsamp_edge.nii.gz -overwrite

	# extract csf edge
	3dmask_tool -dilate_inputs 1 -prefix all_upsamp_thr_d1.nii.gz -overwrite -inputs all_upsamp_thr.nii.gz
	3dcalc -a all_upsamp_thr.nii.gz -b all_upsamp_thr_d1.nii.gz -c gm_upsamp_thr_d1.nii.gz -expr '(step(b)-step(a))*step(c)' \
		-prefix all_upsamp_edge.nii.gz -overwrite

	# generate cortical layermask
	3dcalc -a all_upsamp_edge.nii.gz -b wm_upsamp_edge.nii.gz -c gm_upsamp_thr.nii.gz \
		-expr "(step(a)+2*step(b)+3*step(c)*iszero(a+b))" \
		-prefix LayerMask4smooth.nii.gz -overwrite

	# generate csf layermask
	3dmask_tool -dilate_inputs 1 -prefix all_upsamp_thr_d1.nii.gz -overwrite -inputs all_upsamp_thr.nii.gz
	3dmask_tool -dilate_inputs 2 -prefix all_upsamp_thr_d2.nii.gz -overwrite -inputs all_upsamp_thr.nii.gz
	3dcalc -a all_upsamp_thr_d2.nii.gz -b all_upsamp_thr_d1.nii.gz -expr "(step(a)-step(b))" \
		-prefix csfMask_edge2.nii.gz -overwrite

	3dmask_tool -dilate_inputs 5 -prefix all_upsamp_thr_d5.nii.gz -overwrite -inputs all_upsamp_thr.nii.gz
	3dcalc -a all_upsamp_thr_d5.nii.gz -b all_upsamp_thr_d2.nii.gz -expr "(step(a)-step(b))" \
		-prefix csfMask.nii.gz -overwrite

	3dmask_tool -dilate_inputs 6 -prefix all_upsamp_thr_d6.nii.gz -overwrite -inputs all_upsamp_thr.nii.gz
	3dcalc -a all_upsamp_thr_d6.nii.gz -b all_upsamp_thr_d5.nii.gz -expr "(step(a)-step(b))" \
		-prefix csfMask_edge1.nii.gz -overwrite

	3dcalc -a csfMask_edge1.nii.gz -b csfMask_edge2.nii.gz -c csfMask.nii.gz -d LayerMask4smooth.nii.gz \
		-expr "(step(a)+2*step(b)+3*step(c)*iszero(a+b))*iszero(d)" \
		-prefix csfMask4smooth.nii.gz -overwrite

	# generate wm layermask
	3dmask_tool -dilate_inputs -1 -prefix wm_upsamp_thr_e1.nii.gz -overwrite -inputs wm_upsamp_thr.nii.gz
	3dmask_tool -dilate_inputs -2 -prefix wm_upsamp_thr_e2.nii.gz -overwrite -inputs wm_upsamp_thr.nii.gz
	3dmask_tool -dilate_inputs -4 -prefix wm_upsamp_thr_e4.nii.gz -overwrite -inputs wm_upsamp_thr.nii.gz
	3dmask_tool -dilate_inputs -5 -prefix wm_upsamp_thr_e5.nii.gz -overwrite -inputs wm_upsamp_thr.nii.gz

	3dmask_tool -dilate_inputs 5 -prefix gm_upsamp_thr_d5.nii.gz -overwrite -inputs gm_upsamp_thr.nii.gz

	3dcalc -a wm_upsamp_thr_e4.nii.gz -b wm_upsamp_thr_e2.nii.gz -expr "(step(b)-step(a))" \
		-prefix wmMask.nii.gz -overwrite

	3dcalc -a wm_upsamp_thr_e5.nii.gz -b wm_upsamp_thr_e4.nii.gz -expr "(step(b)-step(a))" \
		-prefix wmMask_edge2.nii.gz -overwrite

	3dcalc -a wm_upsamp_thr_e2.nii.gz -b wm_upsamp_thr_e1.nii.gz -expr "(step(b)-step(a))" \
		-prefix wmMask_edge1.nii.gz -overwrite

	3dcalc -a wmMask_edge1.nii.gz -b wmMask_edge2.nii.gz -c wmMask.nii.gz -d LayerMask4smooth.nii.gz -e gm_upsamp_thr_d5.nii.gz \
		-expr "(step(a)+2*step(b)+3*step(c)*iszero(a+b))*iszero(d)*step(e)" \
		-prefix wmMask4smooth.nii.gz -overwrite

	rm *thr_d* *thr_e* csfMask.nii.gz wmMask.nii.gz *edge*.nii.gz *upsamp.nii.gz *upsamp_thr.nii.gz

	echo "************** grow cortical layers with RENZO's program ******************"
	if [ ! -f layers3D_4smooth_FR_curvature.nii.gz ]; then
		# -curvature -streamlines -thickness 
		LN2_LAYERS -rim LayerMask4smooth.nii.gz -nr_layers 18 \
			-equivol -incl_borders -equal_counts -curvature \
			-output layers3D_4smooth_FR &
		# layerbins equal_counts

		echo "************** grow wm and csf layers with RENZO's program ******************"
		LN2_LAYERS -rim csfMask4smooth.nii.gz -nr_layers 6 \
			-equivol -incl_borders -equal_counts \
			-output csf3D_4smooth_FR &

		LN2_LAYERS -rim wmMask4smooth.nii.gz -nr_layers 5 \
			-equivol -incl_borders -equal_counts \
			-output wm3D_4smooth_FR &

		wait

		rm *metric* # metric no_smooth

		gzip -v -f layers3D_4smooth_FR*nii wm3D_4smooth_FR*nii csf3D_4smooth_FR*nii
	fi

	echo "************** combine cortical layers, wm and csf ******************"
	3dresample -master ../scaledXYZ_mean.sub_d_dant.beta100.nii.gz -rmode NN \
		-overwrite -prefix scaledXYZ_brain_mask_comb.nii.gz -input ../../brain_mask_comb.nii
	3dcalc -a wm3D_4smooth_FR_layerbins_equivol.nii.gz \
		-b layers3D_4smooth_FR_layerbins_equivol.nii.gz \
		-c csf3D_4smooth_FR_layerbins_equivol.nii.gz \
		-d scaledXYZ_brain_mask_comb.nii.gz \
		-expr "(a+(b+step(b)*5)+(c+step(c)*23))*step(d)" \
		-prefix all3D_4smooth_FR_layerbins_equivol.nii.gz -overwrite

	
}&
done
wait

