%% replace xxx with your own data and subject directory
clear;
close all;
dataDir='/media/yuhui/xxx';


%%%%% vein segmenation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(dataDir);
mag_list = [dir('xxx*/mt.sft/mean.rbold.nii.gz'); ...
			dir('xxx*/mt.sft/mean.bold_mdant.nii.gz')];

for index=1:length(mag_list)
	fprintf('++ Begin analyzing %s in %s\n',mag_list(index).name,mag_list(index).folder);

	cd(mag_list(index).folder);
	[mag,magh] = rest_ReadNiftiImage(mag_list(index).name);

	if contains(mag_list(index).name,'bold_mdant')
		swinm='mean.rbold.nii.gz';
	else
		swinm=mag_list(index).name;
	end

	qsm_folder = replace(mag_list(index).folder,'mt.sft','qsm_postproc.sft');

	if isfolder(qsm_folder)
		cd(qsm_folder);
		[qsm,qsmh] = rest_ReadNiftiImage(swinm);

		uwh_folder = replace(qsm_folder,'qsm_postproc.sft','unwrapped_phase.sft');
		cd(uwh_folder);
		[uw,uwh] = rest_ReadNiftiImage(swinm);

		tph_folder = replace(qsm_folder,'qsm_postproc.sft','tissue_phase.sft');
		cd(tph_folder);
		[tph,tphh] = rest_ReadNiftiImage(swinm);

		mask_folder = replace(qsm_folder,'qsm_postproc.sft','mt.sft');
		cd(mask_folder);
		[mask,~] = rest_ReadNiftiImage('brain_mask_aftmc.nii.gz');
		
		swi_folder = replace(qsm_folder,'qsm_postproc.sft','swi_positive.sft');
		cd(swi_folder);

		if isfile(['veinness.' extractBefore(mag_list(index).name,'.gz')])==0 && isfile(['veinness.' mag_list(index).name])==0
			[swi,swih] = rest_ReadNiftiImage(swinm);

			[nx,ny,nz] = size(mag);
			mag1 = zeros(nx+40,ny+40,nz+40);
			uw1 = zeros(nx+40,ny+40,nz+40);
			tph1 = zeros(nx+40,ny+40,nz+40);
			qsm1 = zeros(nx+40,ny+40,nz+40);
			mask1= zeros(nx+40,ny+40,nz+40);

			mag1(21:end-20,21:end-20,21:end-20) = mag;
			uw1(21:end-20,21:end-20,21:end-20) = uw;
			tph1(21:end-20,21:end-20,21:end-20) = tph;
			qsm1(21:end-20,21:end-20,21:end-20) = qsm;
			mask1(21:end-20,21:end-20,21:end-20) = mask;


			[shearSys]=make_Shear_sys(mag1,4);
			r2=0;
			[vein_seg1,vesselness1]=vessel_seg(mag1,uw1,tph1,qsm1,mask1,1,4, 0,shearSys,0,r2,[0.795 0.795 0.840],[6 6 6],1, mask1,20);

			vein_seg = vein_seg1(21:end-20,21:end-20,21:end-20);
			vesselness = vesselness1(21:end-20,21:end-20,21:end-20);

			rest_Write4DNIfTI(vein_seg,swih,['veinseg.' extractBefore(mag_list(index).name,'.gz')]);
			rest_Write4DNIfTI(vesselness,swih,['veinness.' extractBefore(mag_list(index).name,'.gz')]);
		end
	end
end

%%%%% artery segmenation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
dataDir='/media/yuhui/xxx';
cd(dataDir);


mag_list = [dir('xxx/vaper.sft/mean.sub_d_bold.beta50.nii.gz'); ...
			dir('xxx/mt.sft/mean.vaper.beta50.nii.gz')];


for index=1:length(mag_list)
	fprintf('++ Begin analyzing %s in %s\n',mag_list(index).name,mag_list(index).folder);

	cd(mag_list(index).folder);
	if isfile(['artyseg_d5.' extractBefore(mag_list(index).name,'.gz')])==0 && isfile(['artyseg_d5.' mag_list(index).name])==0
	
		qsm_folder = replace(mag_list(index).folder,'vaper.sft','qsm_postproc_d5.sft');
		qsm_folder = replace(qsm_folder,'mt.sft','qsm_postproc_d5.sft');

		if isfolder(qsm_folder)
			cd(qsm_folder);
			[qsm,qsmh] = rest_ReadNiftiImage('mean.rbold.nii.gz');

			uwh_folder = replace(qsm_folder,'qsm_postproc_d5.sft','unwrapped_phase_d5.sft');
			cd(uwh_folder);
			[uw,uwh] = rest_ReadNiftiImage('mean.rbold.nii.gz');

			tph_folder = replace(qsm_folder,'qsm_postproc_d5.sft','tissue_phase_d5.sft');
			cd(tph_folder);
			[tph,tphh] = rest_ReadNiftiImage('mean.rbold.nii.gz');

			mask_folder = replace(qsm_folder,'qsm_postproc_d5.sft','mt.sft');
			cd(mask_folder);
			[mask,~] = rest_ReadNiftiImage('brain_mask_aftmc_d5.nii.gz');

			swi_folder = replace(qsm_folder,'swi_positive.sft','');
			cd(swi_folder);
			[swi,swih] = rest_ReadNiftiImage('mean.rbold.nii.gz');

			cd(mag_list(index).folder);
			[mag,magh] = rest_ReadNiftiImage(mag_list(index).name);

			[nx,ny,nz] = size(mag);
			mag1 = zeros(nx+40,ny+40,nz+40);
			uw1 = zeros(nx+40,ny+40,nz+40);
			tph1 = zeros(nx+40,ny+40,nz+40);
			qsm1 = zeros(nx+40,ny+40,nz+40);
			mask1= zeros(nx+40,ny+40,nz+40);

			mag1(21:end-20,21:end-20,21:end-20) = mag;
			uw1(21:end-20,21:end-20,21:end-20) = uw;
			tph1(21:end-20,21:end-20,21:end-20) = tph;
			qsm1(21:end-20,21:end-20,21:end-20) = qsm;
			mask1(21:end-20,21:end-20,21:end-20) = mask;


			[shearSys]=make_Shear_sys(mag1,4);
			r2=0;
			[vein_seg1,vesselness1]=vessel_seg(mag1,uw1,tph1,qsm1,mask1,1,4, 1,shearSys,0,r2,[0.795 0.795 0.840],[12 12 12],0, mask1,20);

			vein_seg = vein_seg1(21:end-20,21:end-20,21:end-20);
			vesselness = vesselness1(21:end-20,21:end-20,21:end-20);

			rest_Write4DNIfTI(vein_seg,magh,['artyseg_d5.' extractBefore(mag_list(index).name,'.gz')]);
			rest_Write4DNIfTI(vesselness,magh,['artyness_d5.' extractBefore(mag_list(index).name,'.gz')]);
		end
	end
end


