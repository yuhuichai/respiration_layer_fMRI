%% swi_qsm
%% replace xxx with your own data and subject directory
clear;

subjList={'xxx' 'xxx'};
funcList={'bold1' 'bold2' 'dant1' 'dant2'};

for subj=1:length(subjList)
	for func=1:length(funcList)

		fprintf('++ Begin analyzing %s in %s\n',char(funcList(func)),char(subjList(subj)));

		phanii_file = ['/media/yuhui/xxx/' char(subjList(subj)) '/pha.sft/pharad_' char(funcList(func)) '.nii'];
		magnii_file = ['/media/yuhui/xxx/' char(subjList(subj)) '/mt.sft/' char(funcList(func)) '.nii'];
		mask_file = ['/media/yuhui/xxx/' char(subjList(subj)) '/mt.sft/brain_mask_aftmc_d5.nii'];
		if isfile(mask_file)==0
			system(['gzip -v -d ' mask_file '.gz']);
		end
		if isfile(phanii_file)==0
			system(['gzip -v -d ' phanii_file '.gz']);
		end
		if isfile(magnii_file)==0
			system(['gzip -v -d ' magnii_file '.gz']);
		end

		if isfile(phanii_file)

			%% swi
			input(1).name = phanii_file;
			input(2).name = magnii_file;
			output_basename = ['/media/yuhui/xxx/' char(subjList(subj)) '/swi_postproc.sft/' char(funcList(func))];
			algorParam.swi.m = 4 ;
			algorParam.swi.threshold = pi ;
			algorParam.swi.filterSize = 12 ;
			algorParam.swi.method = 'default' ;
			algorParam.swi.isPositive = 1 ;
			algorParam.swi.isNegative = 1 ;
			algorParam.swi.ismIP = 1 ;
			algorParam.swi.slice_mIP = 4 ;

			if isfile(['/media/yuhui/xxx/' char(subjList(subj)) '/swi_postproc.sft/' char(funcList(func)) '_swi-positive.nii'])==0
				SWIIOWrapper(input,output_basename,algorParam);
			end

			if isfile(['/media/yuhui/xxx/' char(subjList(subj)) '/tissue_phase_d5.sft/' char(funcList(func)) '.nii'])==0
				% qsm
				magnii = load_nii(magnii_file);
				mag = magnii.img;

				phanii = load_nii(phanii_file);
				pha = phanii.img;

				masknii = load_nii(mask_file);
				BrainMask = masknii.img;

				voxelsize=[ 0.795  0.795  0.840 ];
				padsize=[12 12 12];
				smvsize = 12;

				H=[0 0 1];
				B0=7;
				TE=20;
					
				[xn,yn,zn,vn] = size(pha);
				qsm = zeros(xn,yn,zn,vn);
				tpha = zeros(xn,yn,zn,vn);
				uwpha = zeros(xn,yn,zn,vn);

				for vol=1:vn

					fprintf('++ Start with volume %d ... \n',vol);

					rawphase = pha(:,:,:,vol);

					%% [1] The Laplacian-based phase unwrapping
					% Inputs:
					% rawphase: 3D raw phase
					% voxelsize: Spatial resolution
					% padsize: size for padarray to increase numerical accuracy

					[Unwrapped_Phase, Laplacian]=MRPhaseUnwrap(rawphase,'voxelsize',voxelsize,'padsize',padsize);
					% % If optional inputs are not assigned, the following default values will be used.
					% voxelsize=[1 1 1];
					% padsize=[12 12 12];

					%% [3] V-SHARP: background phase removal for 3D GRE scan
					% Inputs:
					% Unwrapped_Phase: 3D unwrapped phase using Laplacian unwrapping
					% BrainMask: brain mask 
					% voxelsize: spatial resoluiton 
					% smvsize: filtering size, default value = 12
					% Output:
					% TissuePhase: 3D Filtered TissuePhase
					% NewMask: eroded mask

					[TissuePhase,NewMask]=V_SHARP(Unwrapped_Phase,BrainMask,'voxelsize',voxelsize,'smvsize',smvsize);
					% % If optional inputs are not assigned, the following default values will be used.
					% smvsize=12;
					% voxelsize=[1 1 1];

					%% [5] fast STAR-QSM (~30 sec): Quantative Susceptibility Mapping
					% Inputs:
					% TissuePhase: tissue phase
					% SpatialRes: Spatial resolution
					% padsize: size for padarray to increase numerical accuracy
					% H: the field direction, e.g. H=[0 0 1];
					Susceptibility = QSM_star(TissuePhase,NewMask,'TE',TE,'B0',B0,'H',H,'padsize',padsize,'voxelsize',voxelsize);
					% % If optional inputs are not assigned, the following default values will be used.
					% H=[0 0 1];
					% voxelsize=[1 1 1];
					% padsize=[12 12 12];
					% B0=3;
					% TE=40; 

					qsm(:,:,:,vol) = Susceptibility;
					tpha(:,:,:,vol) = TissuePhase;
					uwpha(:,:,:,vol) = Unwrapped_Phase;
				end

				if func==1
					mkdir(['/media/yuhui/xxx/' char(subjList(subj)) '/qsm_postproc_d5.sft']);
					mkdir(['/media/yuhui/xxx/' char(subjList(subj)) '/tissue_phase_d5.sft']);
					mkdir(['/media/yuhui/xxx/' char(subjList(subj)) '/unwrapped_phase_d5.sft']);
				end

				cd(['/media/yuhui/xxx/' char(subjList(subj)) '/qsm_postproc_d5.sft']);
				qsmfile = phanii;
				qsmfile.filetype = 64;
				qsmfile.img = flipud(qsm);
				qsmfile.fileprefix = [char(funcList(func))];
				save_nii(qsmfile,[char(funcList(func)) '.nii']);

				
				cd(['/media/yuhui/xxx/' char(subjList(subj)) '/tissue_phase_d5.sft']);
				qsmfile = phanii;
				qsmfile.filetype = 64;
				qsmfile.img = flipud(tpha);
				qsmfile.fileprefix = [char(funcList(func))];
				save_nii(qsmfile,[char(funcList(func)) '.nii']);

				cd(['/media/yuhui/xxx/' char(subjList(subj)) '/unwrapped_phase_d5.sft']);
				qsmfile = phanii;
				qsmfile.filetype = 64;
				qsmfile.img = flipud(uwpha);
				qsmfile.fileprefix = [char(funcList(func))];
				save_nii(qsmfile,[char(funcList(func)) '.nii']);


			end
		end
	end
end




