# respiration_layer_fMRI
Codes used in the article "Improving laminar fMRI specificity by reducing macrovascular bias revealed by respiration effects".

(1) Script used to split the original time series into even (CTRL, can be treated as BOLD signal) and odd (DANTE-prepared images in functional runs or MT-prepared in anatomical runs) time points, and to create a mask for motion correction: https://github.com/yuhuichai/respiration_layer_fMRI/blob/main/split_ctrl_dant.sh This script reads the nifti images of all runs and relies on AFNI program. 

(2) Script used for motion correction: https://github.com/yuhuichai/respiration_layer_fMRI/blob/main/mc_run.m This script reads all functional and anatomical runs, replacing the input in mc_job.m with the corresponding nifti names. It depends on SPM12 and REST (http://restfmri.net/forum/) package. 

(3) Script used to censor time points whenever the Euclidean norm of the motion derivatives exceeded 0.4 mm or when at least 10% of image voxels were seen as outliers from the trend: https://github.com/yuhuichai/respiration_layer_fMRI/blob/main/motion_censor.sh It reads the motion parameters estimated by SPM12 and relys on AFNI programs. 

(4) Script used to generate VAPER time series (sub_d_bold1/2/3/....nii.gz) and the antomical image (mean.sub_d_dant.beta100.denoised.nii.gz): https://github.com/yuhuichai/respiration_layer_fMRI/blob/main/vaper.sh It reads all motion-corrected runs and computes subtraction and ratio operations between CTRL and DANTE or MT time series. For anatomical runs, it computes the mean antomical image and denoise it using ANTs (https://github.com/ANTsX/ANTs) program DenoiseImage. 

(5) Script used to do brain segmentation and cortical surface reconstruction: https://github.com/yuhuichai/respiration_layer_fMRI/blob/main/reconall_mtepi.sh It relies on AFNI and FreeSurfer programs. 

(6) Script used to grow cortical layers: https://github.com/yuhuichai/respiration_layer_fMRI/blob/main/layer_MT.sh With the cortical surface automatically generated by FreeSurfer, we calculated cortical depths based on the equi-volume approach (Waehnert et al., 2014) using the LAYNII software suite (Huber et al., 2020) and divided the cortex into 18 equi-volume layers. Since a voxel in the acquired spatial resolution (0.8mm) can lie across several cortical depths, MT-weighted anatomical image was upsampled by a factor of 4 for the cortical layer reconstruction. 

(7) Script used to reconstruct SWI and QSM images based on the magnitude and phase image series of MT-3D-EPI volumes: https://github.com/yuhuichai/respiration_layer_fMRI/blob/main/swi_qsm.m It uses programs from SEPIA (https://doi.org/10.1016/j.neuroimage.2020.117611) and STI Suite (https://people.eecs.berkeley.edu/~chunlei.liu/software.html) toolboxs

(8) Script used to segment veins and arteries from the MT-EPI averaged magnitude, QSM and SWI volumes: https://github.com/yuhuichai/respiration_layer_fMRI/blob/main/swi_vessel.m It uses programs developed by (Straub et al., 2022) (https://doi.org/10.1016/j.neuroimage.2022.118931). Script of vessel_seg_chai.m is modified from vessel_seg_loc.m in (Straub et al., 2022) and it is used to generate vessel segmentation with a low threshold strategy (vessel_seg2).

(9) Script used to compute the respiration effect during a normal fMRI run (participants breathed naturally): https://github.com/yuhuichai/respiration_layer_fMRI/blob/main/censor_deep_breath.sh It uses multiple common programs from AFNI. 

(10) Script used to do spatial correlation analysis of respiration-related fMRI variations and here we showed an example of correction between deep and normal breath: https://github.com/yuhuichai/respiration_layer_fMRI/blob/main/corr_deepnormbreath.m
