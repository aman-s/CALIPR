function SockeyeCluster_CALIPR_Philips_Brain_FullProcess(fulldat_names)
	%%%%%%% 
	% Input is a cell array containing strings with the full path to the .raw data file to be processed
	% Example:
	% 	fulldat_names = {...
	% 	'/home/bizon/Research/MWI_Development/CALIPR/CALIPR_Publication_Share/CALIPR_Code/test_calipr_repro_c05_calipr_brain_1.raw',...
	% 	'/home/bizon/Research/MWI_Development/Data/Wilman_D_YEG/Wilman_D_YEG_calipr_brain_1.raw',...
	% 	};
	% 	SockeyeCluster_CALIPR_Philips_Brain_FullProcess( fulldat_names );
	% Or input directly as:
	% 	SockeyeCluster_CALIPR_Philips_Brain_FullProcess( 	...
	% 	'/home/bizon/Research/MWI_Development/CALIPR/CALIPR_Publication_Share/CALIPR_Code/test_calipr_repro_c05_calipr_brain_1.raw',...
	% 	'/home/bizon/Research/MWI_Development/Data/Wilman_D_YEG/Wilman_D_YEG_calipr_brain_1.raw',... 
	% 	} );


	% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	% ------------------------------------------------------------------------------ START INPUTS, PARAMETERS, OPTIONS 
	% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	format compact

	% % --------------------------------------- Software modules requires to be loaded:

	% %For Matlab:
	% module load Software_Collection/ARC_2023
	% module load gcc openmpi openblas/0.3.15
	% module load matlab/R2023a

	% %For BART:
	% module load bart/0.7.00-cuda11-3

	% %For dcm2niix, FSL, ANTs (all in containers)
	% module load apptainer

	% % --------------------------------------- Accessing tools inside containers:
	% % Note: very long spin up+down time for exec commands requiring --fakeroot option (~80 seconds total)

	% dcm2niix
	% --->
	% apptainer exec /project/st-kolind-1/ken-arc/cbh-b-tools.sif dcm2niix

	% fslmaths
	% --->
	% apptainer exec /project/st-kolind-1/fsl-6.simg fslmaths

	% antsBrainExtraction.sh
	% --->
	% ???

	% sct_deepseg_sc
	% --->
	% apptainer exec --fakeroot /project/st-kolind-1/ken-arc/cbh-b-tools.sif /opt/spinalcordtoolbox/bin/sct_deepseg_sc


	% ------------------------------------------------------------------------------  SETUP COMPUTE CLUSTER FOR MATLAB
	% REQUIRED SETTINGS
	
	oasis_template_path = '/home/amanstaf/templates/oasis/MICCAI2012-Multi-Atlas-Challenge-Data'; %/project/st-kolind-1/Templates/OASIS

	% Connect:
	c = parcluster;
	% Specify our allocation (!Required):
	c.AdditionalProperties.AllocationCode = 'st-kolind-1';

	% OPTIONAL SETTINGS
	% Specify number of GPUs and memory
	c.AdditionalProperties.GpuMem = '32gb';
	c.AdditionalProperties.GpusPerNode = 1;
	% Specify job placement flexibility (in this case each worker can run on any node, not requiring a specific one):
	c.AdditionalProperties.JobPlacement = 'free';
	c.AdditionalProperties.RequireExclusiveNode = false;
	% Specify memory to use for MATLAB jobs, per core
	% Max 24 CPU cores and 192GB RAM (8GB/core)
	c.AdditionalProperties.MemUsage = '6000mb'
	% Request 24 procs per node (instead of default 8)
	c.AdditionalProperties.ProcsPerNode = 24;
	% Specify the walltime (e.g. 1 Hour)
	% >> c.AdditionalProperties.WallTime = '01:00:00';

	% NOTE: To save changes after modifying AdditionalProperties for the above changes to persist between MATLAB sessions.
	% >> c.saveProfile


	% ------------------------------------------------------------------------------  SETUP DEPENDENCIES


	% --------------------------------------- BART TOOLBOX SETTINGS
	%addpath(fullfile('/arc/software/spack-2023/opt/spack/linux-centos7-skylake_avx512/gcc-9.4.0/bart-0.7.00-en3t3cfmmfbdv2xxzo4d2xnlwpvsw6rs', 'matlab'));
	%setenv('TOOLBOX_PATH', '/arc/software/spack-2023/opt/spack/linux-centos7-skylake_avx512/gcc-9.4.0/bart-0.7.00-en3t3cfmmfbdv2xxzo4d2xnlwpvsw6rs');
	addpath("/project/st-kolind-1/aman_bart/bart-matlab-integration-sockeye-cluster/")
	disp('BART Verson info:'); bart('version -V')
	setenv('OMP_NUM_THREADS','24'); % set based on SLURM job settings
	% export DEBUG_LEVEL=5 % for BART debugging

	% --------------------------------------- GYROTOOLS/MRECON SETTINGS AND LICENSE
	mreconpath = '/arc/software/MRecon-5.1.0';
	addpath(genpath(fullfile(mreconpath)))
	disp('Gyrotools MRecon info:'); MRecon.LicenseInfo
	% MRecon.CheckForUpdates
	% doc MRecon % For MRecon documentation

	% Note on MRecon licensing:
	% Ken ran in matlab with admin privileges, to license all users:
	% MRecon.Activate('<ACTIVATION TOKEN>') 
	% Creates license key in /etc/gyrotools/
	% tested and working for all accounts, but may need to move license key to alternative permanent location

	% --------------------------------------- DK MWI ANALYSIS
	% BEFORE RUNNING THE MWI ANALYSIS CODE: make sure that mex file of L-BFGS-B is compiled:
	% Go to ...\Subcodes\L-BFGS-B-C-master\Matlab and run “compile_mex”
	addpath(genpath('/project/st-kolind-1/DK_MWI_Analysis/'))

	% --------------------------------------- ALSO REQUIRED IN PATH
	% - ANTs
	% - FSL
	% - SCT for cord processing
	% - dcm2niix for converting PAR/REC to nifti
	% - OASIS template used for brain extraction

	% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	% ------------------------------------------------------------------------------ INPUTS 

	% fulldat_names = {...
	% '/home/bizon/Research/MWI_Development/CALIPR/CALIPR_Publication_Share/CALIPR_Code/test_calipr_repro_c05_calipr_brain_1.raw',...
	% '/home/bizon/Research/MWI_Development/Data/Wilman_D_YEG/Wilman_D_YEG_calipr_brain_1.raw',...
	% };


	% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	% ------------------------------------------------------------------------------ HARDWARE CORRECTION

	% ----------
	run_recon = 1;
	ZeroFill = 'Yes'; % 'Yes', 'No'
	save_noiseDAT = 1;
	use_bartPW = 1;
	% ----------
	save_output = 0;
	% ----------
	run_CC = 1;
	% ----------

	% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	% ------------------------------------------------------------------------------ BRAIN RECON
	% clear MR Eall Eksp Etmp
	% datids = [ 1 ]; espids = [ 1 ]; % 3 1 8 6  datids = [ 1 2 3 4 5 6 7 8 9 ]; espids = [ 1 ];

	load_data=0;
	ZeroFill = 'Yes'; % 'Yes', 'No'
	use_bartPW = 1;

	run_espirit=0;

	run_CSrecon = 0;
	save_cs_nii = 1;

	% run_SCrecon = 0;

	run_CALIPRrecon = 0;
	save_calipr_nii = 1;

	analyseMWI = 0;
	save_mwi_nii = 1;


	% -------------------------------------------------------------------- BRAIN

	% ---------------------------------- W Zerofill 240x200x122
	% recon_slices = 119:122; view_slice = 3; batchsize = 4;% 4 slices
	% recon_slices = 120:121; view_slice = 2;  % 2 slices
	% recon_slices = 119:122; view_slice = 3;  % 4 slices
	% recon_slices = 109:132; view_slice = 13;  % 24 slices
	% recon_slices = 99:142; view_slice = 23;  % 44 slices
	% recon_slices = 1:240; view_slice = 124;  batchsize = 6;% ALL 240 slices, batchsize 6 -> 40 batches
	recon_slices = 1:240; view_slice = 124;  batchsize = 24;% ALL 240 slices, batchsize 24 -> 10 batches  USED FOR ALL CALIPR_REPRO RECONS
	% recon_slices = 1:160; view_slice = 80;  batchsize = 20;% ALL 240 slices, batchsize 24 -> 10 batches  USED FOR CALIPR_GEMATCH RECONS
	% recon_slices = 1:240; view_slice = 124;  batchsize = 30;% ALL 240 slices, batchsize 30 -> 8 batches
	% Total Memory During Temporal MWI Analysis: 110Gb/81% (for total entire image, processing within entire threshmask, 24 workers, clear Eall etc before, nT2 40)
	% Total Memory During Spatial MWI Analysis: ?Gb/?% (for total entire image, processing within entire threshmask, 24 workers, clear Eall etc before, nT2 40, DSW 10x10x10)




	% ---- To debug:
	debugi = 0;
	% debugi = debugi +1; disp(sprintf('Debug %u',debugi))
	% ----

	espirit_cmds = {...
	'ecalib -d 0 -r 16 -S -a ',... % -d 0 for quiet version
	};


	% --------------------------------- PREP GENERAL PARAMS
	print_output = 0;
	echoes2show = [ ]; % [ 1 4 10 28 ]
	show_coeff = 1;
	imwindow = [ 0 0.020 ]; % [ 0 0.018 ] [0 900]
	ax_shift = 0; % shift all axial slices up or down
	view_slices_ax = round(1.7*[ 12, 15, 18, 24, 32 ]); %  [ 3, 4, 6, 8, 10 ] [ 4, 6, 8, 10, 12 ] [ 8, 12, 16, 20, 24 ] [ 6, 8, 10, 12, 14 ] [ 28, 44, 60, 76, 92 ] % which 5 axial slices to view
	calc_error = 0;
	% echoes2save = [ 1 ]; %  [ 1 7 14 21 28 56 ];
	% --------------------------------- PREP RECON
	% calipr_cmd_pics = 'pics -R W:70:0:0.1 -i 100 -e -d 5 -S -g '; % FOR OVERWRITING
	%         cs_cmd_pics = 'pics -R W:6:0:0.005 -i 250 -e -d 4 -S -H -g '; % was 250 FOR OVERWRITING
	%         cs_cmd_pics = 'pics -R W:6:0:0.005 -i 100 -e -d 4 -S -H '; %
	%         cs_cmd_pics = 'pics -R W:6:0:0.0001 -i 250 -e -d 4 -S -H -g'; % 
	% cs_cmd_pics = 'pics -R W:6:0:0.002 -i 250 -e -d 4 -S -H -g '; % 
	% cs_cmd_pics = 'pics -R W:6:0:0.0004 -i 250 -e -d 4 -S -H -g'; % 
	cs_cmd_pics = 'pics -R W:6:0:0.004 -i 250 -e -d 4 -S -H -g'; % USED FOR CALIBRATION FOR ALL CALIPR_REPRO BRAIN CS RECONS
	% cs_cmd_pics = 'pics -R W:6:0:0.01 -i 250 -e -d 4 -S -H -g'; % USED FOR CALIPR_GEMATCH RECONS
	% cs_cmd_pics = 'pics -R W:6:0:0.4 -i 1 -e -d 4 -S '; % 
	% cs_cmd_pics = 'pics -R W:6:0:0.0005 -i 250 -e -d 4 -S -H -g'; % 

	%         cs_cmd_pics = 'pics -R W:6:0:0.05 -i 250 -e -d 4 -S -H -g'; % 
	%         cs_cmd_pics = 'pics -l1 0.005 -i 250 -e -d 4 -S -H -g '; % 
	%         cs_cmd_pics = 'pics -i 1 -e -d 5 -S -H -g '; % 
	%         cs_cmd_pics = 'pics -R W:6:0:0.01 -i 100 -w 1269 -e -d 5 -S -H -g '; % TESTING
	%         cs_cmd_pics = 'pics -R T:6:0:0.1 -i 100 -w 1269 -e -d 5 -S -H -g '; % TESTING
	%         cs_cmd_pics = 'pics -R W:6:0:0.005 -i 100 -e -d 5 -S -H -g '; % TESTING
	% calipr_cmd_pics = 'pics -R W:6:0:0.1 -i 30 -e -d 5 -S -g '; % FOR OVERWRITING
	ALL_cmd_pics = {...
	% -------------------- CALIPR Oct22 loop
	% 'pics -R W:6:0:0.0020 -i 250 -e -d 4 -S -H -g ',...
	'pics -R W:6:0:0.0040 -i 250 -e -d 4 -S -H -g',... % USED FOR ALL CALIPR_REPRO BRAIN CALIPR RECONS
	% 'pics -R W:6:0:0.0100 -i 250 -e -d 4 -S -H -g',... % USED FOR ALL CALIPR_GEMATCH RECONS
	% 'pics -R W:6:0:0.0200 -i 250 -e -d 4 -S -H -g',... % USED FOR ALL CALIPR_GEMATCH RECONS
	% 'pics -R W:6:0:0.0008 -i 250 -e -d 4 -S -H -g',...
	% 'pics -R W:6:0:0.0060 -i 250 -e -d 4 -S -H -g',...
	% 'pics -R W:6:0:0.0004 -i 250 -e -d 4 -S -H -g',...
	% 'pics -R W:6:0:0.0080 -i 250 -e -d 4 -S -H -g',...
	% 'pics -R W:6:0:0.0100 -i 250 -e -d 4 -S -H -g',...
	% 'pics -R W:6:0:0.0200 -i 250 -e -d 4 -S -H -g',...
	% 'pics -R W:6:0:0.0001 -i 250 -e -d 4 -S -H -g',...
	% --------------------
	% 'pics -R W:6:0:0.1 -i 100 -e -d 5 -S -g -H ',... % DEFINITELY BETTER! Identical with/without -e
	% 'pics -R W:6:0:0.1 -i 30 -e -d 5 -S -g ',...
	% 'pics -R W:6:0:0.1 -i 100 -e -d 5 -S -g ',...
	% TRIED:
	% 'pics -R L:6:6:0.1 -i 30 -e -d 5 -S ',... % Note: would not run on GPU
	% 'pics -R W:6:0:0.1 -R L:6:6:0.1 -i 30 -e -d 5 -S -b 12',...
	% 'pics -R L:6:6:0.1 -i 30 -e -d 5 -S -b 12 -N',... % takes absolutely ages, had to kill
	% 'pics -R M:6:6:0.1 -i 30 -e -d 5 -S ',... % ERROR would not run
	% 'pics -R W:70:0:0.1 -i 100 -e -d 5 -S -g ',... % Makes coefficient images super noisy for SpatTemp, but weirdly WAY BETTER with SpatOnly sampling
	% 'pics -R W:6:64:0.1 -i 100 -e -d 5 -S -g ',...
	% 'pics -R W:6:8192:0.1 -i 100 -e -d 5 -S -g ',...
	% 'pics -R W:6:6:0.1 -i 100 -e -d 5 -S -g ',... % ERROR will not run
	% 'pics -R W:6:0:0.1 -i 100 -e -d 5 -S -g -a',... % ERROR will not run
	};
	% --------------------------------- PREP SC BASIS
	% TEMPORARY
	% sc_basis_file = '/home/bizon/Research/MWI_Development/mcQT2_Shuffling/CALIPER_Dev_FB1_MWI/Basis/datU_cpx_x80.mat'; %80/88, comb/cpx/imag/mag/real
	sc_basis_file = [ '/home/bizon/Research/MWI_Development/CALIPR/MWI_RetroUS/se56_fs/se56_fs_SCbasis_simU.mat' ];
	% --------------------------------- PREP CALIPR CALIBRATION DATA
	calibrate_dat = 'UScs'; % options:  'tinyUScs' 'tinyUSfft_rss' 'tinyUSfft_fmac' 'UScs'
	calibrate_datamatch = 0; % only for tiny calibrate versions (nonsensical to use 0 for tinyUSfft_rss if differing cal regions)
	tinyespirit_cmd = 'ecalib -d 4 -r 16 -S -a ';
	% tinyespirit_cmd = 'ecalib -d 4 -r 0:16:10 -k 0:6:4 -S -m 1 -a ';
	% tinyespirit_cmd = 'ecalib -d 4 -m 2 -S -a '; % 'ecalib -d 4 -m 2 -S '  only for calibrate_dat = 'tinyUScs' or 'tinyUSfft_fmac'
	tinycmd_pics = 'pics -R W:6:0:0.0004 -i 250 -e -d 4 -S -H -g -w 404 ';
	% tinycmd_pics = 'pics -R W:6:0:0.0004 -i 250 -e -d 4 -S -H -g'; % only for calibrate_dat = 'tinyUScs'
	% --------------------------------- PREP CALIPR BASIS
	Kvals = [ 12 ]; % 4 5  % subspace size
	basis_ids = [ 1 ]; % 1 2 3 4 5 6 7 8 basis IDs to loop over, 1:8 total % NOTE: makes no difference for RSS, since absolute value already
	show_threshold = 1;
	show_subspace = 1;
	show_sig_evol = 1;
	thresh = 1/15 % BRAIN threshold at fraction of max intensity voxel
	% thresh = 1/25 % CORD threshold at fraction of max intensity voxel
	maxN = 25000000; % maximum number of unique signal evolutions to use for SVD
	maxN_disp = 50; % maximum number of unique signal evolutions to display
	% --------------------------------- FINAL OPTIONS NOT OFTEN CHANGED
	espirit_use_sets = 1;
	espirit_combine_sets = 0;
	apply_nii_mask = 0;
	apply_mat_mask = 0;
	mask_borderpad = 2; % depth to fill then cut in mask creation
	% ----- NEW DK MWI params -----
	Resol = [1.000, 1.000, 1.000]; % [0.625, 0.625, 2.5] [1.0, 1.0, 5.0] cord with ZF [1.0, 1.0, 2.0] [2.0, 2.0, 2.0]
	TE=0.006000; % 0.008720, 0.007168, 0.007000
	TR=1.252; % 1.013, 1.667, 1.120, 1.306
	T2Range=[0.00400,10.00000]; % 0.00600, 0.00500
	spwin = [0.00000; 0.04000]; % 0.00100; 0.04000
	mpwin=[0.04000,0.20000]; % [0.04000,0.20000]
	nT2=40; % 120
	nCores=12; % 30 -> ?% of CPU
	nThreads=24; % 30 -> ?% of CPU
	FAE_in_fin_step = [0, 50, 1];
	RegType='Both'; % 'Both', 'Temporal-Only', 'Spatial-Only'
	nSpatIter=2;
	Size_DSW = [10, 10, 10]; % 24x24x12 with 20 workers.    ~25:00 for 6iter/60v30p/10x10x10/mask10slice
	Overlap_DSW = 4;
	SpatRegOpt='Optimize muS'; % 'Optimize muS', 'Use set value'
	SpatRegSetAlpha=3000;
	% Initial values for determining muS (spatial regularizatio const for flip angle inhomogeneities)
	SpatReg_Init_FA=double(logspace(1,5,8)); % 10 to 100,000 logarithmically spaced
	% Initial values for determining muS (spatial regularizatio const for T2 dist. map) 
	SpatReg_Init_T2=[10, 100, 200, 350, logspace(log10(5e2), log10(3.5e3), 10), 5e3, 1e4];    % [100, 200, 500, 1000, 2000, 4000, 7500, 10000];
	% ----- Not inputs for analysis, used afterwards
	scale_MWF = [ 0 0.20 ];
	scale_ALPHA = [ 120 180 ];
	% --------------------------------- 
	% % ----- Set extra PM MWI params -----
	% Reg='yes';
	% RefCon=180;
	% nAngles=8; % 64
	% MinRefAngle=90;
	% Chi2Factor=1.02;
	% nii_brainmask_fullfilename = '/home/bizon/Research/MWI_Development//Data/mcQT2Shuffle_C07/SE56_Rsz_E1_mask_brain.nii';
	% --------------------------------- 


	% --------------------------------- OUTPUT FILES
	R_ids = [ 1 ]; % 3:4  or 1:size(usmasks,5)  to do all
	mas_ids = [ 1 ]; % [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 ] or  1:size(USmasks,4)  to do all


	batchsize
	nbatches = length(recon_slices)/batchsize

	% ----------------------------------------------------------------------------------------------------------


	% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	% ------------------------------------------------------------------------------ END INPUTS, PARAMETERS, OPTIONS
	% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------



	% ---------------------------------- VERSION 2.0 GROUPWISE  IMPORT, HARDWARE CORRECT with ENTIRE dataset (all echoes) at once
	for dat_id = 1:length(fulldat_names)
	    tic
	    
	    [ subject_dir dat_name dat_extension ] = fileparts(fulldat_names{dat_id});
	    subject_dir = [ subject_dir '/' ];
	    out_dir = [ subject_dir dat_name '/' ];
	    
	    mkdir(out_dir)
	    cd(out_dir);
	    
	    out_dir_appended = [ out_dir 'ProUS/']
	    mkdir(out_dir_appended)
	    out_textfile = [ out_dir 'Output_CALIPR_ProUS_Cord_ECCENTRIC_DEV10.txt']; 
	    
	    % ----- For saving output
	    disp('Starting subject:')
	    pwd
	    % ---------------------------------- ASSEMBLE NAMES
	    % Assemble full names
	    dat_labfile = [ subject_dir dat_name '.lab' ];
	    dat_rawfile = [ subject_dir dat_name '.raw' ];
	    dat_sinfile = [ subject_dir dat_name '.sin' ];
	    % ----- For saving noise data
	    noisefile = [ out_dir dat_name '_NOISEdat.mat' ];
	    % ----- For saving pre-whitening matrix/COV data
	    PWmatrixfile = [ out_dir dat_name '_PWmatrix' ];
	    noiseCOVfile = [ out_dir dat_name '_noiseCOV' ];
	    % ----- For saving output data
	    CCmatrixfile = [ out_dir dat_name '_CCmatrix' ];
	    datname_Eall = [ out_dir dat_name '_bPW' num2str(use_bartPW) '_zf' ZeroFill '_wK2IM_EALL.mat' ];
	    datname_wCC_Eall = [ out_dir dat_name '_bPW' num2str(use_bartPW) '_zf' ZeroFill '_wK2IM_CC1_EALL.mat' ];
	    datname_wCC_Eall_crop = [ out_dir dat_name '_bPW' num2str(use_bartPW) '_zf' ZeroFill '_wK2IM_CC1_EALL_crop.mat' ];
	    datname_slice_prefix = [ out_dir 'Slices/' dat_name '_bPW' num2str(use_bartPW) '_zf' ZeroFill '_wK2IM_CC0' ]; % appended with _x60.mat etc
	    datname_slice_wCC_prefix = [ out_dir 'Slices/' dat_name '_bPW' num2str(use_bartPW) '_zf' ZeroFill '_wK2IM_CC1' ]; % appended with _x60.mat etc
	    datname_echo_prefix = [ out_dir 'Echoes/' dat_name '_bPW' num2str(use_bartPW) '_zf' ZeroFill '_wK2IM_CC0' ]; % appended with _x60.mat etc
	    

	    % ----------------------------------------------------------------------------------------------------------------------------------------
	    % ---------------------------------------------------------------------------------------------------------------------------------------- HW Correct
	    % ----------------------------------------------------------------------------------------------------------------------------------------
	    if run_recon == 1

	        clear MR Eall
	        % ----- CREATE RECONSTRUCTION OBJECT, SPECIFY RECON PARAMATERS
	%         MR = MRecon(dat_labfile, dat_rawfile, dat_sinfile); % Initialize recon object ~220s
	        MR = MRecon(dat_labfile, dat_rawfile);
	        nEchoes = length(MR.Parameter.Parameter2Read.echo)

	        MR.Parameter.ExtractPDFFile([ out_dir dat_name '_SequencePDF.gve' ]); % Get sequence PDF file for viewing in GVE

	        if save_noiseDAT == 1
	            MR.Parameter.Parameter2Read.typ = [ 5 ]; % read in standard data (1) phase correction data (3) noise data (5)
	            MR.Parameter.Parameter2Read.Update;
	            MR.ReadData;

	            % ----- Most are inactivate but just to get ReconFlags and dimensions to same state as data during pre-whitening for consistency
	            MR.RandomPhaseCorrection;
	            MR.RemoveOversampling;
	            MR.PDACorrection;
	            MR.DcOffsetCorrection;
	            MR.MeasPhaseCorrection;
	            MR.SortData;
	            MR.GridData;
	            % -----
	            noisedat=MR.Data;

	            save(noisefile, 'noisedat');
	            clearvars noisedat
	        end
	        
	        mkdir([ out_dir 'Echoes/' ]);
	        
	        
	        for eid = 1:nEchoes
	            MR.Parameter.Reset;
	            
	            MR.Parameter.Parameter2Read.typ = [ 1 ]; % read in standard data (1) phase correction data (3) noise data (5)
	            MR.Parameter.Parameter2Read.echo = [ eid-1 ]; % only read in only some echoes   format: [ 0; 1; 2 ];
	            MR.Parameter.Parameter2Read.Update;
	            
	%             disp('tmp')
	            
	%             %%%%% For brain reference scan: zerofill to 1mm (160 x 128 x 80 -> 240 x 200 x 122)
	%             MR.Parameter.Encoding.XRes = 240+zeros(nEchoes,1);
	%             MR.Parameter.Encoding.XReconRes = 240+zeros(nEchoes,1);
	%             MR.Parameter.Encoding.YRes = 200+zeros(nEchoes,1);
	%             MR.Parameter.Encoding.YReconRes = 200+zeros(nEchoes,1);
	%             MR.Parameter.Encoding.ZRes = 122+zeros(nEchoes,1);
	%             MR.Parameter.Encoding.ZReconRes = 122+zeros(nEchoes,1);
	            
	%             %%%%% For cord reference scan: zerofill slice by factor of 2, to 2.5mm slices(6 -> 12)
	%             MR.Parameter.Encoding.ZRes = 12+zeros(nEchoes,1);
	%             MR.Parameter.Encoding.ZReconRes = 12+zeros(nEchoes,1);
	            
	%             %%%%% For all cord CALIPR scans: zerofill slice by factor of 2, to 2.5mm slices(16 -> 32)
	%             %%%%% OR alternatively, zerofill in recon code after loading (to reduce saved data size)
	%             MR.Parameter.Encoding.ZRes = 32+zeros(nEchoes,1);
	%             MR.Parameter.Encoding.ZReconRes = 32+zeros(nEchoes,1);

	%             disp('tmp2')
	            
	            MR.Parameter.Recon.kSpaceZeroFill = ZeroFill;
	            MR.Parameter.Recon.UseMatlabInternals = 'Yes';

	%             % ----- Print recon parameters if desired
	%             MR.Parameter.Parameter2Read
	%             MR.Parameter.Recon
	%             MR.Parameter.ReconFlags
	%             % ----- 

	            % ---------------------------------------------------- READ DATA, PERFORM CORRECTIONS
	            MR.ReadData; % disp('read') % First dimension is measurement, second is all profiles in the order they were acquired
	            MR.RandomPhaseCorrection; % disp('Random Phase Correction done')
	            MR.RemoveOversampling; % FOR SENSE UNFOLD require oversampling in the undersampled direction is not removed
	            MR.PDACorrection;
	            MR.DcOffsetCorrection;
	            MR.MeasPhaseCorrection;

	            % Now sort the data into a matrix where each dimension corresponds to an imaging parameter
	            % 12D array:    x – y – z – coils – dynamics - cardiac phases – echoes – locations – mixes – extr1 – extr2 – averages
	%             disp('presort'); size(MR.Data)
	            MR.SortData;
	%             disp('postsort'); size(MR.Data)
	            nChan = size(MR.Data,4);
	            MR.GridData; % not needed for cartesian
	%             disp('postgrid'); size(MR.Data)
	            % ----------------------------------------------------


	            % ---------------------------------------------------- NOISE PRE-WHITENING
	            if use_bartPW == 1

	                % load noise data
	                load(noisefile);

	                % rearrange data properly for BART PW
	                MR.Data = permute(MR.Data, [ 1 2 3 4 5 7 6 ]);
	    %             size(noisedat)
	    %             size(MR.Data)

	                % --------- Echo 1 separately, save matrix and COV to ensure same is used for all echoes
	%                 disp(sprintf('Echo Pre-Whitening'))
	                if eid == 1
	    %                 [ MR.Data(:,:,:,:,:,:) tmp_PWmatrix tmp_noiseCOV ] = bart('whiten ', MR.Data(:,:,:,:,:,:), noisedat); % overwrite with PW version - NO NORMALIZE
	                    [ MR.Data(:,:,:,:,:,:) tmp_PWmatrix tmp_noiseCOV ] = bart('whiten -n ', MR.Data(:,:,:,:,:,:), noisedat); % overwrite with PW version

	                    % write CFL
	                    writecfl(PWmatrixfile,tmp_PWmatrix);
	                    writecfl(noiseCOVfile,tmp_noiseCOV);

	                    % --------- Calculate and save noise covariance matrix before and after pre-whitening
	                    tmp = permute(squeeze(noisedat), [ 2 1 ]);
	                    psi = (1/(size(tmp,1)-1))*(tmp * tmp');
	                    psiscale = max(abs(psi(:)));
	                    figure(), imshow(abs(psi), [0 psiscale/10 ],'InitialMagnification','fit')
	                    title('Noise Covariance Matrix before Pre-Whitening')
	    %                 saveas(gcf,'PLOT_Error_vs_EchoNum.png')
	                    saveas(gcf,[ noiseCOVfile '_beforePW.png'])

	    %                 tmp = permute(squeeze(bart([ 'whiten -o ' PWmatrixfile ' -c ' noiseCOVfile ], noisedat, noisedat)), [ 2 1 ]); %  - NO NORMALIZE
	                    tmp = permute(squeeze(bart([ 'whiten -n -o ' PWmatrixfile ' -c ' noiseCOVfile ], noisedat, noisedat)), [ 2 1 ]);
	                    psi = (1/(size(tmp,1)-1))*(tmp * tmp');
	                    psiscale = max(abs(psi(:)));
	                    figure(), imshow(abs(psi), [0 psiscale/10 ],'InitialMagnification','fit')
	                    title('Noise Covariance Matrix After Pre-Whitening')
	                    saveas(gcf,[ noiseCOVfile '_afterPW.png'])
	                    % ---------
	                else
	    %                 MR.Data(:,:,:,:,:,:) = bart([ 'whiten -o ' PWmatrixfile ' -c ' noiseCOVfile ], MR.Data(:,:,:,:,:,:), noisedat); % overwrite with PW version - NO NORMALIZE
	                    evalc('MR.Data(:,:,:,:,:,:) = bart([ ''whiten -n -o '' PWmatrixfile '' -c '' noiseCOVfile ], MR.Data(:,:,:,:,:,:), noisedat);'); % overwrite with PW version
	                end

	                % undo data permutation, to return to MRecon format
	                MR.Data = permute(MR.Data, [ 1 2 3 4 5 7 6 ]);
	            end
	            % ----------------------------------------------------        

	            % ---------------------------------------------------- REMAINING RECON, FFT
	            MR.RingingFilter;
	%             disp('Pre ZeroFill Size:'); size(MR.Data)
	            MR.ZeroFill; 
	%             disp('Post ZeroFill Size:'); size(MR.Data)
	            MR.K2IM;
	            MR.EPIPhaseCorrection;
	%             disp('Pre second remove oversampling Size:'); size(MR.Data)
	%             MR.RemoveOversampling; % EXTRA added for testing. Does phase iFFT, crop, then FFT
	%             disp('Pre second remove oversampling Size:'); size(MR.Data)
	%             disp('Done recon'); size(MR.Data)
	            % ----------------------------------------------------

	            % ---------------------------------------------------- SAVE without CC
	%             mkdir([ out_dir 'Echoes/' ]);
	    %         tmpname = [ out_dir 'Echoes/' dat_name '_bPW_noN_wZF_wK2IM' sprintf('_E%u.mat',eid) ]; % - NO NORMALIZE
	            tmpname = [ datname_echo_prefix sprintf('_E%u.mat',eid) ];

	%             tmp_echo = bart('fftmod 6',permute(MR.Data, [ 1 2 3 4 5 7 6 ]));
	            evalc('tmp_echo = bart(''fftmod 6'',permute(MR.Data, [ 1 2 3 4 5 7 6 ]));'); % QUIET VERSION
	            % ----- SAVE without CC
	%             save(tmpname, 'tmp_echo','-v7.3', '-nocompression'); % save(tmpname, 'tmp_echo','-v7.3'); % ?FASTER new matlab file format, ensure compression not an issue
	% %             clear tmpname tmp_echo
	            % -----
	            
	            % ---------------------------------------------------- PERFORM COIL COMPRESSION AND SAVE
	            if run_CC == 1
	                if eid == 1
	                    % ----- Calculate CC matrix using Echo 1
	                    e1ksp = bart('fft -u 1', tmp_echo);

	                    % None
	                    e1_noCC = abs(squeeze(bart('rss 8', bart('fft -u -i 7', e1ksp))));
	                    % size(e1_noCC)

	                    % With CC
	                    % e1ksp_CCg_matrix = bart('cc -p 8 -M -A -G', e1ksp);
	                    e1ksp_CCg_matrix = bart('cc -p 8 -M -A -G', e1ksp); % -A using ALL k-space data
	            %         e1ksp_CCg_matrix = bart('cc -p 8 -M -r 12 -G', e1ksp); % -r 12 using k-space data calibration region of size 12
	                    save(CCmatrixfile, 'e1ksp_CCg_matrix')

	                    % ----- Apply to full dataset echo-by-echo
	                    Etmp = tmp_echo;

	                    Eksp = bart('fft -u 1', Etmp);
	                    Eksp = bart('ccapply -p 8 -G', Eksp, e1ksp_CCg_matrix);
	                    Etmp = bart('fft -u -i 1', Eksp); % Overwrite with coil compressed, fftmod version

	                    Eall = Etmp; % size(Eall)
	                else
	                    Etmp = tmp_echo;

	                    % Eksp = bart('fft -u 1', Etmp);
	                    evalc('Eksp = bart(''fft -u 1'', Etmp);');
	                    % Eksp = bart('ccapply -p 8 -G', Eksp, e1ksp_CCg_matrix);
	                    evalc('Eksp = bart(''ccapply -p 8 -G'', Eksp, e1ksp_CCg_matrix);');
	                    % Etmp = bart('fft -u -i 1', Eksp); 
	                    evalc('Etmp = bart(''fft -u -i 1'', Eksp);'); % Overwrite with coil compressed, fftmod version

	                    Eall(:,:,:,:,1,eid) = Etmp;
	            %         disp(sprintf('Echo %u',eid))
	                end
	            end
	            
	        end
	        
	        disp('Done echo-by-echo hardware correct and coil compression')
	        toc
	        
	        size(Eall)

	% %         % ---------------------------------------------------- SAVE ALL TOGETHER
	% %         save(datname_wCC_Eall, 'Eall','-v7.3', '-nocompression');
	% %         disp('Saved all echoes together')

	%         % ---------------------------------------------------- SAVE SLICES
	%         mkdir([ out_dir 'Slices/' ]);

	%         for x_id = 1:size(Eall,1)
	%             tmpname = [ datname_slice_wCC_prefix sprintf('_x%u.mat',x_id) ];
	%             tmp_slice = squeeze(Eall(x_id,:,:,:,:,:,:));

	%             save(tmpname, 'tmp_slice');
	%         %     disp(sprintf('Slice %u',x_id))
	%             clear tmpname tmp_slice
	%         end

	        disp('Done saving coil compressed slices')
	        dat_name
	        toc
	    end
	    
	    
	    % ----------------------------------------------------------------------------------------------------------------------------------------
	    % ---------------------------------------------------------------------------------------------------------------------------------------- Reconstruction and Analysis
	    % ----------------------------------------------------------------------------------------------------------------------------------------
	    
	    if load_data == 1
	        clear USdat FSsensmaps T2_dist_4D_Temporal UScalipr UScalipr_ALPHA UScalipr_MWF USfft USfft_norm pics_mas_8ch
	        % ------------------------------------------------------------------ LOAD USdat
	    %     datname_wCC_Eall
	    %     load(datname_wCC_Eall)
	    %     Eall = bart('flip 6', Eall);
	    %     Eall = Eall(recon_slices,:,:,:,:,:);
	        % ------------------------------------------------------------------ LOAD USdat ONE SLICE AT A TIME w/CC
	        [ datname_slice_wCC_prefix sprintf('_x%u.mat',recon_slices(1)) ]
	        load([ datname_slice_wCC_prefix sprintf('_x%u.mat',recon_slices(1)) ])
	        Eall = bart('flip 6', permute(tmp_slice, [ 5 1 2 3 6 4 ]));
	        count = 2;
	        for x_id = recon_slices(2:end)
	    %         [ datname_slice_wCC_prefix sprintf('_x%u.mat',x_id) ]
	            load([ datname_slice_wCC_prefix sprintf('_x%u.mat',x_id) ]);
	    %         Eall(count,:,:,:,:,:) = bart('flip 6', permute(tmp_slice, [ 5 1 2 3 6 4 ]));
	            evalc('Eall(count,:,:,:,:,:) = bart(''flip 6'', permute(tmp_slice, [ 5 1 2 3 6 4 ]));'); % quiet version
	            count = count + 1;
	        end
	        % ------------------------------------------------------------------ LOAD USdat ONE SLICE AT A TIME NO CC
	%         [ datname_slice_prefix sprintf('_x%u.mat',recon_slices(1)) ]
	%         load([ datname_slice_prefix sprintf('_x%u.mat',recon_slices(1)) ])
	%         Eall = bart('flip 6', permute(tmp_slice, [ 5 1 2 3 6 4 ]));
	%         count = 2;
	%         for x_id = recon_slices(2:end)
	% %             [ datname_slice_prefix sprintf('_x%u.mat',x_id) ]
	%             load([ datname_slice_prefix sprintf('_x%u.mat',x_id) ]);
	% %             Eall(count,:,:,:,:,:) = bart('flip 6', permute(tmp_slice, [ 5 1 2 3 6 4 ]));
	%             evalc('Eall(count,:,:,:,:,:) = bart(''flip 6'', permute(tmp_slice, [ 5 1 2 3 6 4 ]));'); % quiet version
	%             count = count + 1;
	%         end
	        % ------------------------------------------------------------------

	    %     % -------- CROP TO ACQUIRED 1MM OR 1x5MM  VERSION
	    %     Eall = bart('resize -c 1 148', Eall);
	    %     Eall = bart('resize -c 1 148 2 15', Eall);
	    %     % -------- ZEROFILL TO 1MM ISOTROPIC VERSION
	    %     Eall = bart('resize -c 1 196 2 120', Eall);
	    %     % --------
	        % -------- ZEROFILL SLICE DIRECTION FOR 0.625MM ISOTROPIC VERSION (ZF YES)
	%         Eall = bart('resize -c 2 120', Eall);
	        % --------
	%         Eall = bart('resize -c 2 75', Eall);
	        
	%         % -------- ZEROFILL SLICE DIRECTION FOR 1MM ISOTROPIC VERSION (ZF NO)
	%         Eall = bart('resize -c 2 74', Eall);
	%         % --------
	%         Eall = bart('resize -c 2 20', Eall);
	%         Eall = bart('resize -c 2 15', Eall);
	%         Eall = bart('resize -c 2 16', Eall);
	%         Eall = bart('resize -c 2 32', Eall);
	%         Eall = bart('resize -c 2 10', Eall);
	%         % -------- 0.625 isotropic 120x10 -> 240x80
	%         Eall = bart('resize -c 1 240 2 20', Eall);


	%         USdat = single(Eall); 
	%         clear Eall
	    end
	    toc
	    disp('Data loaded')
	    USdat = single(Eall);
	    Eall = single(Eall);
	    size(USdat)
	    
	%     for espid = espids
	%     for espid = 1:length(espiritfiles)
	%         espiritfile = espiritfiles{espid};
	    for espid = 1:length(espirit_cmds)
	        espirit_cmd = espirit_cmds{espid};
	        espiritfile = strrep(espirit_cmd,' ','_')
	        if run_espirit == 1
	    %         % --------------------------------- LOAD ESPIRiT
	    %         [ out_dir espiritfile ]
	    %         load([ out_dir espiritfile ]);
	    %         sensmaps = bart('flip 6', sensmaps);
	    %         sensmaps = sensmaps(recon_slices,:,:,:,:,:);
	    %         disp('ESPIRiT sensitivities (sensmaps)')
	    %         size(sensmaps)

	%             % --------------------------------- CALCULATE ESPIRiT SLICEWISE
	            % --- Use E1
	            Eksp = Eall(:,:,:,:,:,1); size(Eksp)
	%             % --- Use E1-2avg
	%             Eksp = (Eall(:,:,:,:,:,1) + Eall(:,:,:,:,:,2))/2; size(Eksp) % avg E1 and E2
	%             % --- Use EALLavg
	%             Eksp = sum(Eall,6) / 48; size(Eksp) % avg all
	%             % --- Eall_nzmean
	%             Eksp = sum(Eall,6)./(sum(logical(abs(Eall)),6)+0.000001); % WORKING
	%             Eksp = zeros(size(Eall,1), size(Eall,2), size(Eall,3), size(Eall,4));
	%             % Eksp = sum(Eall,6) ./ abs(sum(Eall~=0,6)); size(Eksp)
	%             % Eksp = Eksp + (sum(Eall,6) ./ sum(Eall~=0,6)); size(Eksp)
	%             Eksp = Eksp + (sum(Eall,6)); size(Eksp)
	%             ksp = sum(Eall,6) ./ sum(Eall~=0,6); % BROKEN?

	%             Eksp = bart('resize -c 1 118 2 72', bart('resize -c 1 20 2 10', Eksp));
	%             Eksp = bart('resize -c 1 240 2 20', bart('resize -c 1 16 2 16', Eksp));
	%             Eksp = bart('resize -c 1 249 2 16', bart('resize -c 1 14 2 16', Eksp));
	%             Eksp = bart('resize -c 1 249 2 32', bart('resize -c 1 14 2 16', Eksp));
	%             Eksp = bart('fft -u 4', bart('fft -u -i 4', Eksp)); % FAKE k-space in slice direction with fft
	%             Eksp = bart('fft 4', bart('fft -i 4', Eksp)); % FAKE k-space in slice direction with fft
	            % ---
	%             Eksp = bart('resize -c 1 20 2 10', Eksp); 
	 
	 
	            [ sensmaps espirit_eigenvals ] = bart(espirit_cmd, Eksp(1,:,:,:,:));
	%             evalc('[ sensmaps espirit_eigenvals ] = bart(espirit_cmd, Eksp(1,:,:,:,:));'); % quiet version
	            for x_id = 2:length(recon_slices)
	    %             [ sensmaps(x_id,:,:,:,:) espirit_eigenvals ] = bart(espirit_cmd, Eksp(x_id,:,:,:,:));
	                evalc('[ sensmaps(x_id,:,:,:,:) espirit_eigenvals ] = bart(espirit_cmd, Eksp(x_id,:,:,:,:));');
	            end
	            
	%             size(sensmaps)
	%             sensmaps = bart('fft 4', bart('resize -c 1 118 2 72', bart('fft -i 4', sensmaps)));
	%             size(sensmaps)
	            
	            
	            % --------------------------------- CALCULATE SENSEMAPS WITH DIRECT FFT OR CALDIR METHOD
	            
	%             % --- Crop k-space if desired
	%             Eksp = bart('resize -c 1 118 2 72', bart('resize -c 1 8 2 8', Eksp));
	% %             Eksp = bart('resize -c 1 236 2 30', bart('resize -c 1 16 2 10', Eksp));
	% %             Eksp = bart('resize -c 1 236 2 30', bart('resize -c 1 16 2 15', Eksp));
	% %             Eksp = bart('resize -c 1 160 2 30', bart('resize -c 1 16 ', Eksp));
	%             % --- Calculate
	%             sensmaps = bart('fft -i -u 6', Eksp(:,:,:,:)) ./ bart('rss 8', bart('fft -i -u 6', Eksp(:,:,:,:)));
	% %             sensmaps = bart('fft -i -u 6', Eksp(:,:,:,:)) ./ abs(bart('rss 8', bart('fft -i -u 6', Eksp(:,:,:,:))));
	            
	%             % --- CALDIR
	%             Eksp = bart('fft -u 1', Eksp);
	% %             Eksp = bart('resize -c 0 42 1 160 2 30', bart('resize -c 0 8 1 8 2 8', Eksp));
	%             Eksp = bart('resize -c 1 160 2 30', bart('resize -c 1 8 2 8', Eksp));
	        
	%             sensmaps = bart('caldir 30:8:8 ', Eksp);
	            
	            toc
	            disp('Sensemaps generated')
	            
	            % SINGLE
	            sensmaps = single(sensmaps);
	            
	            size(sensmaps)

	            % ---------------------------------- VIEW MAPS
	            view_chan = 1:8; % 1:2:30
	            for map_set = 1:1:size(sensmaps,5)
	                count = 0;
	                figure()
	                for val = view_chan
	                    count = count + 1;
	                    subplot(2,4,count), imshow(squeeze(abs(sensmaps(view_slice,:,:,val,map_set))), [ ])
	                    title({'Map' + string(map_set) + ' Chan' + string(val) }) 
	                end
	                saveas(gcf,['ESPIRiT_maps_set' num2str(map_set) '_' espiritfile(end-14:end-4) '.png']);
	            end
	            % ---------------------------------
	        end

	        % --------------------------------- OUTPUT FILES
	%         R_ids = [ 1 ]; % 3:4  or 1:size(usmasks,5)  to do all
	%         mas_ids = [ 1 ]; % [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 ] or  1:size(USmasks,4)  to do all
	%         out_dir_appended = [ out_dir 'ProUS/']
	%         mkdir(out_dir_appended)
	%         out_textfile = [ out_dir 'Output_CALIPR_ProUS_Cord_ECCENTRIC_DEV10.txt']; 
	        USmasks_str = [ dat_name '_bPW' num2str(use_bartPW) '_zf' ZeroFill '_' sprintf('%s%u',calibrate_dat,calibrate_datamatch) espiritfile ];
	        out_info = [];
	        

	        % ----------------------------------------------------------------------------------------------------------

	        % ------- If applying mask to remove noise/background regions:
	        if apply_mat_mask == 1; FSref = FSref .* maskbrain(:,:,:); end
	        % ------- Create normalized versions for error calculation
	        % Normalize all echoes together because relative intensities matter for quantitative comparison (not just individually scaled echo images)
	        if calc_error == 1; evalc('FSref_norm = bart(''normalize 15'', FSref);'); end
	        % --------------------------- Prep ESPIRiT sensmaps
	        disp('Loading sensemaps')
	        if espirit_use_sets == 1
	            FSsensmaps = sensmaps(:,:,:,:,:); % use multiple sets
	        else
	            FSsensmaps = sensmaps(:,:,:,:,1); % use one set only
	        end
	        % put slices in SLICE dimension 14, echoes in TE (6th)
	        FSsensmaps = permute(FSsensmaps, [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]); % size(FSsensmaps)

	        batchsize
	        nbatches = length(recon_slices)/batchsize
	        
	        % ------------------------------------------------ BEGINLOOP ------------------------------------------------
	        % over the different acceleration factors
	        for R_id = R_ids % 1 % 3:4  or 1:size(USmasks,5)  to do all
	            cd(out_dir_appended);
	            % over the 10 different temporal sampling schemes
	            for mas_id = mas_ids % [ 1 ] % [ 25 24 23 22 21 30 29 28 27 26 ]  or  1:size(USmasks,4)  to do all

	%                 cmd_mas = string(USmasks_str) + '_mID' + string(mas_id) + '_rID' + R_id
	                cmd_mas = string(USmasks_str)
	                mkdir( out_dir_appended );
	                mkdir( out_dir_appended + string(cmd_mas) );
	                cd( out_dir_appended + string(cmd_mas) );
	                
	                % get accel info
	%                 Rfactor = size(USdat,1)*size(USdat,2)/2*size(USdat,3)/2*size(USdat,4)*size(USdat,6)/nnz(USdat); % account for zero-fill
	                Rfactor = size(USdat,1)*114*72*size(USdat,4)*size(USdat,6)/nnz(USdat); % account for zero-fill
	                acqtime = datestr(seconds(114*72*1.252/Rfactor),'MM:SS ') % TR1277->10245 TR1252->10045
	                Rfactor = num2str(Rfactor,3)

	                % ----- Generate sampling image before clear (from ACTUAL k-space data, after masking)
	                show = 1:12; figure(),montage(squeeze(abs(USdat(1,:,:,1,:,show)).^0.0001),'Indices', show, 'size',[2 length(show)/2], 'DisplayRange', [0 0.5]);
	    %             show = 1:12; figure(),montage(squeeze(USmask(1,:,:,1,1,:)),'Indices', show, 'size',[2 length(show)/2], 'DisplayRange', [0 0.5]);
	                title({'USmask Mid:' + string(mas_id) + '  R: ' + string(Rfactor)})
	                print('USmask.png', '-dpng', '-r900');

	                % FFT ZF recon comparison for US data
	%                 USfft = squeeze(bart('rss 8',bart('fft -u -i 6',USdat)));
	                
	                disp('BATCH USfft generation')
	                for batch_id = 0:nbatches-1
	%                     USfft(1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1,:,:,:) = squeeze(bart('rss 8',bart('fft -u -i 6',USdat(1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1,:,:,:,:,:))));
	                    % QUIET+SINGLE VERSION
	                    evalc('USfft(1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1,:,:,:) = single(squeeze(bart(''rss 8'',bart(''fft -u -i 6'',USdat(1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1,:,:,:,:,:)))));');
	                end
	                
	%                 % FOR ONE SLICE ONLY
	%                 USfft = permute(USfft, [ 4 1 2 3 ]);
	%                 save('USfft.mat', 'USfft');

	                evalc('USfft_norm = bart(''normalize 15'', USfft);');
	                

	                % put slices in SLICE dimension 14, echoes in TE (6th)
	                USdat = permute(USdat, [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]); % size(USdat)
	        %         FSsensmaps = permute(FSsensmaps, [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]); % size(FSsensmaps)
	                % ------------------ END NOT RUNNING ------------------


	                cd( out_dir_appended + string(cmd_mas) );
	                cmd_pics = cs_cmd_pics;
	%                     cmd_pics = cell_cmd_pics{:};
	                cmd_pics_clean = strrep(cmd_pics,' ','_')

	                out_info = [out_info,...
	                    '\n-----------------------------------------------------------------------------------------', ...
	                    '\n---------------------------------- FOR PICS CMD, MASK ------------------------------------', ...
	                    '\n-----------------------------------------------------------------------------------------', ...
	                    '\nRecon Command:                   ',cmd_pics,...
	                    '\nUndersampling Mask:              ',cmd_mas{1},...
	                    '\nUndersampling ID:                ',num2str(mas_id),...
	                    '\nAcquisition Time:                ',acqtime,...
	                    'with acceleration factor R=',Rfactor,...
	                    '\nSensitivities:                   ',num2str(size(FSsensmaps,5)),' | ',num2str(espirit_combine_sets),'      #Sets | Combined(0/1)',...
	                    '\nImage Matlab Mask Info:          ',num2str(apply_mat_mask),' | ',num2str(thresh),' | ',num2str(mask_borderpad),'      Apply | ThreshFrac | BorderPad',...
	                    ];

	                K = 0;
	                % ------------------------------------------------ START CS RECON ------------------------------------------------
	                if run_CSrecon == 1
	                    cd( out_dir_appended + string(cmd_mas) );

	                    K = 0;
	                    % directory for K value
	                    RECONdirstr = string(string(cmd_pics_clean) + 'K_' + string(K));
	                    mkdir(RECONdirstr);
	                    cd(RECONdirstr);
	                    
	%                     % -------- PI ALL VERSION
	%                     USdat(:,2:2:end,:,:,:,:) = 0;
	%                     % -------- CROPPED VERSION
	%                     USdat = bart('resize -c 1 240 2 20',bart('resize -c 1 16 2 10', USdat));
	%                     USdat = bart('resize -c 1 104 2 64',bart('resize -c 1 52 2 32', USdat));
	%                     USdat = bart('crop 1 64',bart('crop 2 40', Eall));
	%                     FSsensmaps = bart('crop 1 64',bart('crop 2 40', FSsensmaps));
	%                     USdat = bart('resize -c 1 128 2 80',bart('resize -c 1 64 2 40', Eall));

	%                     disp('US kspace Data        (USdat)')
	%                     size(USdat)

	%                     pics_mas_8ch_noSqueeze = bart('fft -i 6',USdat);
	%                     pics_mas_8ch_noSqueeze = bart('rss 8', FSsensmaps.*pics_mas_8ch_noSqueeze);
	%                     pics_mas_8ch_noSqueeze = bart('fmac', FSsensmaps, pics_mas_8ch_noSqueeze);
	%                     pics_mas_8ch_noSqueeze = bart('rss 8', bart('fft -i 6',USdat));

	                    % ------------------ START NOT RUNNING ------------------
	% %                         disp('Size of inputs:')
	% %                         size(USdat)
	% %                         size(FSsensmaps)
	%                     out_pics_mas_8ch = evalc( 'pics_mas_8ch_noSqueeze = bart(''' + string(cmd_pics) + ''',USdat, FSsensmaps);' )
	%     %                 disp('Size of output'); size(pics_mas_8ch_noSqueeze)

	                    % --- BATCH MODE
	                    
	%                     pics_mas_8ch_noSqueeze = permute( zeros(size(Eall,1),size(Eall,2),size(Eall,3),1,size(FSsensmaps,5),size(Eall,6)), [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]);
	                    % SINGLE VERSION
	                    pics_mas_8ch_noSqueeze = single(permute( zeros(size(Eall,1),size(Eall,2),size(Eall,3),1,size(FSsensmaps,5),size(Eall,6)), [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]));
	                    
	                    
	%                     size(pics_mas_8ch_noSqueeze)
	                    
	                    % --- BATCH MODE
	                    disp('BATCH UScs generation')
	                    
	                    for batch_id = 0:nbatches-1
	%                         disp('BATCH:')
	%                         batch_id
	                        
	%                         batch_USdat = USdat(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1);
	%                         batch_FSsensmaps = FSsensmaps(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1);
	%                         size(batch_USdat)
	%                         size(batch_FSsensmaps)
	                        
	                        
	% %                         size(pics_mas_8ch_noSqueeze(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1))
	%                         out_pics_mas_8ch = evalc( 'pics_mas_8ch_noSqueeze(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1) = bart(''' + string(cmd_pics) + ''',batch_USdat, batch_FSsensmaps);' )
	%                         out_pics_mas_8ch = evalc( 'batch_pics_mas_8ch_noSqueeze = bart(''' + string(cmd_pics) + ''',batch_USdat, batch_FSsensmaps);' )
	%                         pics_mas_8ch_noSqueeze(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1) = batch_pics_mas_8ch_noSqueeze;
	%                         clear batch_pics_mas_8ch_noSqueeze
	                        
	%                         out_pics_mas_8ch = evalc( 'pics_mas_8ch_noSqueeze(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1) = bart(''' + string(cmd_pics) + ''',USdat(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1), FSsensmaps(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1));' )
	                        % QUIET AND SINGLE VERSION
	                        out_pics_mas_8ch = evalc( 'pics_mas_8ch_noSqueeze(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1) = single(bart(''' + string(cmd_pics) + ''',USdat(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1), FSsensmaps(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1)));' );
	                        
	%                         out_pics_mas_8ch = evalc( 'tmp = bart(''' + string(cmd_pics) + ''',USdat(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1), FSsensmaps(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1));' )
	%                         size(tmp)
	                    end


	                    % ovrwite CSrecon specific changes 
	%                     USdat = Eall;

	                    % put slices back in READOUT (1) from SLICE dimension (14)
	                    pics_mas_8ch_noSqueeze = permute(pics_mas_8ch_noSqueeze, [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]); % size(pics_mas_8ch_noSqueeze)

	                    % ------  For multiple ESPIRiT maps, need to combine or throw out others
	                    if (espirit_use_sets == 1) && (espirit_combine_sets == 1)
	                        pics_mas_8ch_noSqueeze = bart('rss 16', pics_mas_8ch_noSqueeze); % combine multiple sets along ESPIRiT dimension
	                    else
	                        pics_mas_8ch_noSqueeze = pics_mas_8ch_noSqueeze(:,:,:,:,1,:,:,:,:); % just take first (or only) set
	                    end

	                    % Squeeze back to convenient size
	                    UScs = squeeze(pics_mas_8ch_noSqueeze); % size(pics_mas_8ch)
	                    clear pics_mas_8ch_noSqueeze 

	                    save('UScs.mat', 'UScs', '-v7.3', '-nocompression');
	                    
	                    toc
	                    disp('Done CS Recon')

	                    USdat = permute(Eall, [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]);

	                    % ------------------ END NOT RUNNING ------------------

	                    % ------------------ SAVE CS output as par/rec and nii, if desired
	                    if save_cs_nii == 1
	                        disp('Generating and Saving Par/Rec Files')

	                        % ----- IF NOT RUNNING
	                        load('UScs.mat');
	                        % UScs = abs(UScs);

	                        % ----------------- CS
	                        clear tmpMR
	                        tmpMR = MRecon(dat_labfile, dat_rawfile);
	                        tmpMR.Parameter.Parameter2Read.typ = [ 1 ]; % read in standard data (1) phase correction data (3) noise data (5)
	                        tmpMR.Parameter.Parameter2Read.Update;

	                        tmpMR.Data = bart('flip 6', abs(UScs)); % 
	                        tmpMR.Data = permute( tmpMR.Data, [ 1 2 3 5 6 7 4 ]);

	                        tmpMR.Parameter.ReconFlags.isread = 1;
	                        tmpMR.Parameter.ReconFlags.issorted = 1;
	                        tmpMR.Parameter.ReconFlags.ispartialfourier = [0 0];
	                        tmpMR.Parameter.ReconFlags.isgridded = 0;
	                        tmpMR.Parameter.ReconFlags.isimspace = [1 1 1];
	                        tmpMR.Parameter.ReconFlags.iscombined = 1;
	                        tmpMR.Parameter.ReconFlags.isoversampled = [0 0 1];
	                        tmpMR.Parameter.ReconFlags.isreadparameter = 1;
	                        tmpMR.Parameter.ReconFlags.israndphasecorr = 0;
	                        tmpMR.Parameter.ReconFlags.ispdacorr = 1;
	                        tmpMR.Parameter.ReconFlags.isdcoffsetcorr = 1;
	                        tmpMR.Parameter.ReconFlags.isdepicorr = 0;
	                        tmpMR.Parameter.ReconFlags.ismeasphasecorr = 0;
	                        tmpMR.Parameter.ReconFlags.isnonlincorr = 0;
	                        tmpMR.Parameter.ReconFlags.isunfolded = 0;
	                        tmpMR.Parameter.ReconFlags.iszerofilled = [0 0];
	                        tmpMR.Parameter.ReconFlags.isrotated = 0;
	                        tmpMR.Parameter.ReconFlags.isconcomcorrected = 0;
	                        tmpMR.Parameter.ReconFlags.isgeocorrected = 0;
	                        tmpMR.Parameter.ReconFlags.issegmentsdivided = 0;
	                        tmpMR.Parameter.ReconFlags.isecc = 0;
	                        tmpMR.Parameter.ReconFlags.isaveraged = 1;

	                        tmpMR.GridderNormalization;
	                        tmpMR.SENSEUnfold;
	                        tmpMR.PartialFourier;
	                        tmpMR.ConcomitantFieldCorrection;
	                        tmpMR.DivideFlowSegments;
	                        tmpMR.Average;
	                        tmpMR.GeometryCorrection;
	                        tmpMR.RemoveOversampling;
	                        tmpMR.ZeroFill;
	                        tmpMR.RotateImage;

	                        tmpMR.WritePar( [ 'CS.par' ] );
	                        tmpMR.WriteRec( [ 'CS.rec' ] );
	                        !dcm2niix -b n -m y -s y -p y -v n -z y -w 1 -f CS CS.par >/dev/null
	                        !fslmerge -t CS.nii.gz $( for i in {1..56}; do printf "CS_e${i}.nii.gz "; done )
	                        !mv CS_e1.nii.gz CS_E1.nii.gz
	                        !rm CS_e*.nii.gz

	                        system(sprintf("antsBrainExtraction.sh -d 3 -k 0 -z 0 -a CS_E1.nii.gz -e %s/T_template0.nii.gz -m %s/T_template0_BrainCerebellumProbabilityMask.nii.gz -f %s/T_template0_BrainCerebellumRegistrationMask.nii.gz -o CS_E1 &>/dev/null",oasis_template_path,oasis_template_path,oasis_template_path))
	                        !gunzip CS.nii.gz
	                        !gunzip CS_E1BrainExtractionMask.nii.gz
	                        clear tmpMR
	                    end


	                    % ------------------------------------------------ VIEW+COMPARE
	                    % IMAGE VIEW
	                    for eid=echoes2show
	                        % normalize individual slice+echo for viewing with consistent window
	                        evalc('FSref_img = bart(''normalize 7'', squeeze(abs(FSref(view_slice,:,:,eid))));');
	                        evalc('USfft_img = bart(''normalize 7'', squeeze(abs(USfft(view_slice,:,:,eid))));');
	                        evalc('UScs_img = bart(''normalize 7'', squeeze(abs(UScs(view_slice,:,:,eid))));');
	                        % diff and SSIM
	                        diff_UScs = abs(FSref_img - UScs_img); [UScs_E1_ssimval, UScs_E1_ssimmap] = ssim(100*UScs_img,100*FSref_img); 
	                        % phase image, scaled so -pi->pi to fill imwindow, or unwrapped version
	                        UScs_phase = (squeeze(angle(UScs(view_slice,:,:,eid)))).*((imwindow(2)-imwindow(1))/(2*pi)) + ((imwindow(2)-imwindow(1))/2); 
	                        % display
	                        figure(), imshow([ FSref_img UScs_img 10*diff_UScs UScs_E1_ssimmap.*max(imwindow) UScs_phase USfft_img ], imwindow,'InitialMagnification','fit'); 
	                        title({' FSref | UScs | 10x(FFSref-UScs) | SSIMmap | Phase UScs | USfft ', 'Acceleration R = ' + string(Rfactor) + '     Echo: ' + string(eid) + '    K: ' + string(K) + '  MaskID: ' + num2str(mas_id), ' PICS: ' + string(cmd_pics) })
	                        print(sprintf('UScs_E%u.png',eid), '-dpng', '-r600');
	                    end

	                    UScs_phase = (squeeze(angle(UScs(view_slice,:,:,1)))).*((scale_MWF(2)-scale_MWF(1))/(2*pi)) + ((scale_MWF(2)-scale_MWF(1))/2);
	                    evalc('USfft_img = (scale_MWF(2)/imwindow(2))*bart(''normalize 7'', squeeze(abs(USfft(view_slice,:,:,1))));');
	                    evalc('UScs_img1 = (scale_MWF(2)/imwindow(2))*bart(''normalize 7'', squeeze(abs(UScs(view_slice,:,:,1))));');
	                    evalc('UScs_img4 = (scale_MWF(2)/imwindow(2))*bart(''normalize 7'', squeeze(abs(UScs(view_slice,:,:,4))));');
	                    evalc('UScs_img10 = (scale_MWF(2)/imwindow(2))*bart(''normalize 7'', squeeze(abs(UScs(view_slice,:,:,10))));');
	                    evalc('UScs_img28 = (scale_MWF(2)/imwindow(2))*bart(''normalize 7'', squeeze(abs(UScs(view_slice,:,:,28))));');

	                    % view
	                    figure(), imshow([ USfft_img UScs_phase UScs_img1 UScs_img4 UScs_img10 UScs_img28 ], scale_MWF,'InitialMagnification','fit');
	                    title({'USfft | UScs Phase | UScs E1/E4/E10/E28 ', 'Acceleration R = ' + string(Rfactor) + '    K: ' + string(K) + '  MaskID: ' + num2str(mas_id), ' PICS: ' + string(cmd_pics) })
	                    print('UScs_Echoes', '-dpng', '-r600');
	                    
	                    
	                    % Axial views
	                    evalc('UScs_img1_1 = (scale_MWF(2)/imwindow(2))*bart(''normalize 7'', squeeze(abs(UScs(:,:,view_slices_ax(1)+ax_shift,1))));');
	                    evalc('UScs_img1_2 = (scale_MWF(2)/imwindow(2))*bart(''normalize 7'', squeeze(abs(UScs(:,:,view_slices_ax(2)+ax_shift,1))));');
	                    evalc('UScs_img1_3 = (scale_MWF(2)/imwindow(2))*bart(''normalize 7'', squeeze(abs(UScs(:,:,view_slices_ax(3)+ax_shift,1))));');   
	                    evalc('UScs_img1_4 = (scale_MWF(2)/imwindow(2))*bart(''normalize 7'', squeeze(abs(UScs(:,:,view_slices_ax(4)+ax_shift,1))));');
	                    evalc('UScs_img1_5 = (scale_MWF(2)/imwindow(2))*bart(''normalize 7'', squeeze(abs(UScs(:,:,view_slices_ax(5)+ax_shift,1))));'); 

	                    figure(), imshow([ UScs_img1_1; UScs_img1_2; UScs_img1_3; UScs_img1_4; UScs_img1_5 ], scale_MWF,'InitialMagnification','fit');
	                    title({'UScs E1 ', 'Acceleration R = ' + string(Rfactor) + '    K: ' + string(K) + '  MaskID: ' + num2str(mas_id), ' PICS: ' + string(cmd_pics) })
	                    print('UScs_ALL_AXIAL', '-dpng', '-r600');


	    %                 % ------------------------------------------------ APPEND TO OUT INFO
	    %                 out_info = [out_info,...
	    %                     '\n-----------------------------------------------------------------------------------------', ...
	    %                     '\nRecon Type:              ','COMPRESSED SENSING',...
	    %                     '\nRecon Command:           ',cmd_pics,...
	    %                     '\nTemporal Basis:          ',num2str(K),...
	    %                     ];
	    %                 % add error metrics
	    %                 if (calc_error == 1)
	    %                     if analyseMWI == 0
	    %                         out_info = [out_info,...
	    %                             '\n\nError Metrics:          < Eall >   < E1 >   < E13 >  < E28 >   < E56 >',...
	    %                             '\nNRMSE:                   ',num2str(nrmse(1),'%.4f'),'    ',num2str(nrmse(2),'%.4f'),'   ',num2str(nrmse(3),'%.4f'),'   ',num2str(nrmse(4),'%.4f'),'    ',num2str(nrmse(5),'%.4f'),...
	    %                             '\nSSIM:                    ','  N/A ','    ',num2str(ssimval(2),'%.4f'),'   ',num2str(ssimval(3),'%.4f'),'   ',num2str(ssimval(4),'%.4f'),'    ',num2str(ssimval(5),'%.4f'),...
	    %                             ];
	    %                     else
	    %                         out_info = [out_info,...
	    %                             '\n\nError Metrics:          < Eall >   < E1 >   < E13 >  < E28 >   < E56 >    < MWF >   < ALPHA >',...
	    %                             '\nNRMSE:                   ',num2str(nrmse(1),'%.4f'),'    ',num2str(nrmse(2),'%.4f'),'   ',num2str(nrmse(3),'%.4f'),'   ',num2str(nrmse(4),'%.4f'),'    ',num2str(nrmse(5),'%.4f'),'     ',num2str(nrmse(6),'%.4f'),'     ',num2str(nrmse(7),'%.4f'),...
	    %                             '\nSSIM:                    ','  N/A ','    ',num2str(ssimval(2),'%.4f'),'   ',num2str(ssimval(3),'%.4f'),'   ',num2str(ssimval(4),'%.4f'),'    ',num2str(ssimval(5),'%.4f'),'     ',num2str(ssimval(6),'%.4f'),'     ',num2str(ssimval(7),'%.4f'),...
	    %                             '\nStatistics:             < Mean >   < Med >  < Max >  < Std >   < Var >    < IQR >  ',...
	    %                             '\nMWF (Ref, No 0):         ',num2str(mean(nonzeros(FSref_MWF)),'%.4f'),'    ',num2str(median(nonzeros(FSref_MWF)),'%.4f'),'   ',num2str(max(nonzeros(FSref_MWF)),'%.4f'),'   ',num2str(std(nonzeros(FSref_MWF)),'%.4f'),'    ',num2str(var(nonzeros(FSref_MWF)),'%.4f'),'     ',num2str(iqr(nonzeros(FSref_MWF)),'%.4f'),' ',...
	    %                             '\nMWF (US, No 0):          ',num2str(mean(nonzeros(UScs_MWF)),'%.4f'),'    ',num2str(median(nonzeros(UScs_MWF)),'%.4f'),'   ',num2str(max(nonzeros(UScs_MWF)),'%.4f'),'   ',num2str(std(nonzeros(UScs_MWF)),'%.4f'),'    ',num2str(var(nonzeros(UScs_MWF)),'%.4f'),'     ',num2str(iqr(nonzeros(UScs_MWF)),'%.4f'),' ',...
	    %                             '\nALPHA (Ref, No 0/180):   ',num2str(mean(ref),'%.2f'),'    ',num2str(median(ref),'%.2f'),'   ',num2str(max(ref),'%.2f'),'   ',num2str(std(ref),'%.4f'),'    ',num2str(var(ref),'%.3f'),'     ',num2str(iqr(ref),'%.4f'),' ',...
	    %                             '\nALPHA (US, No 0/180):    ',num2str(mean(proj),'%.2f'),'    ',num2str(median(proj),'%.2f'),'   ',num2str(max(proj),'%.2f'),'   ',num2str(std(proj),'%.4f'),'    ',num2str(var(proj),'%.3f'),'     ',num2str(iqr(proj),'%.4f'),' ',...
	    %                             ];
	    %                     end
	    %                 end

	                end
	                % ------------------------------------------------ END CS RECON ------------------------------------------------            
	                
	                % ------------------------------------------------------------------------
	                % ------------------ START CALIPR BASIS SET CALIBRATION ------------------
	                if run_CALIPRrecon == 1
	                    disp('Beginning CALIPR Basis Set Calibration ')
	                    cd( out_dir_appended + string(cmd_mas) );
	                    cmd_pics = cs_cmd_pics;
	%                         cmd_pics = cell_cmd_pics{:};
	                    cmd_pics_clean = strrep(cmd_pics,' ','_')
	                    % directory for K value
	                    RECONdirstr = string(string(cmd_pics_clean) + 'K_' + string(K));
	                    mkdir(RECONdirstr); cd(RECONdirstr);

	                    if strcmp(calibrate_dat, 'UScs')
	                        % ------------------------------------------------------------- USE PREVIOUSLY RECONSTRUCTED UScs (full resolution)
	                        Dat = UScs;
	%                         Dat = UScs(:,70-20:70+20,:,:);
	%                         Dat = bart('homodyne -I 0 1', UScs);
	                        
	%                         % ----------- N4 Bias Field Correction
	%                         UStmp = Dat;

	%                         % UStmp_nii = make_nii(double(abs(UStmp(:,:,:,1)))); 
	%                         UStmp_nii = make_nii(double(abs(UStmp(:,:,:,10)))); 
	%                         % size(UStmp_nii.img)
	%                         save_nii(UStmp_nii, 'UStmp_nii.nii');

	%                         % mask_nii = make_nii(abs(thresh_mask)); 
	%                         % mask_nii = make_nii(repmat(abs(thresh_mask), [ 1 1 1 56 ])); 
	%                         % size(mask_nii.img)
	%                         % save_nii(mask_nii, 'mask_nii.nii');

	%                         % NOTE: needed to remove link to matlab version of /usr/local/MATLAB/R2019b/sys/os/glnxa64/libstdc++.so.6    to force use of system version used for ANTs installation
	%                         % TO RUN:
	%                         system('N4BiasFieldCorrection -d 3 -i UStmp_nii.nii -o [UStmp_nii_N4.nii,UStmp_nii_N4BiasField.nii] ');

	%                         UStmp_nii_N4BiasField = load_untouch_nii('UStmp_nii_N4BiasField.nii');
	%                         UStmp_N4BiasField = UStmp_nii_N4BiasField.img;
	%                         % size(UStmp_N4BiasField)

	%                         % Divide each image by bias field to get corrected version
	%                         Dat = UScs ./ UStmp_N4BiasField;
	%                         % UScalipr = UStmp;
	%                         % size(UStmp)
	%                         % -----------
	                        
	                        
	                        % -------------------------------------------------------------
	                    else % all others use tinydat   
	                        % ------------------------------------------------------------- CREATE tinyUSdat, option to create tinyFSsensmaps, option to crop for matching data across echo train
	%                         tinyUSdat = bart('resize -c 1 8 2 8', bart('resize -c 1 8 2 8', permute(USdat, [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]))); % resize twice for zero-fill, if desired
	%                         tinyUSdat = bart('resize -c 1 118 2 15', bart('resize -c 1 118 2 15', permute(USdat, [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]))); % resize twice for zero-fill, if desired
	                        tinyUSdat = permute(USdat, [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]);
	                        
	                        if (strcmp(calibrate_dat, 'tinyUScs') || strcmp(calibrate_dat, 'tinyUSfft_fmac'))
	                            % --------------------------------- CALCULATE ESPIRiT SLICEWISE
	                            tinyEksp = tinyUSdat(:,:,:,:,:,1);
	                            % Eksp = (Eall(:,:,:,:,:,1) + Eall(:,:,:,:,:,2))/2; % use average of first two echoes
	                            tinyFSsensmaps = bart(tinyespirit_cmd, tinyEksp(1,:,:,:,:));
	%                             evalc('tinyFSsensmaps = bart(tinyespirit_cmd, tinyEksp(1,:,:,:,:));'); % quiet version
	                            for x_id = 2:length(recon_slices)
	                            %             [ sensmaps(x_id,:,:,:,:) espirit_eigenvals ] = bart(tinyespirit_cmd, Eksp(x_id,:,:,:,:));
	                                evalc('tinyFSsensmaps(x_id,:,:,:,:) = bart(tinyespirit_cmd, tinyEksp(x_id,:,:,:,:));');
	                            end
	                            % ---------------------------------- VIEW MAPS
	                            view_chan = 1:8; % 1:2:30
	                            for map_set = 1:1:size(tinyFSsensmaps,5)
	                                count = 0;
	                                figure()
	                                for val = view_chan
	                                    count = count + 1;
	                                    subplot(2,4,count), imshow(squeeze(abs(tinyFSsensmaps(view_slice,:,:,val,map_set))), [ ])
	                                    title({'Map' + string(map_set) + ' Chan' + string(val) }) 
	                                end
	                                saveas(gcf,['tinyESPIRiT_maps_set' num2str(map_set) '_' espiritfile(end-14:end-4) '.png']);
	                            end
	                        end
	                        if calibrate_datamatch
	                            % --------------------------------- CROP FOR MATCHING DATAPOINTS ACROSS ECHOES
	                            tic
	                            datmatch = ones(1, size(tinyUSdat,2), size(tinyUSdat,3));
	                            for eid = 1:size(tinyUSdat,6)
	                                tmpdat = tinyUSdat(size(tinyUSdat,1)/2,:,:,size(tinyUSdat,4)/2,1,eid);
	                            %     indices = find(tmpdat==0);
	                            %     datmatch(indices)=0;
	                                datmatch(find(tmpdat==0))=0;
	                            end
	                            disp(sprintf('Before calibration data matching across echoes:  %u', nnz(tinyUSdat(:))))
	                            tinyUSdat = tinyUSdat(:,:,:,:,1,:) .* datmatch(1,:,:,1,1,1); % mask to only include datapoints acquired at all echoes
	                            disp(sprintf('After calibration  data matching across echoes:  %u', nnz(tinyUSdat(:))))
	                            % ---------------------------------
	                       end
	                       
	%                        % To match PI for all:
	%                         tinyUSdat(:,1:2:end,:,:,:,:) = 0;
	                        
	                        figure, imshow( [ squeeze(abs(tinyUSdat(2,:,:,4,1,1))).^0.1 squeeze(abs(tinyUSdat(2,:,:,4,1,2))).^0.1 squeeze(abs(tinyUSdat(2,:,:,4,1,3))).^0.1 squeeze(abs(tinyUSdat(2,:,:,4,1,4))).^0.1] , [],'InitialMagnification','fit')
	                        title('tinyUSdat')
	                        % ------------------------------------------------------------- 

	                        if strcmp(calibrate_dat, 'tinyUScs')
	                        
	                            % ------------------------------------------------------------- USE tinyUScs FOR CALIBRATION: perform pics recon (requires tinyFSsensmaps)
	                            tinyUSdat = permute(tinyUSdat, [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]);
	                            tinyFSsensmaps = permute(tinyFSsensmaps, [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]);
	                            evalc( 'tinypics_mas_8ch_noSqueeze = bart(''' + string(tinycmd_pics) + ''',tinyUSdat, tinyFSsensmaps);' )
	                            tinyUSdat = permute(tinyUSdat, [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]);
	                            tinyFSsensmaps = permute(tinyFSsensmaps, [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]);
	                            tinypics_mas_8ch_noSqueeze = permute(tinypics_mas_8ch_noSqueeze, [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]);

	                            tinyUScs = squeeze(tinypics_mas_8ch_noSqueeze(:,:,:,:,1,:,:,:,:)); % just take first (or only) set
	                            Dat = tinyUScs;
	                            figure, imshow( [ squeeze(abs(Dat(2,:,:,1))) squeeze(abs(Dat(2,:,:,2))) squeeze(abs(Dat(2,:,:,3))) squeeze(abs(Dat(2,:,:,4))) squeeze(abs(Dat(2,:,:,5))) squeeze(abs(Dat(2,:,:,6))) squeeze(abs(Dat(2,:,:,7))) squeeze(abs(Dat(2,:,:,8))) ], [],'InitialMagnification','fit')
	                            title(calibrate_dat)
	                            print(sprintf('calibration_%s%u',calibrate_dat,calibrate_datamatch), '-dpng', '-r600');
	                            % -------------------------------------------------------------
	                        elseif strcmp(calibrate_dat, 'tinyUSfft_rss')
	                            % ------------------------------------------------------------- USE tinyUSfft_rss FOR CALIBRATION
	                            tinyUSfft_rss = squeeze(bart('rss 8', bart('fft -i -u 6', tinyUSdat)));
	                            Dat = tinyUSfft_rss;
	                            figure, imshow( [ squeeze(abs(Dat(2,:,:,1))) squeeze(abs(Dat(2,:,:,2))) squeeze(abs(Dat(2,:,:,3))) squeeze(abs(Dat(2,:,:,4))) squeeze(abs(Dat(2,:,:,5))) squeeze(abs(Dat(2,:,:,6))) squeeze(abs(Dat(2,:,:,7))) squeeze(abs(Dat(2,:,:,8))) ], [],'InitialMagnification','fit')
	                            title(calibrate_dat)
	                            print(sprintf('calibration_%s%u',calibrate_dat,calibrate_datamatch), '-dpng', '-r600');
	                            % -------------------------------------------------------------
	                        elseif strcmp(calibrate_dat, 'tinyUSfft_fmac')
	                            % ------------------------------------------------------------- USE tinyUSfft_fmac FOR CALIBRATION (requires tinyFSsensmaps)
	                            % % tinyUSfft = bart('fmac -s 8 -C',bart('fft -i -u 6',tinyUSdat), tinysensmaps);
	                            % tinyUSfft = squeeze(bart('fmac -s 8 -C',bart('fft -i -u 6',tinyUSdat(:,:,:,:,1,:)), tinyFSsensmaps(:,:,:,:,1,1)));
	                            tinyUSfft_fmac = bart('fmac -s 8 -C',bart('fft -i -u 6',tinyUSdat(:,:,:,:,1,:)), tinyFSsensmaps(:,:,:,:,:,1)); % allow for multiple maps
	                            Dat = squeeze(tinyUSfft_fmac(:,:,:,:,1,:)); % take first map only
	                            figure, imshow( [ squeeze(abs(Dat(2,:,:,1))) squeeze(abs(Dat(2,:,:,2))) squeeze(abs(Dat(2,:,:,3))) squeeze(abs(Dat(2,:,:,4))) squeeze(abs(Dat(2,:,:,5))) squeeze(abs(Dat(2,:,:,6))) squeeze(abs(Dat(2,:,:,7))) squeeze(abs(Dat(2,:,:,8))) ], [],'InitialMagnification','fit')
	                            title(calibrate_dat)
	                            print(sprintf('calibration_%s%u',calibrate_dat,calibrate_datamatch), '-dpng', '-r600');
	                            % -------------------------------------------------------------
	                        end
	                        UScs = USfft; % for including in figures later, to avoid error if no UScs ran
	                    end
	                    % REPLACE WITH INTENSITY THRESHOLD
	    %                 Dat = Dat .* maskbrain;

	                    % ----------------- CREATE SIGNAL EVOLUTION MATRIX with column of all echo times, column of all voxels 
	                    clear iflat_ALL iflat sig_evol_ALL sig_evol sig_evol_disp datU
	                    for echo = 1:size(UScs,4)
	                        clear tmpdat
	    %                     evalc('tmp_i_flat = bart(''flatten'', Dat(:,:,echo));'); % flatten
	                        evalc('tmp_i_flat = bart(''flatten'', Dat(:,:,:,echo));'); % flatten
	                        iflat_ALL(:,echo) = tmp_i_flat; % add to all sig evol
	                    end
	                    % ----------------- 
	                    % ----------------- remove threshold voxels
	    %                 disp('Threshold:')
	    %                 max(abs(iflat_ALL(:,1)))*thresh
	    %                 size(iflat_ALL)
	                    % iflat_ALL = iflat_ALL(iflat_ALL(:,1) >= 200, :);
	                    iflat_ALL = iflat_ALL(abs(iflat_ALL(:,1)) >= max(abs(iflat_ALL(:,1)))*thresh, :); % size(iflat_ALL)
	                    iflat = iflat_ALL;
	                    sig_evol_ALL = permute(iflat_ALL, [ 2 1 ]); % rearrange for SVD

	                    thresh_img = Dat(:,:,:,1);    
	                    thresh_img(abs(thresh_img) < max(abs(iflat_ALL(:,1)))*thresh) = 0; % size(thresh_img)
	                    thresh_mask = thresh_img;
	                    thresh_mask(abs(thresh_mask) >= max(abs(iflat_ALL(:,1)))*thresh) = 1;

	                    if show_threshold == 1
	                        figure(), imshow([ squeeze(abs(Dat(view_slice,:,:,1))) squeeze(abs(Dat(view_slice,:,:,2))) squeeze(abs(Dat(view_slice,:,:,3))) 1000*squeeze(abs(thresh_img(view_slice,:,:))) ], [ 0 1*max(max(abs(Dat(view_slice,:,:,1)))) ],'InitialMagnification','fit');
	                        title('Calibration Data Echo 1 | Echo 2 | Echo 3 | Threshold Mask')
	                        print('calibration_threshold_mask', '-dpng', '-r300');
	                    end

	                    % ----------------- if more than maxN are given, randomly choose
	                    if length(iflat) > maxN
	                        clear sig_evol
	                        idx = randperm(length(iflat));
	                        iflat_maxN = iflat(idx(1:maxN),:);
	                        % maxN and ALL are separate
	                        sig_evol = permute(iflat_maxN, [ 2 1 ]); % maxN and ALL are different
	                        disp('Length sig_evol > maxN)')
	                    else
	                        sig_evol = permute(iflat_ALL, [ 2 1 ]); % maxN and ALL are the same
	                        disp('Length sig_evol < maxN)')
	                    end

	                    % ----------------- smaller number just for display
	                    % randomly choose if more than N are given
	                    idx_disp = randperm(length(iflat));
	                    iflat_maxN_disp = iflat(idx_disp(1:maxN_disp),:);
	                    sig_evol_disp = permute(iflat_maxN_disp, [ 2 1 ]);

	                    % ----------------- SVD
	                    % doc svd()
	                    % X is size (echo train length , number of signal evolutions)
	                    [datU_mag, ~, ~] = svd(abs(sig_evol), 'econ'); % COMPLEX MAGNITUDE
	                    [datU_real, ~, ~] = svd(real(sig_evol), 'econ'); % REAL
	                    [datU_imag, ~, ~] = svd(imag(sig_evol), 'econ'); % IMAG
	                    [datU_cpx, ~, ~] = svd(sig_evol, 'econ'); % COMPLEX
	                    datU_comb = complex(datU_real, datU_imag); % COMBINED REAL+IMAG (COMPLEX)
	                    [datU_apnd, ~, ~] = svd([ real(sig_evol) imag(sig_evol) ], 'econ'); % APPEND REAL AND IMAG SIGNALS

	                    % use:
	    %                 datU = datU_mag;
	    %                 datU = datU_real;
	    %                 datU = datU_imag;
	    %                 datU = datU_cpx;
	    %                 datU = datU_comb;
	    %                 datU = real(datU_cpx);

	                    datU_ALL = datU_mag;
	                    datU_ALL(:,:,2) = datU_real;
	                    datU_ALL(:,:,3) = datU_imag;
	                    datU_ALL(:,:,4) = datU_cpx;
	                    datU_ALL(:,:,5) = real(datU_cpx);
	                    datU_ALL(:,:,6) = abs(datU_cpx);
	                    datU_ALL(:,:,7) = datU_comb;
	                    datU_ALL(:,:,8) = datU_apnd;
	                    
	                    toc
	                    disp('Done CALIPR Basis Set Calibration ')

	                end
	                % ------------------ END CALIPR BASIS SET CALIBRATION ------------------
	                % ----------------------------------------------------------------------
	                
	                % ------------------------------------------------ BEGINLOOP ------------------------------------------------
	                for cell_cmd_pics = ALL_cmd_pics
	                    % ---------------------------------------------- START KVAL LOOP ----------------------------------------------
	                    for K = Kvals
	                        % ---------------------------------------------- START CALIPR RECON ----------------------------------------------
	                        if run_CALIPRrecon == 1
	                            disp('Beginning CALIPR Recon ')
	                            cd( out_dir_appended + string(cmd_mas) );
	        %                     cmd_pics = calipr_cmd_pics;
	                            cmd_pics = cell_cmd_pics{:};
	                            cmd_pics_clean = strrep(cmd_pics,' ','_')

	                            % to loop over bases
	                            for bid = basis_ids % 1:7
	                                cd( out_dir_appended + string(cmd_mas) );
	                                % directory for K value
	                                RECONdirstr = string(string(cmd_pics_clean) + 'K_' + string(K)) + '_CALIPR_B' + string(bid);
	                                mkdir(RECONdirstr); cd(RECONdirstr);
	        %                         % directory for K value
	        %                         RECONdirstr = string(string(cmd_pics_clean) + 'K_' + string(K)) + '_CALIPR';
	        %                         mkdir(RECONdirstr); cd(RECONdirstr);
	                                datU = datU_ALL(:,:,bid);


	                                % ----------------- View Temporal Basis
	                                xlims = [1 56];
	                                ylims = [-0.8 0.8];

	                                if show_subspace == 1
	                                    Kshow = 1:4;
	                                    figure();
	                                    plot(imag(datU(:,Kshow)), 'linewidth', 3);
	                                    xlim(xlims);
	                                    ylim(ylims);
	                                    title({'Imaginary Temporal Subspace'}, 'FontSize', 16)
	                                    legend();
	                                    % faxis;
	                                    print('UScalipr_Subspace_Im', '-dpng', '-r600');

	                                    figure();
	                                    plot(real(datU(:,Kshow)), 'linewidth', 3);
	                                    xlim(xlims);
	                                    ylim(ylims);
	                                    title({'Real Temporal Subspace'}, 'FontSize', 16)
	                                    legend();
	                                    print('UScalipr_Subspace_Re', '-dpng', '-r600');
	                                    % faxis;

	                                    Kshow = 5:8;
	                                    figure();
	                                    plot(real(datU(:,Kshow)), 'linewidth', 3);
	                                    xlim(xlims);
	                                    ylim(ylims);
	                                    title({'Real Temporal Subspace'}, 'FontSize', 16)
	                                    legend();
	                                    print('UScalipr_Subspace_Re5to8', '-dpng', '-r600');
	                                    % faxis;

	                %                     figure();
	                %                     plot(datU(:,Kshow), 'linewidth', 3);
	                %                     xlim(xlims);
	                %                     ylim(ylims);
	                %                     ftitle({'Entire Temporal Subspace'}, 16)
	                %                     legend();
	                %                     print('UScalipr_Subspace_Mag', '-dpng', '-r600');
	                %                     % faxis;
	                                end
	                                % ----------------- View Signal Evolutions
	                                if show_sig_evol == 1
	                                    figure();
	                                    plot(abs(imag(sig_evol_disp)), 'linewidth', 3);
	                                    xlim(xlims);
	                                    title({'Imaginary Sig Evol'}, 'FontSize', 16)
	                                    print('UScalipr_Signal_Im', '-dpng', '-r600');
	                                    % faxis;

	                                    figure();
	                                    plot(abs(real(sig_evol_disp)), 'linewidth', 3);
	                                    xlim(xlims);
	                                    title({'Real Sig Evol'}, 'FontSize', 16)
	                                    print('UScalipr_Signal_Re', '-dpng', '-r600');
	                                    % faxis;

	                                    figure();
	                                    plot(abs(sig_evol_disp), 'linewidth', 3);
	                                    xlim(xlims);
	                                    title({'Complex Mag. Sig Evol '}, 'FontSize', 16)
	                                    print('UScalipr_Signal_Mag', '-dpng', '-r600');
	                                    % faxis;
	                                end

	            %                     U = abs(U);
	                %                 U = real(U);
	                                Phi = permute(datU(:,1:K,:), [ 4 5 6 7 8 1 2 3 ]); % size(Phi)
	                %                     Phi = permute(U(:,1:K), [ 3 4 5 6 7 8 2 9 10 11 1 ]); %size(Phi)
	                                Phi_sq = squeeze(Phi);
	                                writecfl('Phi_save',Phi);
	                                tmp = readcfl('Phi_save');
	            %                     size(tmp)


	                                % ------------------ START NOT RUNNING ------------------
	%                                 disp('Size of input')
	%                                 size(USdat)
	%                                 size(FSsensmaps)
	%                                 size(Phi)

	%                                 % ------ OVERWRITE WITH PICS COMMAND, FOR CHANGING CS RECON COMMAND
	%                                 out_pics_mas_8ch = evalc( 'pics_mas_8ch_noSqueeze = bart(''' + string(cmd_pics) + string(' -B Phi_save ') + ''',USdat, FSsensmaps);' )
	%             %                     disp('Size of output'); % size(pics_mas_8ch_noSqueeze)
	            
	                                pics_mas_8ch_noSqueeze = single(permute( zeros(size(USdat,14),size(USdat,2),size(USdat,3),1,size(FSsensmaps,5),1,K), [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]));
	                    
	                                % --- BATCH MODE
	                                disp('BATCH UScalipr generation')

	                                for batch_id = 0:nbatches-1
	%                                     disp('BATCH:')
	%                                     batch_id
	                                    
	            %                         batch_USdat = USdat(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1);
	            %                         batch_FSsensmaps = FSsensmaps(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1);
	            %                         size(batch_USdat)
	            %                         size(batch_FSsensmaps)


	            % %                         size(pics_mas_8ch_noSqueeze(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1))
	            %                         out_pics_mas_8ch = evalc( 'pics_mas_8ch_noSqueeze(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1) = bart(''' + string(cmd_pics) + ''',batch_USdat, batch_FSsensmaps);' )
	            %                         out_pics_mas_8ch = evalc( 'batch_pics_mas_8ch_noSqueeze = bart(''' + string(cmd_pics) + ''',batch_USdat, batch_FSsensmaps);' )
	            %                         pics_mas_8ch_noSqueeze(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1) = batch_pics_mas_8ch_noSqueeze;
	            %                         clear batch_pics_mas_8ch_noSqueeze

	            %                         out_pics_mas_8ch = evalc( 'pics_mas_8ch_noSqueeze(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1) = bart(''' + string(cmd_pics) + ''',USdat(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1), FSsensmaps(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1));' )
	                                    % QUIET AND SINGLE VERSION
	                                    out_pics_mas_8ch = evalc( 'pics_mas_8ch_noSqueeze(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1) = single(bart(''' + string(cmd_pics) + string(' -B Phi_save ') + ''',USdat(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1), FSsensmaps(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1)));' );

	            %                         out_pics_mas_8ch = evalc( 'tmp = bart(''' + string(cmd_pics) + ''',USdat(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1), FSsensmaps(:,:,:,:,:,:,:,:,:,:,:,:,:,1+(batch_id*batchsize):1+(batch_id*batchsize)+batchsize-1));' )
	            %                         size(tmp)
	                                end
	            

	                                % put slices back in READOUT (1) from SLICE dimension (14)
	                                pics_mas_8ch_noSqueeze = permute(pics_mas_8ch_noSqueeze, [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]); % size(pics_mas_8ch_noSqueeze)

	                                % ------  For multiple ESPIRiT maps, need to combine or throw out others
	                                if (espirit_use_sets == 1) && (espirit_combine_sets == 1)
	                                    pics_mas_8ch_noSqueeze = bart('rss 16', pics_mas_8ch_noSqueeze); % combine multiple sets along ESPIRiT dimension
	                                else
	                                    pics_mas_8ch_noSqueeze = pics_mas_8ch_noSqueeze(:,:,:,:,1,:,:,:,:); % just take first (or only) set
	                                end

	                                % Squeeze back to convenient size
	                                pics_mas_8ch = squeeze(pics_mas_8ch_noSqueeze); % size(pics_mas_8ch)      

	                                % Project temporal basis to image space
	                                UScalipr = Phi_sq * reshape(pics_mas_8ch, size(pics_mas_8ch,1)*size(pics_mas_8ch,2)*size(pics_mas_8ch,3),size(Phi_sq,2)).';
	                                UScalipr = reshape(UScalipr.', size(pics_mas_8ch,1), size(pics_mas_8ch,2), size(pics_mas_8ch,3), size(Phi_sq,1));
	                                
	                                clear pics_mas_8ch_noSqueeze % save pics_mas_8ch to make coeff images
	                                
	                                % UScalipr = abs(UScalipr);
	                                save('UScalipr.mat', 'UScalipr', '-v7.3', '-nocompression');
	                                save('UScalipr_Phi.mat', 'Phi', '-v7.3', '-nocompression');

	                                % delete temporary cfl/hdr version
	                                delete Phi_save.cfl Phi_save.hdr

	                                toc
	                                disp('Done CALIPR Recon')
	                                % ------------------ END NOT RUNNING ------------------
	                                
	                                if save_calipr_nii == 1
	                                    disp('Generating and Saving Par/Rec Files')
	                                    
	                                    % ----- IF NOT RUNNING
	                                    load('UScalipr.mat');

	                                    % ----------------- CALIPR
	                                    clear tmpMR
	                                    tmpMR = MRecon(dat_labfile, dat_rawfile);
	                                    tmpMR.Parameter.Parameter2Read.typ = [ 1 ]; % read in standard data (1) phase correction data (3) noise data (5)
	                                    tmpMR.Parameter.Parameter2Read.Update;

	                                    tmpMR.Data = bart('flip 6', abs(UScalipr)); % 
	                                    tmpMR.Data = permute( tmpMR.Data, [ 1 2 3 5 6 7 4 ]);

	                                    tmpMR.Parameter.ReconFlags.isread = 1;
	                                    tmpMR.Parameter.ReconFlags.issorted = 1;
	                                    tmpMR.Parameter.ReconFlags.ispartialfourier = [0 0];
	                                    tmpMR.Parameter.ReconFlags.isgridded = 0;
	                                    tmpMR.Parameter.ReconFlags.isimspace = [1 1 1];
	                                    tmpMR.Parameter.ReconFlags.iscombined = 1;
	                                    tmpMR.Parameter.ReconFlags.isoversampled = [0 0 1];
	                                    tmpMR.Parameter.ReconFlags.isreadparameter = 1;
	                                    tmpMR.Parameter.ReconFlags.israndphasecorr = 0;
	                                    tmpMR.Parameter.ReconFlags.ispdacorr = 1;
	                                    tmpMR.Parameter.ReconFlags.isdcoffsetcorr = 1;
	                                    tmpMR.Parameter.ReconFlags.isdepicorr = 0;
	                                    tmpMR.Parameter.ReconFlags.ismeasphasecorr = 0;
	                                    tmpMR.Parameter.ReconFlags.isnonlincorr = 0;
	                                    tmpMR.Parameter.ReconFlags.isunfolded = 0;
	                                    tmpMR.Parameter.ReconFlags.iszerofilled = [0 0];
	                                    tmpMR.Parameter.ReconFlags.isrotated = 0;
	                                    tmpMR.Parameter.ReconFlags.isconcomcorrected = 0;
	                                    tmpMR.Parameter.ReconFlags.isgeocorrected = 0;
	                                    tmpMR.Parameter.ReconFlags.issegmentsdivided = 0;
	                                    tmpMR.Parameter.ReconFlags.isecc = 0;
	                                    tmpMR.Parameter.ReconFlags.isaveraged = 1;

	                                    tmpMR.GridderNormalization;
	                                    tmpMR.SENSEUnfold;
	                                    tmpMR.PartialFourier;
	                                    tmpMR.ConcomitantFieldCorrection;
	                                    tmpMR.DivideFlowSegments;
	                                    tmpMR.Average;
	                                    tmpMR.GeometryCorrection;
	                                    tmpMR.RemoveOversampling;
	                                    tmpMR.ZeroFill;
	                                    tmpMR.RotateImage;

	                                    tmpMR.WritePar( [ 'CALIPR.par' ] );
	                                    tmpMR.WriteRec( [ 'CALIPR.rec' ] );
	                                    !dcm2niix -b n -m y -s y -p y -v n -z y -w 1 -f CALIPR CALIPR.par >/dev/null
	                                    !fslmerge -t CALIPR.nii.gz $( for i in {1..56}; do printf "CALIPR_e${i}.nii.gz "; done )
	                                    !mv CALIPR_e1.nii.gz CALIPR_E1.nii.gz
	                                    !rm CALIPR_e*.nii.gz
	                                    
	                                    system(sprintf("antsBrainExtraction.sh -d 3 -k 0 -z 0 -a CALIPR_E1.nii.gz -e %s/T_template0.nii.gz -m %s/T_template0_BrainCerebellumProbabilityMask.nii.gz -f %s/T_template0_BrainCerebellumRegistrationMask.nii.gz -o CALIPR_E1 &>/dev/null",oasis_template_path,oasis_template_path,oasis_template_path))
	                                    !gunzip CALIPR.nii.gz
	                                    !gunzip CALIPR_E1BrainExtractionMask.nii.gz
	                                    clear tmpMR
	                                end
	                                

	                                
	                                evalc('USfft_img = (scale_MWF(2)/imwindow(2))*bart(''normalize 7'', squeeze(abs(USfft(view_slice,:,:,1))));');
	                                evalc('UScalipr_img1 = (scale_MWF(2)/imwindow(2))*bart(''normalize 7'', squeeze(abs(UScalipr(view_slice,:,:,1))));');
	                                evalc('UScalipr_img4 = (scale_MWF(2)/imwindow(2))*bart(''normalize 7'', squeeze(abs(UScalipr(view_slice,:,:,4))));');
	                                evalc('UScalipr_img10 = (scale_MWF(2)/imwindow(2))*bart(''normalize 7'', squeeze(abs(UScalipr(view_slice,:,:,10))));');
	                                evalc('UScalipr_img28 = (scale_MWF(2)/imwindow(2))*bart(''normalize 7'', squeeze(abs(UScalipr(view_slice,:,:,28))));');
	                                UScalipr_phase = (squeeze(angle(UScalipr(view_slice,:,:,1)))).*((scale_MWF(2)-scale_MWF(1))/(2*pi)) + ((scale_MWF(2)-scale_MWF(1))/2);

	                                % view
	                                figure(), imshow([ USfft_img UScalipr_phase UScalipr_img1 UScalipr_img4 UScalipr_img10 UScalipr_img28 ], scale_MWF,'InitialMagnification','fit');
	                                title({'USfft | UScalipr Phase | UScalipr E1/E4/E10/E28 | UScalipr MWF | UScalipr ALPHA ', 'Acceleration R = ' + string(Rfactor) + '    K: ' + string(K) + '  MaskID: ' + num2str(mas_id), ' PICS: ' + string(cmd_pics) })
	                                print('UScalipr_Echoes', '-dpng', '-r600');

	                                if show_coeff == 1
	                                    evalc('coeffimg = squeeze(bart(''normalize 7'', abs(pics_mas_8ch(view_slice,:,:,:))));');
	                                    figure()
	                                    for coefid = 1:size(coeffimg,3)
	                                        if coefid < 16
	                                            subplot(3,5,coefid), imshow([ coeffimg(:,:,coefid) ], imwindow,'InitialMagnification','fit');
	                                            title({' UScalipr C' + string(coefid) })
	                                        end
	                                    end
	                                    print(sprintf('UScalipr_Coeffs.png'), '-dpng', '-r900');
	                                end
	                            end
	                        end
	                        % ---------------------------------------------- END CALIPR RECON ----------------------------------------------
	                        
	                        % ---------------------------------------------- START MWI ANALYSIS ----------------------------------------------
	                        if analyseMWI == 1
	                            disp('Beginning MWI Analysis ')
	                            cd( out_dir_appended + string(cmd_mas) );
	        %                     cmd_pics = calipr_cmd_pics;
	                            cmd_pics = cell_cmd_pics{:};
	                            cmd_pics_clean = strrep(cmd_pics,' ','_')

	                            % to loop over bases
	                            for bid = basis_ids % 1:7
	                                cd( out_dir_appended + string(cmd_mas) );
	                                % directory for K value
	                                RECONdirstr = string(string(cmd_pics_clean) + 'K_' + string(K)) + '_CALIPR_B' + string(bid);
	                                mkdir(RECONdirstr); cd(RECONdirstr);
	%         %                         % directory for K value
	%         %                         RECONdirstr = string(string(cmd_pics_clean) + 'K_' + string(K)) + '_CALIPR';
	%         %                         mkdir(RECONdirstr); cd(RECONdirstr);
	%                                 datU = datU_ALL(:,:,bid);
	                                
	                                %----- TEMPORARY TO SAVE MEMORY: Clear largest files not needed (for not looping multiple recons) to make space for DKspat
	                                % USdat largest in memory but needed, temporarily write to file
	                                USdat = single(USdat);
	                                save('tmp_USdat.mat', 'USdat', '-v7.3', '-nocompression');
	                                clear USdat Eall Eksp Dat UScs iflat iflat_ALL sig_evol sig_evol_ALL sensmaps tmp_slice

	%                                 % ----- Older DK_NNLS analysis for matlab files
	% %                                 UScalipr = UScalipr .* tmpmask;
	%                                 delete(gcp('nocreate')) % resolves occasional issues with matlab parpool
	%                                 DK_NNLS(abs(UScalipr(:,:,:,:)), tmpmask, Resol, TE, TR, T2Range, spwin, nT2, nCores, nThreads, FAE_in_fin_step, RegType, nSpatIter, Size_DSW, Overlap_DSW, SpatRegOpt, SpatRegSetAlpha, SpatReg_Init_FA, SpatReg_Init_T2)

	                                % ----- Updated DK_nii_NNLS analysis
	                                delete(gcp('nocreate')) % resolves occasional issues with matlab parpool
	                                DK_nii_NNLS('CALIPR.nii', 'CALIPR_E1BrainExtractionMask.nii', Resol, TE, TR, T2Range, spwin, nT2, nCores, nThreads, FAE_in_fin_step, RegType, nSpatIter, Size_DSW, Overlap_DSW, SpatRegOpt, SpatRegSetAlpha, SpatReg_Init_FA, SpatReg_Init_T2)

	                                %--- Add USdat back into memory, post-MWI analysis
	                                load('tmp_USdat');
	                                USdat = double(USdat);
	                                delete('tmp_USdat.mat')

	                                if strcmp(RegType, 'Temporal-Only') % 'Temporal-Only'
	                                    nSpatIter=0;
	                                    mwioutfile = dir(sprintf('Temporal_*.mat')); mwioutfile = mwioutfile.name; load(mwioutfile );
	                                    UScalipr_IET2=zeros(size(FAErrorMap_3D));
	%                                     UScalipr_GMT2=zeros(size(FAErrorMap_3D));
	%                                     UScalipr_SNR=zeros(size(FAErrorMap_3D));
	%                                     UScalipr_FNR=zeros(size(FAErrorMap_3D));

	                                    mp=T2_Scale>=mpwin(1) & T2_Scale<=mpwin(2);
	                                    UScalipr_MWF = sum(T2_dist_4D_Temporal(:,:,:,find(T2_Scale >= spwin(1) &T2_Scale <= spwin(2))),4)./(sum(T2_dist_4D_Temporal,4)+eps);
	                                    for row=1:size(FAErrorMap_3D,1)
	                                        for col=1:size(FAErrorMap_3D,2)
	                                            for slice=1:size(FAErrorMap_3D,3)
	                                                dist=squeeze(T2_dist_4D_Temporal(row,col,slice,:));
	                                                UScalipr_IET2(row,col,slice)=exp(dot(dist(mp),log(T2_Scale(mp)))./(sum(dist(mp)+eps)));
	%                                                 UScalipr_GMT2(row,col,slice)=exp(dot(dist,log(T2_Scale))./(sum(dist)+eps));
	                                            end
	                                        end
	                                    end    
	                                else % catch for 'Both' or 'Spatial-Only'
	                                    mwioutfile = dir(sprintf('Spatial_*_Iter_%u.mat',nSpatIter)); mwioutfile = mwioutfile.name; load(mwioutfile );
	                                    UScalipr_IET2=zeros(size(FAErrorMap_3D));
	%                                     UScalipr_GMT2=zeros(size(FAErrorMap_3D));
	%                                     UScalipr_SNR=zeros(size(FAErrorMap_3D));
	%                                     UScalipr_FNR=zeros(size(FAErrorMap_3D));

	                                    mp=T2_Scale>=mpwin(1) & T2_Scale<=mpwin(2);
	                                    UScalipr_MWF = sum(T2_dist_Spa_4D(:,:,:,find(T2_Scale >= spwin(1) &T2_Scale <= spwin(2))),4)./(sum(T2_dist_Spa_4D,4)+eps);
	                                    for row=1:size(FAErrorMap_3D,1)
	                                        for col=1:size(FAErrorMap_3D,2)
	                                            for slice=1:size(FAErrorMap_3D,3)
	                                                dist=squeeze(T2_dist_Spa_4D(row,col,slice,:));
	                                                UScalipr_IET2(row,col,slice)=exp(dot(dist(mp),log(T2_Scale(mp)))./sum(dist(mp)));
	%                                                 UScalipr_GMT2(row,col,slice)=exp(dot(dist,log(T2_Scale))./sum(dist));
	                                            end
	                                        end
	                                    end
	                                end
	                                UScalipr_ALPHA = (180-(180*FAErrorMap_3D));

	%                                 % save structures with MWI metrics of interest
	% %                                 save('UScalipr_MWI.mat', 'UScalipr_MWF', 'UScalipr_ALPHA', 'UScalipr_IET2', 'UScalipr_GMT2',  '-v7.3', '-nocompression');
	%                                 save('UScalipr_MWI.mat', 'UScalipr_MWF', 'UScalipr_ALPHA',  'UScalipr_IET2', '-v7.3', '-nocompression');

	                                if save_mwi_nii == 1
	                                    disp('Generating and Saving MWI Par/Rec Files')

	                                    copyfile('CALIPR.nii', [ string(dat_name) + '.nii' ]);

	                                    nii = load_untouch_nii('CALIPR_E1BrainExtractionMask.nii');

	                                    % BART flip fixes LR flip when loading into matlab
	%                                     nii.img = bart('flip ', UScalipr_MWF);
	                                    nii.img = UScalipr_MWF;
	                                    save_untouch_nii(nii, [ string(dat_name) + '_DK' + string(nSpatIter) + '_MWF.nii' ]);
	%                                     nii.img = bart('flip 1', UScalipr_ALPHA);
	                                    nii.img = UScalipr_ALPHA;
	                                    save_untouch_nii(nii, [ string(dat_name) + '_DK' + string(nSpatIter) + '_ALPHA.nii' ]);
	%                                     nii.img = bart('flip 1', UScalipr_IET2);
	                                    nii.img = UScalipr_IET2;
	                                    save_untouch_nii(nii, [ string(dat_name) + '_DK' + string(nSpatIter) + '_IET2.nii' ]);
	%                                     nii.img = bart('flip 1', UScalipr_GMT2);
	%                                     save_untouch_nii(nii, [ string(dat_name) + '_DK' + string(nSpatIter) + '_IET2.nii' ]);
	%                                     nii.img = bart('flip 1', UScalipr_SNR);
	%                                     save_untouch_nii(nii, [ string(dat_name) + '_DK' + string(nSpatIter) + '_SNR.nii' ]);
	%                                     nii.img = bart('flip 1', UScalipr_FNR);
	%                                     save_untouch_nii(nii, [ string(dat_name) + '_DK' + string(nSpatIter) + '_FNR.nii' ]);
	                                end
	                            end
	                        end
	                        % ---------------------------------------------- END MWI ANALYSIS ----------------------------------------------
	                    end
	                    % ---------------------------------------------- END KVAL LOOP ----------------------------------------------
	                end
	                % ---------------------------------------------- END cmdPICS LOOP ----------------------------------------------
	                % swap back in case not reloading data
	%                 USdat = permute(USdat, [ 14 2 3 4 5 6 7 8 9 10 11 12 13 1 ]); % size(USdat)
	            end
	            % ---------------------------------------------- END Mid LOOP ----------------------------------------------
	        end
	        % ---------------------------------------------- END Rid LOOP ----------------------------------------------
	%         % write info to txt file
	%         outfileID = fopen(out_textfile, 'a+');
	%         fprintf(outfileID, out_info);
	%         fclose(outfileID);

	%         % If txt output is desired:
	%         if print_output == 1 disp(sprintf(out_info)), end
	    end
	end

	toc
	disp('Done all')

end
