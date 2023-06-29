function DK_NNLS(T2Data, matmask, Resol, TE, TR, T2Range, spwin, nT2, nCores, nThreads, FAE_in_fin_step, RegType, nSpatIter, Size_DSW, Overlap_DSW, SpatRegOpt, SpatRegSetAlpha, SpatReg_Init_FA, SpatReg_Init_T2)

    % clear DataFolder DataFileName WOR MaskFileName Spatial_ProcessMode RegTag1 MatPriorityTag Capsule1 muS_Opt_Tag
     %% Data input
    % Input data tag
    InputDataMode               = 'Nifti-4D-Kumar-Format';       %'Raw Dicom Folder', 'Nifti-4D-Kumar-Format', '.mat-4D-Kumar-Format'};
    ImageTypeTag               = 'Philips_Magn: ORIGINAL\PRIMARY\M_SE\M\SE';        % 'Philips_Magn: ORIGINAL\PRIMARY\M_SE\M\SE', 'Siemens_Magn: ORIGINAL\PRIMARY\M\ND', 'Philips_Magn_2Slabs: ORIGINAL\PRIMARY\M_IR\M\IR'    
    
    %% Closing matlabpool if that is open 
    %delete(gcp('nocreate'));
    CloseMatPool();
    

    % time
    starttime = datetime('now');
    disp([ 'Start:                 ' + string(starttime) ]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data prep: Take input T2 data and mask and save in .nii format
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Foldername = pwd;
    nEchoes = size(T2Data,4);
    save_prefix = 'tmp_DK_NNLSdat';
    esp_ms = TE*1e+3;
    TR_ms = TR*1e+3;
    % TE=esp_ms:esp_ms:(esp_ms*nEchoes);

    nii = make_nii(T2Data); % size(nii.img)
    nii_mask = make_nii(matmask); % size(nii_mask.img)

    % Use NiiStr to save echo spacing, TR values
    NiiStr = ['Echo Spacing in ms = ',num2str(esp_ms),'//;' 'Total Echoes = ',num2str(nEchoes), '//;', 'TR in ms= ',num2str(TR_ms), '//'];       


    % For 4D data
    nii.hdr.dime.bitpix     = 32; 
    nii.hdr.dime.datatype   = 16;        
    nii.hdr.dime.dim(1)     = 4;
    nii.hdr.dime.dim(2)     = size(T2Data, 1);
    nii.hdr.dime.dim(3)     = size(T2Data, 2);
    nii.hdr.dime.dim(4)     = size(T2Data, 3);
    nii.hdr.dime.dim(5)     = size(T2Data, 4);
    nii.hdr.hist.descrip = NiiStr;  % Extra string containing echo info
    % For 3D data
    nii_mask.hdr.dime.bitpix     = 32; 
    nii_mask.hdr.dime.datatype   = 16;        
    nii_mask.hdr.dime.dim(1)     = 4;
    nii_mask.hdr.dime.dim(2)     = size(matmask, 1);
    nii_mask.hdr.dime.dim(3)     = size(matmask, 2);
    nii_mask.hdr.dime.dim(4)     = size(matmask, 3);    
    nii_mask.hdr.dime.dim(5)     = size(matmask, 4);   
    nii_mask.hdr.hist.descrip = NiiStr;  % Extra string containing echo info

    save_nii(nii, [ Foldername '/' save_prefix '.nii' ]);
    save_nii(nii_mask, [ Foldername '/' save_prefix '_mask.nii' ]);

        % Clear nii to save memory
    clear nii nii_mask T2Data matmask

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data folder: You can feed mutiple data sets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DataFiles_wt_Dir_Cell = {'/home/bizon/Research/MWI_Development/Data/CALIPER_Dev_FB1/DICOM_SE56_FS/DK_MWI_Test10Slices_lrgDSW2/SE56_FS.nii'}; % OLD
    DataFiles_wt_Dir_Cell = {[ Foldername '/' save_prefix '.nii' ]};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  BRAIN MASK
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % It's advisable to use the masked brain data for processing. Else, noise only part of brain would become computationally expensive 
    UseMask                 = 'Use Mask - *.NII file'                      % 'Use Mask - automatic', 'Use Mask - *.NII file', 'Use Mask - *.ROI file', 'No Mask';
    % Mask could be either automatically constructed, if you have supplied
    % dicoms or you can pass *nii/*.roi mask after creating it.
    % The mask filename must be '*mask.nii*';        '*mask.roi' 
    % MaskExt     = '*mask.nii';                                                      % could be .nii or nii.gz
    MaskExt     = [ save_prefix '_mask.nii' ];                                % could be .nii or nii.gz
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Used if dicom files are inputs: MAY NOT BE WORKING AS OF NOW; NII input recommended.
    % if dicom data is passed, then automatic mask is prepared using BET command of 'fsl'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    EchoNo_4_Mask = 1;  % 1st echo is used to make nifti mask using home-made code, similar to BET mask of FSL
    DicomFileExt    = '';       %'*.dcm', '*.ima', '*.MR'; 
    
    %% Adding appropriate directory etc.
   	[CodeFolder_Parent, CodeFilename, ~] = fileparts(mfilename('fullpath'));
    addpath(genpath(CodeFolder_Parent));   % Add all subfolders
   	JavaPathStr1 = Convert_2_Forwardslash([CodeFolder_Parent, '\Subcodes\ParforProgMonv2\java']); 
    Capsule1.JavaPathStr1 = JavaPathStr1;

    %% T2-scale and FAE scale (in degree)
    % Capsule1.FAE_in = 0;
    % Capsule1.FAE_Fin  = 50;
    % Capsule1.FAE_step = 1;
    Capsule1.FAE_in = FAE_in_fin_step(1);
    Capsule1.FAE_Fin  = FAE_in_fin_step(2);
    Capsule1.FAE_step = FAE_in_fin_step(3);
    


    % Capsule1.Myelin_Cutoff      = [4; 20]*1e-3;   % in sec
    Capsule1.Myelin_Cutoff      = spwin;   % in sec   
    % Capsule1.T2_Initial              = 4;   % in ms   
    % Capsule1.T2_Final               = 2000;   % in ms   
    % Capsule1.nT2                       = 40;   % in ms
    Capsule1.T2_Initial              = T2Range(1)*1e+3;   % in ms   
    Capsule1.T2_Final               = T2Range(2)*1e+3;   % in ms   
    Capsule1.nT2                       = nT2;   % in ms
    
    %% DSW (Data selecting window size)
    %assigning the values to proper fields
    %NOTE: 
    % 16x16x20 previously
    % use larger when limited to 10slice mask
    %       tried 16^3:         maybe ~37 minutes?
    %       tried 22^3:         39 minutes
    %       tried 40x10x40:     60 minutes
    %       tried 10x10x10:     30 minutes!
    % SizeX1_DSW = 10;
    % SizeY1_DSW = 10;
    % SizeZ1_DSW = 10;

    % WOR               = 4;
    WOR               = Overlap_DSW;

    % Size_DSW = [SizeX1_DSW, SizeY1_DSW, SizeZ1_DSW];

    % Resolution of data
    % ResX                    = 1.5;
    % ResY                    = 1.5;
    % ResZ                    = 1.5;
    % Resol                   = [ResX, ResY, ResZ]; 

    % %%
    % Capsule1.nMatlabPoolWorkers_Conv = 60; %ncVirtual;      % One can use virtual cores, which is usually twice the number of physical cores. 
    % Capsule1.nMatlabPoolWorkers_Spat = 30; %ncPhysical;      %  Use only physical core. Using the number equals to the virtual cores does not result in any speed up
    Capsule1.nMatlabPoolWorkers_Conv = nThreads; %ncVirtual;      % One can use virtual cores, which is usually twice the number of physical cores. 
    Capsule1.nMatlabPoolWorkers_Spat = nCores; %ncPhysical;      %  Use only physical core. Using the number equals to the virtual cores does not result in any speed up
    
    Excit_FA    = 90;       % in degree
    Refoc_FA  = 180;   % in degree
    Capsule1.alphaEx_NominalRad      = (double(Excit_FA)/90)*(pi/2);
    Capsule1.alphaRef_NominalRad     = (double(Refoc_FA)/180)*pi;

    %%
   % Previously, Spatial_ProcessMode has onlu '3D' and '2D' option.
   % '2D_3D_Mixed' was later on added
    Spatial_ProcessMode         =  '3D';                                               % '3D', '2D_3D_Mixed', '2D';    
    % Capsule1.No_of_Iterations                = 6;
    Capsule1.No_of_Iterations                = nSpatIter;
    Capsule1.No_of_MixedIterations      = 0;
    Capsule1.MixedIter.WinSz2D            = [20, 20];                          % If chosen mixed iteration, choose this window size
    % Regularization tag
    % RegTag                                           = 'Both';                                   % 'Both', 'Temporal-Only', 'Spatial-Only'
    RegTag                                           = RegType;                                   % 'Both', 'Temporal-Only', 'Spatial-Only'
    
    % If all echoes should be used for the data processing.
    % Used 'AllEchoes' options, other two options may not be working
    Capsule1.SelectEchoTag                = 'AllEchoes';                       %  'AllEchoes', 'EvenEchoes', 'AllEchoes_Except_1st'};    
    %
    % Spatial regularization option
    % Capsule1.OptimizatioTag    = 'Optimize muS';                              % 'Optimize muS', 'Use set value'
    % Capsule1.alpha1                  = 3000;                                              % Only used if  Capsule1.OptimizatioTag = 'Use set value'
    Capsule1.OptimizatioTag    = SpatRegOpt;                              % 'Optimize muS', 'Use set value'
    Capsule1.alpha1                  = SpatRegSetAlpha;                                              % Only used if  Capsule1.OptimizatioTag = 'Use set value'
    %
    % For Capsule1.FAECorrection_Tag, use 'WITH FAE correction'. The other option may not be working.
    Capsule1.FAECorrection_Tag      = 'WITH FAE correction';       %'WITH FAE correction', 'WITHOUT any FAE correction'
    % Refinement level for muS and muFAE
    Capsule1.muS_OptRefine_Tag           = 'Sufficient';                     % 'Sufficient', 'Refined'  ; recommended choice 'Sufficient' 
    Capsule1.muFAI_OptRefine_Tag         = 'Sufficient';                    % 'Sufficient', 'Refined'  ; recommended choice 'Sufficient' 

    %%     %UseMask has to be converted into 'Yes' and 'No' to be consistent with nonGUI code
    if strcmp(UseMask, 'Use Mask - automatic')
        UseMask     = 'Yes';
        MaskType    = 'nii';                % 'nii',  'roi', 'none'  
        MaskFileNamePrefix = '';
        MaskExt     = '*mask.nii*';             % could be .nii or nii.gz    
    elseif strcmp(UseMask, 'Use Mask - *.ROI file')
        UseMask = 'Yes';
        MaskType = 'roi';                % 'nii',  'roi', 'none'    
        MaskFileNamePrefix = '';
        MaskExt     = '*mask.roi';
    elseif strcmp(UseMask, 'Use Mask - *.NII file')
        UseMask     = 'Yes';
        MaskType    = 'nii';                % 'nii',  'roi', 'none'  
        MaskFileNamePrefix = '';
        MaskExt     = '*mask.nii*';             % could be .nii or nii.gz
    else
        UseMask     = 'No';
        MaskType    = 'none';                % 'nii',  'roi', 'none' 
    end
    
    %   ImageType_Desired should be changed based on scanner and sequence type  
    if strcmp(ImageTypeTag, 'Philips_Magn: ORIGINAL\PRIMARY\M_SE\M\SE')
        ImageType_Desired = 'ORIGINAL\PRIMARY\M_SE\M\SE';       % Single slab expt:  'ORIGINAL\PRIMARY\M_SE\M\SE', 'ORIGINAL\PRIMARY\PHASE MAP\P\SE' for philips data
    elseif strcmp(ImageTypeTag, 'Siemens_Magn: ORIGINAL\PRIMARY\M\ND')
        ImageType_Desired = 'ORIGINAL\PRIMARY\M\ND';            % Single slab expt:  'ORIGINAL\PRIMARY\M\ND' for siemens data
    elseif strcmp(ImageTypeTag, 'Philips_Magn_2Slabs: ORIGINAL\PRIMARY\M_IR\M\IR')
        ImageType_Desired = 'ORIGINAL\PRIMARY\M_IR\M\IR';       % DOUBLE slab expt:  'ORIGINAL\PRIMARY\M_IR\M\IR', 'ORIGINAL\PRIMARY\PHASE MAP\P\IR' for philips data   
    else
        % This ImageType_Desired NOT possible as chosen from dropdown menu.
    end
    %
    Capsule1.SeqInfo.b = [esp_ms, nEchoes, TR_ms];
    % Str = ['Echo Spacing in ms = ',num2str(esp_ms),'//;' 'Total Echoes = ',num2str(nEchoes), '//', 'TR in ms= ',num2str(TR), '//'];
    %%
    for i = 1:length(DataFiles_wt_Dir_Cell)
        Filename_wt_Folder = squeeze(DataFiles_wt_Dir_Cell{i});
        
        if(~isempty(Filename_wt_Folder))
            [DataFolder, DataFileName, Ext] = fileparts(Filename_wt_Folder);
            DataFileName    = [DataFileName, Ext];
            MaskFolder      = DataFolder; 
            AnalysisFolder  = DataFolder;            

            % I am assuming that any file with '.roi' is the MaskFileName
            if strcmp(UseMask, 'No')    %string UseMask has already been overwritten
                MaskFileName = '';
            else
                MaskFileName = findFile_with_roi_OR_nii_Ext(MaskFolder, MaskFileNamePrefix, MaskExt);
                %if no file with .roi extension was found
                if isempty(MaskFileName)
                    display('NO file with .roi extension was found');
                    display('No mask will be used');
                    UseMask         = 'No Mask';
                    MaskFileName    = '';
                end
            end

            %  Prefix of the file based on SelectEchoTag
            if strcmp(Capsule1.SelectEchoTag, 'EvenEchoes')
                Temporal_Prefix = 'Temporal_EvenEchoes_';   
                Spatial_Prefix = 'Spatial_EvenEchoes_'; 
            elseif strcmp(Capsule1.SelectEchoTag, 'AllEchoes_Except_1st') 
                Temporal_Prefix = 'Temporal_Exc_1stEchoes_';   
                Spatial_Prefix = 'Spatial_Exc_1stEchoes_';                
            else  % strcmp(Capsule1.SelectEchoTag, 'AllEchoes') 
                Temporal_Prefix = 'Temporal_AllEchoes_';   
                Spatial_Prefix = 'Spatial_AllEchoes_'; 
            end

            AnalysisFolder = DataFolder;        
            MaskFolder = DataFolder;
            RegularizationTag = RegTag; 
            
            MatPriorityTag = 'normal';  % 'normal', 'low', 'below normal', 'above normal', 'high priority', 'real time'
            
            % Some variables in GUI fields are found to be cleared when
            % multiple data are processed; so, reading everything back:
            % GetAllPar_in_ProcessingParameterPanel(handles);

            % EvenEchoes tag not allowed for stimulated correction.
            if strcmp(Capsule1.SelectEchoTag, 'EvenEchoes')
                error('Error by DK. \n EvenEchoes tag not allowed for stimulated correction. Please correct your selection and then run')
            end

            % Open file for debug logging
            fileID = fopen(fullfile(DataFolder, 'Report_Debug_DK.txt'),'w');
            d = debug_point(fileID);
            Capsule1.fileID = fileID;

            % All parameters related to matpool workers
            MatpoolPropConv.nPoolsNeeded  	= Capsule1.nMatlabPoolWorkers_Conv;
            MatpoolPropConv.MatPriorityTag  = MatPriorityTag;
            %
            MatpoolPropSpat.nPoolsNeeded  	= Capsule1.nMatlabPoolWorkers_Spat;
            MatpoolPropSpat.MatPriorityTag  = MatPriorityTag; 

            %%
            % Adam: you may want to change these two vector based on your prior calculated data
            % Initial values for determining muS (spatial regularizatio const for flip angle inhomogeneities)
            % Capsule1.alphaS_DeltaVec = double(logspace(1,5,8));
            Capsule1.alphaS_DeltaVec = SpatReg_Init_FA;                                                                                  
            % Initial values for determining muS (spatial regularizatio const for T2 dist. map) 
            % Capsule1.alphas_Info.alpha1_AllRange     = [10, 100, 200, 350, logspace(log10(5e2), log10(3.5e3), 10), 5e3, 1e4];    % [100, 200, 500, 1000, 2000, 4000, 7500, 10000];
            Capsule1.alphas_Info.alpha1_AllRange     = SpatReg_Init_T2;
            %%
           % If both temporal and spatial regularization ask for the
           % same number of Matpool worker, then open matpool only once
           if (Capsule1.nMatlabPoolWorkers_Conv == Capsule1.nMatlabPoolWorkers_Spat)
                % Open Matlab pool and set priority
                MatpoolPropSpat = OpenMatPool_N_SetPriority(MatpoolPropSpat);  % Conventional regularization is not that computational demanding; so, speed up is possble using virtual cores.  
           end

           Capsule1.NoRegulMode = 'Yes';               % 'Yes', 'No'
           % Performing temporal regularization ONLY
           debug_point(fileID);
           if (strcmp(RegTag, 'Both') || strcmp(RegTag, 'Temporal-Only'))  %'Temporal', 'Spatial', 'Both'
                % Open Matlab pool and set priority
                MatpoolPropConv = OpenMatPool_N_SetPriority(MatpoolPropConv);  % Conventional regularization is not that computational demanding; so, speed up is possble using virtual cores.  

                RegTag1 = 'Temporal-Only';  % Temporarily rewriting RegTag1
                Capsule1.alphas_Info.alpha1 = 1e-6;    % This is a dummy value, not used for temporal regularizatio. But, the code needs it.
                muS_Opt_Tag = 'No'; 	  % 'Yes', 'No';   This is a dummy value, not used for temporal regularizatio. But, the code needs it.

                % [ResVec, NormSpaVec] = QT2RProcess_wt_ParforMonitor_Stim3D(DataFolder, DataFileName, Size_DSW, WOR, Resol, MaskFileName, Spatial_ProcessMode, RegTag1, MatPriorityTag, Capsule1, muS_Opt_Tag);
                % quiet version
                evalc('[ResVec, NormSpaVec] = QT2RProcess_wt_ParforMonitor_Stim3D(DataFolder, DataFileName, Size_DSW, WOR, Resol, MaskFileName, Spatial_ProcessMode, RegTag1, MatPriorityTag, Capsule1, muS_Opt_Tag);');
           end
           %
           debug_point(fileID);
           save TempAnalysis_DK;   %% Saving temporarily 
           load TempAnalysis_DK; 
           %
            % Performing spatial, along with optimization of regularization parameters
            if (~strcmp(RegTag, 'Temporal-Only'))
               % If both temporal and spatial regularization ask for the different number of Matpool workera, then open matpool separately
               if (Capsule1.nMatlabPoolWorkers_Conv ~= Capsule1.nMatlabPoolWorkers_Spat)
                    % Closing the matlabpool first 
                    CloseMatPool(); 
                    % Open Matlab pool and set priority
                    MatpoolPropSpat = OpenMatPool_N_SetPriority(MatpoolPropSpat);  % Conventional regularization is not that computational demanding; so, speed up is possble using virtual cores.  
               end
               %
               debug_point(fileID);
                RegTag1 = 'Spatial-Only'; % Temporarily rewriting RegTag1
                JavaPathStr1         = Capsule1.JavaPathStr1;
                %
                %% First calculate the residual and prior term for alpha_s = 0
                fprintf('Calculation started: for normalizing factor for spatial smoothing on T2 distr \n')
                Capsule1_temp = Capsule1;
                Capsule1_temp.No_of_Iterations   = 1;
                Capsule1_temp.NoRegulMode = 'Yes';               % 'Yes', 'No'

                Capsule1.JavaPathStr1 = JavaPathStr1;                    
                %
                % muS_Opt_Tag = 'Yes' ensures that there are no optimization over muS_Delta_Opt
                muS_Opt_Tag = 'Yes'; 	% 'Yes', 'No'
               %
                T2_Initial	= Capsule1.T2_Initial;      % in ms   
                T2_Final   	= Capsule1.T2_Final;        % in ms   
                nT2        	= Capsule1.nT2;     
                T2_Scale    = logspace(log10(T2_Initial*1e-3), log10(T2_Final*1e-3), nT2);  

                debug_point(fileID);
               % Load residuals for unregularized case
                if strcmpi(Capsule1.FAECorrection_Tag, 'WITH FAE correction')
                    Temporal_Prefix	= ['Temporal_', 'TP', 'Style_', num2str(length(T2_Scale)), 'pts_'];  % saved as .mat files;
                else
                    Temporal_Prefix	= ['Temporal_noStimCorr_',num2str(length(T2_Scale)), 'pts_'];  % saved as .mat files;
                end                   
               %
               Temp_Filename = [Temporal_Prefix, '3D','.mat']; % DK removing temporarily in 2020
               %
               S = load([Convert_2_Forwardslash(DataFolder), Temp_Filename]);                  
               Res_NoReg = S.Res_NoReg;
               Res_Tempo = S.Res_Tempo;
               Capsule1_temp.alphas_Info.Res_0      = Res_NoReg;
               Capsule1_temp.alphas_Info.Res_Tempo  = Res_Tempo;
                %
                fprintf('Calculation ended  : for normalizing factor for spatial smoothing on T2 distr \n');

                %
                Capsule1.NoRegulMode = 'No';               % 'Yes', 'No'
                debug_point(fileID);

                if strcmpi(Capsule1.OptimizatioTag,'Optimize muS')
                    fprintf('Calculation started: normalizing factor for muS optimization \n')
                    muS_OptRefine_Tag = Capsule1.muS_OptRefine_Tag;
                    %Find optimum alphaS = AlphaS_OPt
                    % muS_Opt_Tag = 'Yes' ensures that there are no optimization over muS_Delta_Opt
                    muS_Opt_Tag = 'Yes'; 	% 'Yes', 'No'
                    %
                    if (strcmpi(Spatial_ProcessMode, '2D_3D_Mixed'))
                        Spatial_ProcessModeTemp = '3D';    % For FindOptAlphaS, we will process using 3D, even when '2D_3D_Mixed'
                    else
                        Spatial_ProcessModeTemp = Spatial_ProcessMode;
                    end
                    %
                    debug_point(fileID);
                    Capsule1_temp.alphaS_DeltaVec = Capsule1.alphaS_DeltaVec;
                    Capsule1_temp.alphas_Info.alpha1 = Capsule1.alphas_Info.alpha1_AllRange;

                    % AlphaS_OPt                  = FindOptAlphaS(DataFolder, DataFileName, Size_DSW, WOR, Resol, MaskFileName, Spatial_ProcessModeTemp, RegTag1, MatPriorityTag, Capsule1_temp, muS_Opt_Tag, muS_OptRefine_Tag, fileID);
                    % quiet version
                    evalc('AlphaS_OPt            = FindOptAlphaS(DataFolder, DataFileName, Size_DSW, WOR, Resol, MaskFileName, Spatial_ProcessModeTemp, RegTag1, MatPriorityTag, Capsule1_temp, muS_Opt_Tag, muS_OptRefine_Tag, fileID);');
                    Capsule1.alphas_Info.alpha1 = AlphaS_OPt;
                    fprintf('Calculation ended  : normalizing factor for muS optimization \n')

                    debug_point(fileID);
                else % Pick up the supplied value of muS
                    Capsule1.alphas_Info.alpha1 = Capsule1.alpha1;
                end

                save TempAnalysis_DK2;   %% Saving temporarily 
                load TempAnalysis_DK2;
                %%
                %
                % Now process entire data set
                % muS_Opt_Tag = 'Yes' ensures that there are no optimization over muS_Delta_Opt
                muS_Opt_Tag = 'No'; 	% 'Yes', 'No';
                Capsule1.OptimizatioTag    = 'Use set value';                              % Optimized muS has already been calculated
                %Capsule1.alphas_Info.alpha1 = 1000;
                d = debug_point(fileID);
                fprintf(fileID,'Optimization of muS parameters done;Still has to perform: muS_Delta optimization+ spatially regularized solution. \n');

                % [ResVec, NormSpaVec] = QT2RProcess_wt_ParforMonitor_Stim3D(DataFolder, DataFileName, Size_DSW, WOR, Resol, MaskFileName, Spatial_ProcessMode, RegTag1, MatPriorityTag, Capsule1, muS_Opt_Tag);
                % quiet version
                evalc('[ResVec, NormSpaVec] = QT2RProcess_wt_ParforMonitor_Stim3D(DataFolder, DataFileName, Size_DSW, WOR, Resol, MaskFileName, Spatial_ProcessMode, RegTag1, MatPriorityTag, Capsule1, muS_Opt_Tag);');
                debug_point(fileID);
            end
            %end
            d = debug_point(fileID);
            fprintf(fileID,'Done completely with file %s, function="%s", line %i\n', d.file, d.name, d.line);
            fclose(fileID);
        end
        
        % Clearing all variables and loading variables that was defined outsie loop
        %clearvars('-except',initialVars_OutsideLoopDKU{:}); 
    end
    
    %% Closing matlabpool if that is open 
    CloseMatPool(); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  SAVE DESIRED OUTPUT MAPS FROM FINAL ITERATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mwioutfile = dir('FS_MWI/Spatial_*_Iter_1.mat'); mwioutfile = mwioutfile.name;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  CLEAN UP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmpfiles = dir('TempAnalysis*.mat');
    tmpfiles = {tmpfiles.name};

    for file = tmpfiles
        delete(string(file));
    end

    % delete temporary .nii files used for analysis
    delete([ Foldername '/' save_prefix '.nii' ]);
    delete([ Foldername '/' save_prefix '_mask.nii' ]);

    %Priority set back to normal (Second way of priority setting.)
    if (~isunix)
        try % Even if it does not run, no big deal
            SetMatlabPriorityWindowDKU('normal'); 
        end
    end
    %clearvars;  % Clear all variables active in workspace

    % finish timing
    finishtime = datetime('now');
    elapsed = diff([ starttime finishtime ]);
    fprintf('\n \n');
    disp([ 'Start:                 ' + string(starttime) ]);
    disp([ 'Finish:                ' + string(finishtime) ]);
    disp([ 'Total elapsed time:    ' + string(elapsed) ]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %% Paper to Cite   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HelperMessage1 = 'To ackonwledge the author of this code, please cite this paper:';    
    % HelperMessage2 = '======================================================================';  
    % Message1 = 'Kumar, D., Hariharan, H., Faizy, T. D., Borchert, P., Siemonsen, S., Fiehler, J., et al. (2018).  ';
    % Message2 = 'Using 3D spatial correlations to improve the noise robustness of multi component analysis of 3D multi echo quantitative T2 relaxometry data. Neuroimage 178, 583–601. doi: 10.1016/j.neuroimage.2018.05.026';    
    % HelperMessage = sprintf([HelperMessage1, '\n' , HelperMessage2, '\n',  Message1, '\n' , Message2]);
    % mb = msgbox(HelperMessage, 'Paper to cite');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%
function MaskFileName = findFile_with_roi_OR_nii_Ext(MaskFolder,  MaskFileNamePrefix, MaskExt)
    %ensure that MaskFolder is ended with filesep
    MaskFolder = Convert_2_Forwardslash(MaskFolder);
    fs = sprintf('%s%s', MaskFolder, MaskFileNamePrefix, MaskExt);
    MaskInfo = dir(fs);
    % if both *.nii and .nii.gz exists, pickup the first one
    MaskFileName = MaskInfo(1).name;  
end 

% function Folder = Check_4_FileSep(Folder)
%     if ~strcmp(Folder(end), filesep)
%         Folder = [Folder, filesep];
%     end
% end

function subFolders = FindSubfolder(ParentFolder)
    % Get a list of all files and folders in this folder.
    files = dir(ParentFolder);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
end