%
% JULY 23- AUG 13 2015:
%
%  1) Two types of masks are included: (a) BrainMask (2) TissueClassMask
%
%  2) Rather than making map for all data @once, it's being done in 2 steps
%      (a) TissueClassMask .*Data       ==> solve 
%      (b)(1-TissueClassMask) .*Data    ==> solve
%       Results from step (a) and (b) can be combined now
%
%   3) Following function is most recent and robust
%    Ds = getSpatialDerivOperator3D_wt_Masking(......)
%
%    It also accounts for tissue class
%
%    4)Fixed confusing notation between X-axis, Y-axis, Vertical Axis,optio
%   Horizontal Axis.   
%   !Vertical axis is sometimes called x
%
%    5) Calculate weightages of residuals and priors
%
%    6) handles flip angle error (FAE) correction as well.
%
%   Change 1:
%           Filter was unbalanced as the central voxel did not have weight
%           of 0.5
%
%   Change 3:
%           Wrong looping was resulting in tiling. Trying to fix
%
%  Change 4: Implemented 3D filter;
%             Neighbours of central voxel in 3D filter have same weights
%
%
%  Change 5:
%             Neighbours of central voxel have weight inversely propotion
%             to distance.
%           Later on, need to fix this tag:
%           Weight_Tag = 'Inv-distance';   % 'Equal', 'Inv-distance'    % TEMPORARY . NOT ACTIVE: WOULD THINK LTER ON THIS
%
%  Change 6: New functions introduced 
%           It turned out that muT_map was converted from 3D-array to
%           1D-Vec in wrong manner. FIXED
%
%  Change 7: 
%           Some exisiting data-processing parameters have been saved, only 
%           for the purpose of reporting in paper. Since these parameters are not needed otherwise, a "keyword" is specially created: named- muS_Opt_Tag  
%
%  Change 8: 
%   Jan reported some bug with inverse distance filter; fixed it.
%
%  Change 9: 
%   Current version, uses lsqnonlin for FAE-refinement; this is done on
%   slice by slice basis. Spatial smoothness is implemented even for FAEmap
%   Older function: CalculateDeltaErrorMap_Voxelwise removed
%
%  Change 10:
%   Matlab priorty set on window using system command.
%
%  Change 11:
%       CLEARING UNNECESSARY VARIABLES and ASSIGNING HIGH PRIORITY
%       : I clear unnecessary variables after keeping necessary ones.
%       If I do not clear the variables, sometimes Matpool worked quits for 
%       reasons such as memory limitation and others unspecified reason. Error not reproducible.
%       For instance, if I increase "number of T2" to 200, these types of
%       error occurs for certain at different stage of execution.
%       CLEARING unnecessary variables and ASSIGNING high priority to Matlab helps.
%
%       REMOVED: Inverse distance filter is removed.
%
%       NEW DEFAULT: For conventional regularization, algorithm similar to
%       Prasloski's has been implemented. This is also the default for
%       CONVENTIONAL regularization.
%
%        MORE SPEEDUP for CONVENTIONAL:
%        In conventional, using all virtual cores results in speed up. This is no while we only use physical core for spatial regularization
%
%Version 5 vs Version 6:
%       I am rewriting "getSpatialDerivOperator3D_wt_Masking" ==> This
%       allows me to use different filtering for T2dist-smoothing, FAE-smoothing
%
%        Version 5 also tried to use lbfgsb for FAE, which were commented out. Now removed in Version 6.   
%
%   TO BE CHANGED LATER:
%       Resol is hardwired at many places.  ==> Check Resol is a row vec
%       FAE correction still uses 2D smoothing, change to 3D smoothing ==>
%       may have to rewrite Invert_4_FAEMAP_FULL
%
function [ResVec, NormSpaVec] = QT2RProcess_wt_ParforMonitor_Stim3D(DataDir, DataFile, Size_DSW, WOR, Resol, MaskFileName, Spatial_ProcessMode, RegTag, MatPriorityTag, Capsule1, muS_Opt_Tag)
    Debug = 0;  % 0 for no, 1 for yes
    
    % Monitor process
% %     if ismac % Code to run on Mac plaform
% %         MonitorProcess = 'Yes_WO_JAVA'; % 'Yes', 'No', 'Yes_WO_JAVA'
% %     elseif isunix % Code to run on Linux plaform
% %         MonitorProcess = 'Yes';       	% 'Yes', 'No', 'Yes_WO_JAVA'
% %     elseif ispc   % Code to run on Windows platform
% %         MonitorProcess = 'Yes';       	% 'Yes', 'No', 'Yes_WO_JAVA'
% %     else
% %         disp('Platform not supported')
% %     end
    %
    MonitorProcess = 'Yes_WO_JAVA';
    % 
    if (isempty(MaskFileName)== 1)
        MaskInfo.UseMask        = 'No'; 
    else
        MaskInfo.UseMask        = 'Yes';
    end
    MaskInfo.MaskFileName   = MaskFileName;
    
    %
    RegRefineTag.muS_OptRefine_Tag = Capsule1.muS_OptRefine_Tag; 
    RegRefineTag.muFAI_OptRefine_Tag = Capsule1.muFAI_OptRefine_Tag;
    %
% % %     TempoalAlgoTag   	= Capsule1.TempoalAlgoTag;
    nPoolsNeeded_Conv	= Capsule1.nMatlabPoolWorkers_Conv;
    nPoolsNeeded_Spat	= Capsule1.nMatlabPoolWorkers_Spat;    
    IterInfo.No_of_Iterations    = Capsule1.No_of_Iterations;
    IterInfo.No_of_Mixed2DIterations = Capsule1.No_of_MixedIterations;
    %alpha1Vec           = Capsule1.alpha1;
    %alpha1Vec           = Capsule1.alphas_Info.alpha1 .* Capsule1.alphas_Info.Factor_s;   
    alpha1Vec           = Capsule1.alphas_Info.alpha1;
        
    JavaPathStr         = Capsule1.JavaPathStr1;
    SelectEchoTag       = Capsule1.SelectEchoTag;
    FAECorrection_Tag   = Capsule1.FAECorrection_Tag; 
    fileID              = Capsule1.fileID;
    
    Myelin_Cutoff    = Capsule1.Myelin_Cutoff;
    %
    EPG_Param.alphaEx_NominalRad    = Capsule1.alphaEx_NominalRad;
    EPG_Param.alphaRef_NominalRad   = Capsule1.alphaRef_NominalRad;
    EPG_Param.T1                    = 1000e-3;  % Constant T1 
    %
    ParamMatpool.nPoolsNeeded_Conv	= nPoolsNeeded_Conv;
    ParamMatpool.nPoolsNeeded_Spat 	= nPoolsNeeded_Spat;
    ParamMatpool.MatPriorityTag     = MatPriorityTag;
    %
    T2_Initial                      = Capsule1.T2_Initial;   % in ms   
    T2_Final                        = Capsule1.T2_Final;   % in ms   
    nT2                             = Capsule1.nT2;    
    
    % Processing condition
    ProcessCond.muS_Opt_Tag = muS_Opt_Tag; 	% 'Yes', 'No';
    ProcessCond.B1_min      = 0.01*Capsule1.FAE_in;
    ProcessCond.B1_max      = 0.01*Capsule1.FAE_Fin;
    ProcessCond.StepSize_LR	= 0.01*Capsule1.FAE_step;
    ProcessCond.muT                 = double(logspace(-5, -3, 50));  % 100 points rather than 140 should be fine.
    ProcessCond.BaselineTag         = 'No';                 	%'Yes', 'No'     
    ProcessCond.UseSigmoidSpatial   = 'No';                     %'Yes', 'No'
    ProcessCond.UseExpVarWeightsTemporal   = 'No';              %'Yes', 'No'  
    ProcessCond.T2_Scale            = logspace(log10(T2_Initial*1e-3), log10(T2_Final*1e-3), nT2);  
    %ProcessCond.T2_Scale            = [linspace(5e-3, 40e-3, 25), linspace(41e-3, 130e-3, 50), linspace(135e-3, 2000e-3, 75)];

    ProcessCond.nFAE_HR             = 101;                      % Points on high resolution FAEVec scale.
    ProcessCond.nTE_FAE             = 32;                       % 16, 32    
    ProcessCond.muS_Delta           = 0.01;                     % Regularization constants  
    ProcessCond.SelectEchoTag   	= SelectEchoTag;    
    ProcessCond.FAECorrection_Tag	= FAECorrection_Tag; 
    ProcessCond.TemporalRegStyle    = 'TP';                     % 'DK' --> Kumar's JMIR, 'TP' --> Thomas Prasloski's stype
    %
    %ProcessCond.alphaS_DeltaVec =  [10, double(logspace(2,3,5)), double(logspace(3.3,5,6))];
    ProcessCond.alphaS_DeltaVec =  Capsule1.alphaS_DeltaVec;
    %
    IterInfo.MixedIter_WinSz2D = Capsule1.MixedIter.WinSz2D;
    %
    if ((~strcmpi(RegTag, 'Temporal-Only')) && (~strcmpi(Spatial_ProcessMode, '3D'))) 
        Filter.FilterDim  = '2D';
        Filter.FilterSize =  [3, 3];     
    else 
        Filter.FilterDim  = '3D';
        %Filter.FilterSize =  [3, 3, 3]; 
        Filter.FilterSize =  [3, 3, 3]; 
    end
    
    Filter.FilterTag  = 'IDW-NonLoc' ;   % 'IDW-NonLoc', 'ExpDecay-NonLoc', 'Equal-Local' 
    %
    %ConditionStr  = [num2str(length(ProcessCond.T2_Scale)), 'pts_', '_UseExpVarWeightsTemporal_', ProcessCond.UseExpVarWeightsTemporal, '_FilterTag_', Filter.FilterTag, '_'];
    ConditionStr  = [num2str(length(ProcessCond.T2_Scale)), 'pts_',  'FilterTag_', Filter.FilterTag, '_'];
    %
    if strcmpi(FAECorrection_Tag, 'WITH FAE correction')
        Temporal_Prefix	= ['Temporal_', ProcessCond.TemporalRegStyle, 'Style_', num2str(length(ProcessCond.T2_Scale)), 'pts_'];  % saved as .mat files;
        if strcmpi(Spatial_ProcessMode, '3D')
            Spatial_Prefix  = ['Spatial_UsingFAI6_', ConditionStr, num2str(Filter.FilterSize(1)), 'By', num2str(Filter.FilterSize(2)), 'By', num2str(Filter.FilterSize(3))]; 
        else
            Spatial_Prefix  = ['Spatial_UsingFAI6_', ConditionStr, num2str(Filter.FilterSize(1)), 'By', num2str(Filter.FilterSize(2))]; 
        end
    else
        Temporal_Prefix	= ['Temporal_noStimCorr_',num2str(length(ProcessCond.T2_Scale)), 'pts_'];  % saved as .mat files;
        if strcmpi(Spatial_ProcessMode, '3D')
            Spatial_Prefix  = ['Spatial_noStimCorr_', ConditionStr, num2str(Filter.FilterSize(1)), 'By', num2str(Filter.FilterSize(2)), 'By', num2str(Filter.FilterSize(3))];  
        else
            Spatial_Prefix  = ['Spatial_noStimCorr_', ConditionStr, num2str(Filter.FilterSize(1)), 'By', num2str(Filter.FilterSize(2))]; 
        end
    end
    
% %     for ia = 1:Len_a
% %         alpha1_0 = alpha1Vec(ia);
% %         [ResVec(1,ia), NormSpaVec(1,ia)] = Process_EntireData(DataFile, DataDir, JavaPathStr, ParamMatpool, alpha1_0, Temporal_Prefix, Spatial_Prefix, Filter, Resol, Spatial_ProcessMode, MaskInfo, Size_DSW, WOR, RegTag, IterInfo, EPG_Param, Myelin_Cutoff, ProcessCond, Debug, MonitorProcess, RegRefineTag, fileID);  
% %     end
 
    [ResVec, NormSpaVec] = Process_EntireData(DataFile, DataDir, JavaPathStr, ParamMatpool, alpha1Vec, Temporal_Prefix, Spatial_Prefix, Filter, Resol, Spatial_ProcessMode, MaskInfo, Size_DSW, WOR, RegTag, IterInfo, EPG_Param, Myelin_Cutoff, ProcessCond, Debug, MonitorProcess, RegRefineTag, fileID);   
end

function [ResVec, NormSpaVec] = Process_EntireData(DataFileName, DataDir, JavaPathStr, ParamMatpool, alpha1Vec, Temporal_Prefix, Spatial_Prefix, Filter, Resol, Spatial_ProcessMode, MaskInfo, DSW, WOR, RegTag, IterInfo, EPG_Param, Myelin_Cutoff, ProcessCond, Debug, MonitorProcess, RegRefineTag, fileID)    
    [StartTime, ~] = What_s_TheTime();
    display(['Data processingg started @ ', StartTime]);
    TemporalRegStyle = ProcessCond.TemporalRegStyle; % 'DK' --> Kumar's JMIR, 'TP' --> Thomas Prasloski's stype
    
    % Reserving space for these two variables
    RES0 = NaN;         CSP0 = NaN;
    Len_a = length(alpha1Vec);
    ResVec = nan(1,Len_a);
    NormSpaVec = nan(1,Len_a);
    %
    FilterSize = Filter.FilterSize;
    %
    muS_OptRefine_Tag   = RegRefineTag.muS_OptRefine_Tag; 
    muFAI_OptRefine_Tag = RegRefineTag.muFAI_OptRefine_Tag ;
    
    % Parameters for data processing
    B1_min              = ProcessCond.B1_min;
    B1_max              = ProcessCond.B1_max;
    muT                 = ProcessCond.muT;                	% 100 points rather than 140 should be fine.
    T2_Scale            = ProcessCond.T2_Scale;   
    nFAE_HR             = ProcessCond.nFAE_HR;          	% Points on high resolution FAEVec scale.
    muS_Delta           = ProcessCond.muS_Delta;         	% Regularization constants      
    muS_Opt_Tag         = ProcessCond.muS_Opt_Tag; 	% 'Yes', 'No';
    BaselineTag         = ProcessCond.BaselineTag;      	%'Yes', 'No';
    UseSigmoidSpatial   = ProcessCond.UseSigmoidSpatial;   	% 'yes', 'no' 
    nTE_FAE             = ProcessCond.nTE_FAE;              % For FAE optimization, we do not full TE
    UseExpVarWeightsTemporal = ProcessCond.UseExpVarWeightsTemporal; % 'Yes', 'No' initially no
% % %     TempoalAlgoTag     	= ProcessCond.TempoalAlgoTag;   % 'Prasloski_Type_LCurve_Fast', 'DK_LCurve_Slow' 
    alphaS_DeltaVec = ProcessCond.alphaS_DeltaVec;
    SelectEchoTag = ProcessCond.SelectEchoTag;
    FAECorrection_Tag = ProcessCond.FAECorrection_Tag; 
    
    UseMask             = MaskInfo.UseMask;
    MaskFileName        = MaskInfo.MaskFileName;
    alphaEx_NominalRad  = EPG_Param.alphaEx_NominalRad;
    alphaRef_NominalRad = EPG_Param.alphaRef_NominalRad;
    T1 = EPG_Param.T1;  % Constant T1
    
    %Size of data selecting window
    SizeV =  DSW(1)+WOR;  SizeH = DSW(2)+WOR;   SizeS = DSW(3)+WOR;
    
    % Low resolution FAEVec for conventional regularization
    B1_min = ProcessCond.B1_min;
    B1_max = ProcessCond.B1_max;
    StepSize_LR = ProcessCond.StepSize_LR;
    
    % Low resolution FAEVec
    nFAE_LR = (B1_max-B1_min)/StepSize_LR +1;  
    FAEVec_LR = linspace(B1_min,B1_max,nFAE_LR);
    % High resolution FAEVec for refinement
    FAEVec_HR = linspace(B1_min,B1_max,nFAE_HR);

    % Relevant folders
    CodeFolder      = DataDir;    
    AnalysisFolder  = DataDir;
    MaskFolder      = DataDir;
    
    %%
    B1_Info.B1_min = B1_min;
    B1_Info.B1_max = B1_max;
    B1_Info.nFAE_LR = nFAE_LR;
  
    %% Load data
    % Read data from Nifti/ mat data file
    Data_Filename_wt_Folder = [Convert_2_Forwardslash(DataDir), DataFileName];
    DataFile_Ext = Find_the_Ext(Data_Filename_wt_Folder);
    if strcmp(DataFile_Ext, '.nii')        % Nifti data file 
        nii =   load_untouch_nii(Data_Filename_wt_Folder);
        T2Data4D =   double(nii.img); 
        b = Extract_Numbers_in_nii_DescripField(nii.hdr.hist.descrip);        
        esp = b(1);
        TE = (1:1:b(2)).*esp;   % echo time in msec
        TR = b(3);              % repetition time in msec
    else
        S = load(Data_Filename_wt_Folder);
        T2Data4D    =  S.T2Data;
        TE = S.TE;
        %TE = TE';   %Non simulation
        esp = TE(2)-TE(1);
        %esp = S.esp;  
    end
    
    %% TEMPORARILY REDUCING DATA SIZE FOR FAST PROCESSING
    
    
    %%
    % Choose only few/all echoes for processing based on the choice of SelectEchoTag
    if strcmp(SelectEchoTag, 'EvenEchoes')
        TE = TE(2:2:end);
        T2Data4D = T2Data4D(:,:,:,2:2:end);
    elseif strcmp(SelectEchoTag, 'AllEchoes_Except_1st') 
        TE = TE(2:1:end);
        T2Data4D = T2Data4D(:,:,:,2:1:end);               
    else  % strcmp(SelectEchoTag, 'AllEchoes') 
        % Do nothing
    end

    T2Data4D(isnan(T2Data4D))= 0.0;  % REPLACING nan by zeros
    %
    %% Read the ROI/nii mask file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For the time being, BrainMask, TissueClassMask are same. 
    % In general, I intend BrainMask to be BET mask and "TissueClassMask" to
    % classify WM, lesions, GM=CSF (everything else)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp = size(T2Data4D);
    dims = temp(1:3); 
    [~, ~, MaskExt] = fileparts(MaskFileName); 
    
    if strcmp(UseMask, 'Yes')  
        MaskFilename_with_Dir = [Convert_2_Forwardslash(MaskFolder), MaskFileName]; 
        if (strcmp(MaskExt, '.nii') || strcmp(MaskExt, '.gz'))  
            if strcmp(MaskExt, '.gz')
                gunzip(MaskFilename_with_Dir);
                MaskFilename_with_Dir = strrep(MaskFilename_with_Dir, '.gz', '');
            end
            
            nii =   load_untouch_nii(MaskFilename_with_Dir);
            Mask3D = nii.img;
        elseif strcmp(MaskExt, '.roi')
            Mask3D = double(readroi_DK(MaskFilename_with_Dir, dims'));
        else 
            display('Exiting now; Mask only three formats accepted: .nii, .nii.gz, .roi');
        end

        % If input DataFile_Ext is '.mat', then rotate data to have
        % consistent orientation as if data were saved nifti
        if strcmpi(DataFile_Ext, '.mat')
            SzM = size(Mask3D);
            Mask3D_New = zeros(SzM(2), SzM(1), SzM(3));
            T2Data4D_New = zeros([size(T2Data4D)]);
            for l = 1:size(T2Data4D,4)
                for k = 1:size(T2Data4D,3)
                    T2Data4D_New(:,:,k,l) = flipdim(rot90(squeeze(T2Data4D(:,:,k,l)),-1),1);
                end
            end
            T2Data4D = T2Data4D_New;
        end
        
    else
        Mask3D = ones([dims]);
    end 
    Mask3D = double(Mask3D);   % Making sure that Mask3D is double.

    %Combining mask with data
    for k = 1:dims(3) 
        for l = 1:length(TE)
            T2Data4D(:,:,k,l) = squeeze(T2Data4D(:,:,k,l)).*squeeze(Mask3D(:,:,k));
        end
    end
    
    %% Mask out non brain part
    Sz_T2Data       = size(T2Data4D);
    Sz_T2Dist       = [Sz_T2Data(1), Sz_T2Data(2), Sz_T2Data(3), length(T2_Scale)];
    Sz3D0        	= [Sz_T2Data(1), Sz_T2Data(2), Sz_T2Data(3)];
    [Vert1, Vert2]  = Find_ROI_Coord_4_MaskedData(Mask3D);
    T2Data4D        = T2Data4D(Vert1(1):Vert2(1),Vert1(2):Vert2(2),Vert1(3):Vert2(3),:);
    Mask3D          = Mask3D(Vert1(1):Vert2(1),Vert1(2):Vert2(2),Vert1(3):Vert2(3));
%%
    % If data consists of a single slice, then reshape the data to be 4D
    Sz0 = size(T2Data4D);
    if (length(size(T2Data4D))==3)
        T2Data4D = reshape(T2Data4D, [Sz0(1), Sz0(2), 1, Sz0(3)]);   % Single slice
    end  
    %% Relevant data portion to be processed
    TotSlices = size(T2Data4D,3);
    SliceSeries = 1:1:TotSlices;
    nSlices = length(SliceSeries);
    %
    temp = size(T2Data4D);
    dims = temp(1:3);      % Recalculate dims
    % Testing the orientation of mask wrt image using the middle slice
    ClickTag = 'No';  % 'Yes', 'No'
    SliceNo2Check = SliceSeries(round(length(SliceSeries)/2));
    Test_Mask_Orient_WRT_Img(T2Data4D, TE, Mask3D, dims, SliceNo2Check, ClickTag);   

    %%
    % For now, I am making "BrainMask" and "TissueClassMask" identical. In
    % general, "TissueClassMask" can be modified to incorporate tissue edge
    % info.
    BrainMask       = Mask3D;
    TissueClassMask = Mask3D;
    NoOfVoxel       = sum(BrainMask(:));
    %% Ensuring TE is inputed correctly
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (size(TE,1) >1) % Ensure TE to be a row vector
        TE = TE';         
    end 
    LenTE = length(TE); LenT2 = length(T2_Scale);
    if max(TE) >=100 %it means that TE is in msec.
        display('TE was in "msec"; converted to "sec"')
        TE = TE/1000;
        esp = esp/1000;
    end
    %
    LenT2_Offset = (LenT2+1).*strcmp(BaselineTag, 'Yes') + (LenT2)*(~strcmp(BaselineTag, 'Yes'));  % Eqv to concise if statement     
    %
    A = exp(-kron(TE', 1./T2_Scale));
    A_with_Offset = [A ones(size(A,1), 1)];   
    %
    %% Different number of malabpool workers for conventional and spatial
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % In conventional, we can use all virtual cores; while we only use physical core for spatial regularization
    MatpoolPropConv.nPoolsNeeded  	= ParamMatpool.nPoolsNeeded_Conv;
    MatpoolPropConv.MatPriorityTag  = ParamMatpool.MatPriorityTag;
    %
    MatpoolPropSpat.nPoolsNeeded  	= ParamMatpool.nPoolsNeeded_Spat;
    MatpoolPropSpat.MatPriorityTag  = ParamMatpool.MatPriorityTag;    
    %
    %% Simulate System matrix: SysMat_HR and part of Jacobian: JMat_HR
    TE_FAE = [1:1:nTE_FAE]* esp;
    LenTE_FAE = length(TE_FAE); 
    % Calculating system matrix and Jacobian
    SysMat_HR = cell(length(FAEVec_HR),1); % This system matrix would be used for FAE opt.
    JMat_HR = cell(length(SysMat_HR),1);
    if strcmpi(FAECorrection_Tag, 'WITH FAE correction')
        dx_HR = FAEVec_HR(2)-FAEVec_HR(1); 
        % Precalculate system matrxi and Jacobian at predefined grid
        parfor k = 1:length(FAEVec_HR)
            DeltaVox = FAEVec_HR(k);
            SysMat_HR{k,1} = Given_FAE_ij_Find_SysMatrixWithOffset(DeltaVox, alphaEx_NominalRad, alphaRef_NominalRad, TE_FAE, esp, T1, T2_Scale, BaselineTag); 
        end
        % Calculating Jacobian
        parfor k = 2:length(FAEVec_HR)
            JMat_HR{k,1} = (SysMat_HR{k,1}-SysMat_HR{k-1,1})./dx_HR;
        end
        JMat_HR{1,1} = JMat_HR{2,1};
        save([Convert_2_Forwardslash(AnalysisFolder), 'SysMat_n_Others.mat'], 'SysMat_HR','JMat_HR');  
    end
    %
    %% Perform temporal regularization
    Temp_Filename = [Temporal_Prefix, '3D','.mat']; 
    initialVars_b4_TempoRegul = who;                                % Store a list of the names of all the variables currently in the workspace
    if (strcmp(RegTag, 'Both') || strcmp(RegTag, 'Temporal-Only'))  %'Temporal', 'Spatial', 'Both'
        % RUNNING conventional regularization
      	muT_Info.muTMax = 1e-3;
        muT_Info.muTMin = 1e-5; 
        if strcmpi(TemporalRegStyle, 'DK')
            mu     = double(logspace(log10(muT_Info.muTMin), log10(muT_Info.muTMax), 50));            
            [FAErrorMap_3D, T2_dist_4D_Temporal, muT_Map_3D, Res_NoReg, Res_Tempo] = estimate_FAEMap_N_TemporalRegul_DK(T2Data4D, T2_Scale, TE, BrainMask, esp, alphaEx_NominalRad, alphaRef_NominalRad, T1, FAEVec_LR, MatpoolPropConv, JavaPathStr, BaselineTag, mu, UseExpVarWeightsTemporal, MonitorProcess) ; 
        else % TP style
            [FAErrorMap_3D, T2_dist_4D_Temporal, muT_Map_3D, Res_NoReg, Res_Tempo] = estimate_FAEMap_N_TemporalRegul   (T2Data4D, T2_Scale, TE, BrainMask, esp, alphaEx_NominalRad, alphaRef_NominalRad, T1, FAEVec_LR, MatpoolPropConv, JavaPathStr, BaselineTag, muT_Info, FAECorrection_Tag, UseExpVarWeightsTemporal, MonitorProcess); 
        end
        T2_dist_4D_Temporal = AddMarginToProcessedData(T2_dist_4D_Temporal, Sz3D0, Vert1, Vert2);
        muT_Map_3D          = AddMarginToProcessedData(muT_Map_3D, Sz3D0, Vert1, Vert2);
        if strcmpi(FAECorrection_Tag, 'WITH FAE correction')
            FAErrorMap_3D       = AddMarginToProcessedData(FAErrorMap_3D, Sz3D0, Vert1, Vert2);
            save([Convert_2_Forwardslash(AnalysisFolder),Temp_Filename], 'T2_dist_4D_Temporal', 'FAErrorMap_3D', 'muT_Map_3D', 'T2_Scale', 'BaselineTag', 'Res_NoReg', 'Res_Tempo');   
        else
            save([Convert_2_Forwardslash(AnalysisFolder),Temp_Filename], 'T2_dist_4D_Temporal', 'muT_Map_3D', 'T2_Scale', 'BaselineTag', 'Res_NoReg', 'Res_Tempo');  
        end
    end
    clearvars('-except',initialVars_b4_TempoRegul{:})  % Clearing all variable that were generated during conventional regul.
    
    %% Loading Temporal Analysis
    S = load([Convert_2_Forwardslash(AnalysisFolder),Temp_Filename]);
    T2_Scale = S.T2_Scale;
    T2_dist_4D_Temporal = S.T2_dist_4D_Temporal(Vert1(1):Vert2(1),Vert1(2):Vert2(2),Vert1(3):Vert2(3),:); 
    T2_dist_4D_Guess    = T2_dist_4D_Temporal;
    muT_Map_3D          = S.muT_Map_3D(Vert1(1):Vert2(1),Vert1(2):Vert2(2),Vert1(3):Vert2(3)); 
    if strcmpi(FAECorrection_Tag, 'WITH FAE correction')
        FAErrorMap_3D   = S.FAErrorMap_3D(Vert1(1):Vert2(1),Vert1(2):Vert2(2),Vert1(3):Vert2(3));
    else
        FAErrorMap_3D   = zeros([size(muT_Map_3D)]);  % Dummy FAErrorMap_3D with zeros
    end

    %% SMoothing out the first guess
% % % % % % %     FAErrorMap_3D = Interpolate_3D_Based_on_IDW_Fn_Parloop(FAErrorMap_3D, Resol,BrainMask);
% % % %     FAErrorMap_3D  = Filter_KSpace_3D(FAErrorMap_3D);   % Performing k-space filter April 2018: commented out
    
    %% Initial step in the estimation of the spatial regularization constant
    % Also erode mask bit agressively to find out the muS ==> ensures that regions outside brain are not included.
    nLayers = 6;
    Mask3D_Eroded       = ErodeMask3D(Mask3D, nLayers);
    muT_Map_3D_Eroded   = muT_Map_3D .*Mask3D_Eroded;
    Vec                 = muT_Map_3D_Eroded(:);
    muT_Map_3D_ErodedVec = Vec(abs(Vec)>eps); % removing non-zero values
    median_muT_Map      = median(muT_Map_3D_ErodedVec(:));
    
    %% Decide which loop to run
    if strcmpi(Spatial_ProcessMode, '2D')
        response_LoopSBS = true;
        response_Loop4D = false;
        No_of_Iterations_SBS = IterInfo.No_of_Iterations;
    elseif strcmpi(Spatial_ProcessMode, '3D')
        response_LoopSBS = false;
        response_Loop4D = true;  
        No_of_Iterations_4D = IterInfo.No_of_Iterations;
    elseif (strcmpi(Spatial_ProcessMode, '2D_3D_Mixed') && strcmpi(muS_Opt_Tag, 'Yes'))
        response_LoopSBS = false;
        response_Loop4D = true; 
        No_of_Iterations_SBS = 0;
        No_of_Iterations_4D = IterInfo.No_of_Iterations;   
    else % (strcmpi(Spatial_ProcessMode, '2D_3D_Mixed') && ~strcmpi(muS_Opt_Tag, 'Yes'))
        response_LoopSBS = true;
        response_Loop4D = true;
        No_of_Iterations_SBS = IterInfo.No_of_Mixed2DIterations;
        No_of_Iterations_4D = IterInfo.No_of_Iterations-IterInfo.No_of_Mixed2DIterations;
    end
    
 %%   
     muS_Delta_Opt_Tag = 'No'; 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Spatial 2D-SBS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %% Open Matlab pool and set priority for spatial regulairzation
% %     MatpoolPropSpat = OpenMatPool_N_SetPriority(MatpoolPropSpat);  % For regularization, using virtual core does not result in speed up. Usually some worker aborts resulting in error.
    Vars_4_Spat2DRegul = who;   % List all variables available here.
    muS = 0.0;                  % Dummy value
    %% Spatial 2D : slice-by-slice      
    if ((~strcmp(RegTag, 'Temporal-Only')) && (response_LoopSBS == true)) 
        %ProcessingMode = '2D';
        ProcessingMode = Spatial_ProcessMode;
        % Overwriting 
        if strcmpi(Spatial_ProcessMode, '2D_3D_Mixed')
            SizeV = IterInfo.MixedIter_WinSz2D(1);
            SizeH = IterInfo.MixedIter_WinSz2D(2);   
            Filter.FilterDim  = '2D';
            Filter.FilterSize =  [3, 3]; 
        end       
        
        %% Window size for FAI smoothing..
        SizeV1 = 20;  SizeH1 = 20;  SizeS1 = 1;    WOR1 = 4;   % DK think how to pass that from GUI nicely; hardwired values using for now  
        LengthV1 = SizeV1-WOR1; LengthH1 = SizeH1-WOR1;  LengthS1 = SizeS1-0; 
        DSWInfo2D_FAI.LengthV1  = LengthV1;
        DSWInfo2D_FAI.LengthH1  = LengthH1;
        DSWInfo2D_FAI.LengthS1  = LengthS1;
        DSWInfo2D_FAI.WOR1      = WOR1;       
        %% 
        disp('Optimzation for muS_Delta (FAE-mapping) started');
        % muS_Opt_Tag = 'Yes' is set whenever we solve for spatially regularized soln
% % %         if ((~strcmpi(muS_Opt_Tag, 'Yes')) && strcmpi(FAECorrection_Tag, 'WITH FAE correction') && strcmpi(muS_Delta_Opt_Tag, 'No')) 
% % %             % OPtimization will be done over small volume "SV"
% % %             % notice that mus and mus_delta optimization do not run simultaneously
% % %             T2Data4D_SV         = DefineOptimizationVolume(T2Data4D);
% % %             T2_dist_4D_Guess_SV = DefineOptimizationVolume(T2_dist_4D_Guess);
% % %             FAErrorMap_3D_SV    = DefineOptimizationVolume(FAErrorMap_3D);
% % %             BrainMask_SV        = DefineOptimizationVolume(BrainMask);
% % %             % DK: overwrting for '2D_3D_Mixed'                                                                  %% DK TEMPORARY CHANGES??
% % %             if (strcmpi(Spatial_ProcessMode, '2D_3D_Mixed'))
% % %                 Filter.FilterDim = '3D';
% % %                 Filter.FilterSize = [3 3 3];
% % %             end
% % %             % muS_Delta_Opt is already multiplied by an appropriate factor
% % %             [alphaS_Delta_Opt, Factor_Delta] = Optimize_muS_Delta_4_FAESmoothing(FAErrorMap_3D_SV, T2Data4D_SV, T2_dist_4D_Guess_SV, BrainMask_SV, FAEVec_HR, SysMat_HR, JMat_HR, LenTE_FAE, LenT2, alphaS_DeltaVec, Debug, BaselineTag, Resol, Filter, MatpoolPropSpat, JavaPathStr, DSWInfo2D_FAI, AnalysisFolder, MonitorProcess, muFAI_OptRefine_Tag, fileID);
% % %             muS_Delta_Opt                       = alphaS_Delta_Opt*Factor_Delta;  
% % % % % %         else
% % % % % %             Factor_Delta = 3000;
% % % % % %             muS_Delta_Opt = 2500*Factor_Delta;            
% % %         end


%%
        Factor_Delta = 3000;
        muS_Delta_Opt = 2500*Factor_Delta;
%%
        % DK:reverting BACK for '2D_3D_Mixed'
        if (strcmpi(Spatial_ProcessMode, '2D_3D_Mixed'))                                                       %% DK TEMPORARY CHANGES?
            Filter.FilterDim = '2D';
            Filter.FilterSize = [3 3];
        end
        %
        muS_Delta_Opt_Tag = 'Yes';    % We have already performed
        %% Again using smaller part of data fro optimization; not using _SV tag though
        % notice that mus and mus_delta optimization do not run simultaneously
        if strcmpi(muS_Opt_Tag, 'Yes') % 
            T2Data4D            = DefineOptimizationVolume(T2Data4D);
            T2_dist_4D_Guess    = DefineOptimizationVolume(T2_dist_4D_Guess);
            FAErrorMap_3D       = DefineOptimizationVolume(FAErrorMap_3D);
            BrainMask           = DefineOptimizationVolume(BrainMask); 
            %No_of_StepsS        = size(T2Data4D,3);
        end
        %
        Sz3 = size(BrainMask);
        DSWInfo2D_T2D.SizeV = SizeV;
        DSWInfo2D_T2D.SizeH = SizeH;
        DSWInfo2D_T2D.SizeS = SizeS;
        DSWInfo2D_T2D.WOR = WOR;
        %
        disp('Optimzation for muS_Delta (FAE-mapping) DONE');
        %
        %% Saving some necessary variables, before clearing remaining.
        VarsNeeded_2Run_Spat2D = cell(15,1);
        VarsNeeded_2Run_Spat2D{1,1} = 'ProcessingMode'; 
        VarsNeeded_2Run_Spat2D{2,1} = 'p';              
        VarsNeeded_2Run_Spat2D{3,1} = 'q';              VarsNeeded_2Run_Spat2D{4,1} = 'muS'; 
        VarsNeeded_2Run_Spat2D{5,1} = 'RES0';       	VarsNeeded_2Run_Spat2D{6,1}= 'CSP0';
        VarsNeeded_2Run_Spat2D{7,1} ='DSWInfo2D_FAI';      VarsNeeded_2Run_Spat2D{8,1} = 'Filter';       VarsNeeded_2Run_Spat2D{9,1} = 'muS_Delta_Opt'; 
        VarsNeeded_2Run_Spat2D{10,1} = 'FAECorrection_Tag';
        VarsNeeded_2Run_Spat2D{11,1} = 'muS_Delta_Opt_Tag';
        VarsNeeded_2Run_Spat2D{12,1} = 'DSWInfo2D_T2D';
        % VarsNeeded_2Run_Spat2D{12,1} = 'LengthV';        VarsNeeded_2Run_Spat2D{13,1} = 'LengthH';       VarsNeeded_2Run_Spat2D{14,1} = 'LengthS';
        % VarsNeeded_2Run_Spat2D{15,1} = 'No_of_StepsV';   VarsNeeded_2Run_Spat2D{16,1} = 'No_of_StepsH';  VarsNeeded_2Run_Spat2D{17,1} = 'No_of_StepsS';
        
        % Remove empty cells in cell array
        VarsNeeded_2Run_Spat2D = FindEmptyElementsInCell_n_RemoveThose(VarsNeeded_2Run_Spat2D);
        % Concatenate two cell array along columns
        Vars_4_Spat2DRegul = cat(1 , Vars_4_Spat2DRegul, VarsNeeded_2Run_Spat2D);        
        
        %% BACK TO ACTUAL PROCESSING: Arrange to 4D
        T2_dist_Spa_4D 	= zeros(size(T2Data4D,1), size(T2Data4D,2), nSlices, LenT2_Offset);
        PriorSpa_4D    	= zeros(size(T2Data4D,1), size(T2Data4D,2), nSlices, LenT2_Offset);
        y_fit_4D        = zeros(size(T2Data4D,1), size(T2Data4D,2), nSlices, LenTE); 
        for p = 1:length(alpha1Vec) 
            alpha1 = alpha1Vec(p);
            muS = alpha1*median_muT_Map;  
            for q = 1:No_of_Iterations_SBS

                Resol2D = [Resol(1), Resol(2)];
                %
                if (~strcmpi(muS_Opt_Tag, 'Yes') && (q ~= 1))
                % if (~strcmpi(muS_Opt_Tag, 'Yes'))    % TEMPORARILY
                    % DK: overwrting for '2D_3D_Mixed'                                                                
                    if (strcmpi(Spatial_ProcessMode, '2D_3D_Mixed'))
                       Filter.FilterDim = '3D';
                       Filter.FilterSize = [3 3 3];               
                       %%
                       if strcmpi(FAECorrection_Tag, 'WITH FAE correction')
                            SizeV1 = 12; SizeH1 = 12; SizeS1 = 16; WOR1 = 4;
                            LengthV1 = SizeV1-WOR1; LengthH1 = SizeH1-WOR1;  LengthS1 = SizeS1-WOR1; 
                            DSWInfo.LengthV1 = LengthV1;
                            DSWInfo.LengthH1 = LengthH1;
                            DSWInfo.LengthS1 = LengthS1;
                            DSWInfo.WOR1 = WOR1;  
                            disp('Refining FAErrorMap 3D') 
                            [FAErrorMap_3D, ~, ~]   = Invert_4_FAEMAP_3D_FULL(FAErrorMap_3D,T2Data4D, T2_dist_4D_Guess, BrainMask, FAEVec_HR, SysMat_HR, JMat_HR, LenTE_FAE, LenT2, muS_Delta_Opt, Debug, BaselineTag, Resol, Filter, MatpoolPropSpat, JavaPathStr, DSWInfo, MonitorProcess, fileID);                 
                       else
                            SzD3t = size(T2Data4D);
                            if (length(SzD3t) == 4)
                                FAErrorMap_3D  = zeros(SzD3t(1), SzD3t(2), SzD3t(3)); % 3D volume/multislice
                            else
                                FAErrorMap_3D  = zeros(SzD3t(1), SzD3t(2)); % SINGLE SLICE
                            end
                       end                          
                    else
                       disp('Refining FAErrorMap slice by slice')  
                       if strcmpi(FAECorrection_Tag, 'WITH FAE correction')
                            [FAErrorMap_3D, ~, ~] = Invert_4_FAEMAP_SBS_FULL(FAErrorMap_3D, T2Data4D, T2_dist_4D_Guess, BrainMask, FAEVec_HR, SysMat_HR, JMat_HR, LenTE_FAE, LenT2, muS_Delta_Opt, Debug, BaselineTag, Resol2D, Filter, MatpoolPropSpat, JavaPathStr, DSWInfo2D_FAI,  MonitorProcess);
                       else
                            SzD3t = size(T2Data4D);
                            if (length(SzD3t) == 4)
                                FAErrorMap_3D  = zeros(SzD3t(1), SzD3t(2), SzD3t(3)); % 3D volume/multislice
                            else
                                FAErrorMap_3D  = zeros(SzD3t(1), SzD3t(2)); % SINGLE SLICE
                            end
                       end                       
                    end
                    FAErrorMap_3D  = FAErrorMap_3D.*BrainMask;
                    disp('Refining FAErrorMap DONE');
                end
                
                % DK:reverting BACK for '2D_3D_Mixed'
                if (strcmpi(Spatial_ProcessMode, '2D_3D_Mixed'))                                                      
                    Filter.FilterDim = '2D';
                    Filter.FilterSize = [3 3];
                end
        
                %%
                % For debugging
                d = debug_point(fileID);
                fprintf(fileID, 'Starting: iteration: %s out of %s', q, No_of_Iterations_SBS);
                %
                [T2_dist_Spa_4D, TimeTaken_4D_Spatial, ConvScoreInfo, y_fit_4D,PriorSpa_4D] = ...
                        SpatialRegul_3D_SBS(T2Data4D, T2_dist_4D_Guess, FAErrorMap_3D, TissueClassMask, BrainMask, muT_Map_3D, DSWInfo2D_T2D, T2_Scale, muS, UseSigmoidSpatial, BaselineTag, AnalysisFolder, ProcessingMode, JavaPathStr, MatpoolPropSpat, FilterSize, Resol, TE, esp, T1, Debug,UseExpVarWeightsTemporal, muS_Opt_Tag, Filter, FAECorrection_Tag, MonitorProcess);        
                %
                fprintf(fileID, 'Ended: iteration: %s out of %s', q, No_of_Iterations_SBS);

                % Calculate residual and norm of spatial prior for an appropriate flag
                if strcmpi(muS_Opt_Tag, 'Yes') % If this condition is not true, then both output variales: RES0, CSP0 would be NaNs.
                   [RES0, CSP0] = CalculatePrior_N_Res_3Dt(T2Data4D, y_fit_4D, PriorSpa_4D, BrainMask);  
                end

                % Saving analysis files if muS_Opt_Tag is not 'Yes'
                if ~strcmpi(muS_Opt_Tag, 'Yes')
                    if (strcmpi(Spatial_ProcessMode, '2D_3D_Mixed'))
                        Spatial_Filename_2D = [Spatial_Prefix, 'alpha1_', num2str(alpha1), '_SliceBySlice_3Dt', '_MixedIter_', num2str(q),'.mat']; 
                    else
                        Spatial_Filename_2D = [Spatial_Prefix, 'alpha1_', num2str(alpha1), '_SliceBySlice_3Dt', '_Iter_', num2str(q),'.mat']; 
                    end
                else
                    if (strcmpi(Spatial_ProcessMode, '2D_3D_Mixed'))
                        Spatial_Filename_2D = ['RegOpt_', 'alpha1_', num2str(alpha1), '_', Spatial_Prefix, 'alpha1_', num2str(alpha1), '_SliceBySlice_3Dt', '_MixedIter_', num2str(q),'.mat']; 
                    else
                        Spatial_Filename_2D = ['RegOpt_', 'alpha1_', num2str(alpha1), '_',Spatial_Prefix, 'alpha1_', num2str(alpha1), '_SliceBySlice_3Dt', '_Iter_', num2str(q),'.mat']; 
                    end
                    
                end
                    % Without margin values stored
                    T2_dist_Spa_4D0         = T2_dist_Spa_4D;
                    T2_dist_4D_Temporal0    = T2_dist_4D_Temporal;
                    muT_Map_3D0             = muT_Map_3D;
                    FAErrorMap_3D0          = FAErrorMap_3D;
                    % Adding margin before saving
                    T2_dist_Spa_4D      = AddMarginToProcessedData(T2_dist_Spa_4D, Sz3D0, Vert1, Vert2);
                    T2_dist_4D_Temporal = AddMarginToProcessedData(T2_dist_4D_Temporal, Sz3D0, Vert1, Vert2);
                    muT_Map_3D          = AddMarginToProcessedData(muT_Map_3D, Sz3D0, Vert1, Vert2);                    
                    if strcmpi(FAECorrection_Tag, 'WITH FAE correction')
                        FAErrorMap_3D = AddMarginToProcessedData(FAErrorMap_3D, Sz3D0, Vert1, Vert2);
                        save([Convert_2_Forwardslash(AnalysisFolder),Spatial_Filename_2D], 'T2_dist_Spa_4D', 'T2_dist_4D_Temporal', 'muT_Map_3D', 'FAErrorMap_3D','muS', 'muT', 'T2_Scale','BaselineTag');            
                    else
                        save([Convert_2_Forwardslash(AnalysisFolder),Spatial_Filename_2D], 'T2_dist_Spa_4D', 'T2_dist_4D_Temporal', 'muT_Map_3D', 'muS', 'muT', 'T2_Scale','BaselineTag');              
                    end
                    % For debugging
                    d = debug_point(fileID);
                    
                    % Restoring without margin values
                    T2_dist_Spa_4D = T2_dist_Spa_4D0;
                    T2_dist_4D_Temporal = T2_dist_4D_Temporal0;
                    muT_Map_3D = muT_Map_3D0;
                    FAErrorMap_3D = FAErrorMap_3D0;
                    % For debugging
                    d = debug_point(fileID);
                % end    % DK temporarily commenting out
                
                
                
                % For debugging
                d = debug_point(fileID);
                %
                if ((No_of_Iterations_SBS > 1) || (length(alpha1Vec)>1))
                    T2_dist_4D_Guess = T2_dist_Spa_4D;
                end
                
% % %                 clearvars('-except',Vars_4_Spat2DRegul{:})   % Clearing all variable that were generated spatial-3D regul; T2_dist_4D_Guess is also saved.
% % %                 Vars_4_Spat2DRegul = who;  % AGAIN Listing all variables available here.
            end
            % For debugging
            d = debug_point(fileID);
            fprintf(fileID, 'Outside q-loop');
            %
            ResVec(1,p)        = RES0;
            NormSpaVec(1,p)    = CSP0;
            clearvars('-except',Vars_4_Spat2DRegul{:})   % Clearing all variable that were generated spatial-3D regul; T2_dist_4D_Guess is also saved.
            Vars_4_Spat2DRegul = who;  % AGAIN Listing all variables available here.
        end
        % For debugging
        d = debug_point(fileID);
        fprintf(fileID, 'Outside p-loop');
    end   

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Spatial 3D      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %% Spatial 3D: process in 3D mode
    Vars_4_Spat3DRegul = who;  % Listing all variables needed for spat3D.
    if ((~strcmp(RegTag, 'Temporal-Only')) && (response_Loop4D == true))
        display('Computation started in 3D mode');
        % ProcessingMode = '3D';
        ProcessingMode = Spatial_ProcessMode;
        %
        if strcmpi(Spatial_ProcessMode, '2D_3D_Mixed')
            Filter.FilterDim  = '3D';
            Filter.FilterSize =  [3, 3, 3];  
        end
        %% For FAI related processing
        % It will be used twice: i) Optimize_muS_Delta_4_FAESmoothing ii) Invert_4_FAEMAP_3D_FULL
        % Number of steps will be calculated inside the callinf code
        SizeV1_FAI = 12; SizeH1_FAI = 12; SizeS1_FAI = 16; WOR1_FAI = 4;% <-------------------------- DK look into this
        %
        LengthV1_FAI = SizeV1_FAI-WOR1_FAI; LengthH1_FAI = SizeH1_FAI-WOR1_FAI;  LengthS1_FAI = SizeS1_FAI-WOR1_FAI; 
        DSWInfo_FAI.LengthV1    = LengthV1_FAI;
        DSWInfo_FAI.LengthH1    = LengthH1_FAI;
        DSWInfo_FAI.LengthS1    = LengthS1_FAI;
        DSWInfo_FAI.WOR1        = WOR1_FAI;        
        %

        if ((~strcmpi(muS_Opt_Tag, 'Yes')) && strcmpi(FAECorrection_Tag, 'WITH FAE correction') && strcmpi(muS_Delta_Opt_Tag, 'No'))
            FAErrorMap_3D = Filter_KSpace_3D(FAErrorMap_3D.*BrainMask);
            
            %Optimization will be done over small volume "SV"
            %notice that mus and mus_delta optimization do not run simultaneously
            T2Data4D_SV         = DefineOptimizationVolume(T2Data4D);
            T2_dist_4D_Guess_SV = DefineOptimizationVolume(T2_dist_4D_Guess);
            FAErrorMap_3D_SV    = DefineOptimizationVolume(FAErrorMap_3D);
            BrainMask_SV        = DefineOptimizationVolume(BrainMask);

            %muS_Delta_Opt is already multiplied by an appropriate factor
            [alphaS_Delta_Opt, Factor_Delta]    = Optimize_muS_Delta_4_FAESmoothing(FAErrorMap_3D_SV,  T2Data4D_SV, T2_dist_4D_Guess_SV, BrainMask_SV, FAEVec_HR, SysMat_HR, JMat_HR, LenTE_FAE, LenT2, alphaS_DeltaVec, Debug, BaselineTag, Resol, Filter, MatpoolPropSpat, JavaPathStr, DSWInfo_FAI, AnalysisFolder, MonitorProcess, muFAI_OptRefine_Tag, fileID);
            muS_Delta_Opt                       = alphaS_Delta_Opt*Factor_Delta;  
% %         else  
% %             %muS_Delta_Opt =  0.0;
% %             Factor_Delta = 1.8135e+03;
% %             muS_Delta_Opt = 1642*Factor_Delta;    
        end     
        
      
        %% Again using smaller part of data for optimization; not using _SV tag though
        % notice that mus and mus_delta optimization do not run simultaneously
        if strcmpi(muS_Opt_Tag, 'Yes') % 
            T2Data4D            = DefineOptimizationVolume(T2Data4D);
            T2_dist_4D_Guess    = DefineOptimizationVolume(T2_dist_4D_Guess);
            FAErrorMap_3D       = DefineOptimizationVolume(FAErrorMap_3D);
            BrainMask           = DefineOptimizationVolume(BrainMask); 
            % Masked version with margin off
            Sz_T2Data       = size(T2Data4D);
            Sz3D0            = [Sz_T2Data(1), Sz_T2Data(2), Sz_T2Data(3)];
            [Vert1, Vert2]  = Find_ROI_Coord_4_MaskedData(BrainMask);
            T2Data4D        = T2Data4D(Vert1(1):Vert2(1),Vert1(2):Vert2(2),Vert1(3):Vert2(3),:);
            BrainMask          = BrainMask(Vert1(1):Vert2(1),Vert1(2):Vert2(2),Vert1(3):Vert2(3));                 
        end
        Sz3 = size(BrainMask);
        DSWInfo3D_T2D.SizeV = SizeV;
        DSWInfo3D_T2D.SizeH = SizeH;
        DSWInfo3D_T2D.SizeS = SizeS;
        DSWInfo3D_T2D.WOR = WOR;
    
% %         DSWInfo3D_T2D.Vertex1 = [1 1 1];
% %         DSWInfo3D_T2D.Vertex2 = [Sz3(1) Sz3(2) Sz3(3)];       
        
        % Saving some necessary variables, before clearing remaining.
        VarsNeeded_2Run_Spat3D       = cell(15,1);
        VarsNeeded_2Run_Spat3D{1,1}  = 'ProcessingMode'; 
        VarsNeeded_2Run_Spat3D{2,1} = 'p';            
        VarsNeeded_2Run_Spat3D{3,1}  = 'q';              VarsNeeded_2Run_Spat3D{4,1}= 'muS';
        VarsNeeded_2Run_Spat3D{5,1} = 'RES0';       	 VarsNeeded_2Run_Spat3D{6,1}= 'CSP0';
        VarsNeeded_2Run_Spat3D{7,1} = 'DSWInfo_FAI';   	 VarsNeeded_2Run_Spat3D{8,1} = 'Filter';        VarsNeeded_2Run_Spat3D{9,1} = 'muS_Delta_Opt';  
        VarsNeeded_2Run_Spat3D{10,1} = 'FAECorrection_Tag';    
        VarsNeeded_2Run_Spat2D{11,1} = 'muS_Delta_Opt_Tag';
        VarsNeeded_2Run_Spat3D{12,1} = 'DSWInfo3D_T2D'; 
        %VarsNeeded_2Run_Spat3D{12,1}  = 'LengthV';        VarsNeeded_2Run_Spat3D{13,1} = 'LengthH';        VarsNeeded_2Run_Spat3D{14,1} = 'LengthS';
        %VarsNeeded_2Run_Spat3D{15,1}  = 'No_of_StepsV';   VarsNeeded_2Run_Spat3D{16,1} = 'No_of_StepsH';   VarsNeeded_2Run_Spat3D{17,1} = 'No_of_StepsS';
 
        % Remove empty cells in cell array
        VarsNeeded_2Run_Spat3D = FindEmptyElementsInCell_n_RemoveThose(VarsNeeded_2Run_Spat3D);
        % Concatenate two cell array along columns
        Vars_4_Spat3DRegul = cat(1 , Vars_4_Spat3DRegul, VarsNeeded_2Run_Spat3D);

        display('Spatial 3D processing starting ...')
        for p = 1:length(alpha1Vec) 
            alpha1 = alpha1Vec(p);
            muS = alpha1*median_muT_Map;
            for q = 1:No_of_Iterations_4D 
                %initialVars_b4_SpatialRegul_3D = who;      %Store a list of the names of all the variables currently in the workspace.                
                %Refining FAErrorMap_3D
                FAECorrectionRun_Tag = strcmpi(Spatial_ProcessMode, '2D_3D_Mixed') || (q ~= 1) ;
                if (FAECorrectionRun_Tag == 1)                
                    if (~strcmpi(muS_Opt_Tag, 'Yes'))
                        display('Refining FAE-map ...')  
                        if strcmpi(FAECorrection_Tag, 'WITH FAE correction')
                           [FAErrorMap_3D, ~, ~]= Invert_4_FAEMAP_3D_FULL(FAErrorMap_3D,  T2Data4D, T2_dist_4D_Guess, BrainMask, FAEVec_HR, SysMat_HR, JMat_HR, LenTE_FAE, LenT2, muS_Delta_Opt, Debug, BaselineTag, Resol, Filter, MatpoolPropSpat, JavaPathStr, DSWInfo_FAI, MonitorProcess, fileID);
                        else
                            SzD3t = size(T2Data4D);
                            if (length(SzD3t) == 4)
                                FAErrorMap_3D  = zeros(SzD3t(1), SzD3t(2), SzD3t(3)); % 3D volume/multislice
                            else
                                display('Single slice can not be processed in 3D mode');
                                exit
                            end
                        end
                    end
                    FAErrorMap_3D = FAErrorMap_3D.*BrainMask; 
                end
                % For debugging
                d = debug_point(fileID);
                fprintf(fileID, 'Starting: iteration: %s out of %s', q, No_of_Iterations_4D);
                %
                % y_fit_4D, PriorSpa_4D are meaningful only when strcmpi(muS_Opt_Tag, 'Yes') 
                [T2_dist_Spa_4D, TimeTaken_4D_Spatial, ConvScoreInfo, y_fit_4D, PriorSpa_4D] ...
                    = SpatialRegul_3D(T2Data4D, T2_dist_4D_Guess, FAErrorMap_3D, TissueClassMask, BrainMask, muT_Map_3D, DSWInfo3D_T2D, T2_Scale, muS, UseSigmoidSpatial, BaselineTag, AnalysisFolder, ProcessingMode, JavaPathStr, MatpoolPropSpat, FilterSize, Resol, TE, esp, T1, Debug, UseExpVarWeightsTemporal, muS_Opt_Tag, Filter, FAECorrection_Tag, MonitorProcess);
                %
                fprintf(fileID, 'Ended: iteration: %s out of %s', q, No_of_Iterations_4D);
                %
                if strcmpi(muS_Opt_Tag, 'Yes') % If this condition is not true, then both output variables: RES0, CSP0 would be NaNs.
                    RES0 = ConvScoreInfo.ResOverall_SpatReg;
                    CSP0 = ConvScoreInfo.CSP_Overall_SpatReg;  
                    [RES0, CSP0, ResMap, NormPriorMap] = CalculatePrior_N_Res_3Dt(T2Data4D, y_fit_4D, PriorSpa_4D, BrainMask);  
                end
                
                cd(CodeFolder); %Going back to folder containing code helps when running on parallel cluster.
                if ((No_of_Iterations_4D > 1) || (length(alpha1Vec)>1))
                    T2_dist_4D_Guess = T2_dist_Spa_4D;
                end
                %
                % Save analysis files if muS_Opt_Tag is not 'Yes'
                if ~strcmpi(muS_Opt_Tag, 'Yes')
                    if (strcmpi(Spatial_ProcessMode, '2D_3D_Mixed'))
                        Spatial_Filename_3D = [Spatial_Prefix, 'alpha1_', num2str(alpha1), '_3Dt', '_MixedIter_', num2str(No_of_Iterations_SBS+q),'.mat']; 
                    else
                        Spatial_Filename_3D = [Spatial_Prefix, 'alpha1_', num2str(alpha1), '_3Dt', '_Iter_', num2str(q),'.mat']; 
                    end                     
                    %
                    % Without margin values stored
                    T2_dist_Spa_4D0         = T2_dist_Spa_4D;
                    T2_dist_4D_Temporal0    = T2_dist_4D_Temporal;
                    muT_Map_3D0             = muT_Map_3D;
                    FAErrorMap_3D0          = FAErrorMap_3D;
                    % Adding margin before saving
                    T2_dist_Spa_4D  = AddMarginToProcessedData(T2_dist_Spa_4D, Sz3D0, Vert1, Vert2);
                    muT_Map_3D      = AddMarginToProcessedData(muT_Map_3D, Sz3D0, Vert1, Vert2);                    
                    if strcmpi(FAECorrection_Tag, 'WITH FAE correction')
                        FAErrorMap_3D = AddMarginToProcessedData(FAErrorMap_3D, Sz3D0, Vert1, Vert2);
                        save([Convert_2_Forwardslash(AnalysisFolder),Spatial_Filename_3D], 'T2_dist_Spa_4D', 'FAErrorMap_3D', 'muT_Map_3D', 'ConvScoreInfo', 'muS', 'muT', 'T2_Scale','TimeTaken_4D_Spatial', 'BaselineTag'); 
                    else
                       save([Convert_2_Forwardslash(AnalysisFolder),Spatial_Filename_3D], 'T2_dist_Spa_4D', 'muT_Map_3D', 'ConvScoreInfo',  'muS', 'muT', 'T2_Scale','TimeTaken_4D_Spatial', 'BaselineTag');                         
                    end 
                    % Restoring without margin values
                    T2_dist_Spa_4D = T2_dist_Spa_4D0;
                    T2_dist_4D_Temporal = T2_dist_4D_Temporal0;
                    muT_Map_3D = muT_Map_3D0;
                    FAErrorMap_3D = FAErrorMap_3D0;
                end
% % %                 clearvars('-except',Vars_4_Spat3DRegul{:});
% % %                 Vars_4_Spat3DRegul = who;  % AGAIN Listing all variables needed for spat3D.                
            end
            % For debugging
            d = debug_point(fileID);
            fprintf(fileID, 'Outside q-loop');
            %
            ResVec(1,p)        = RES0;
            NormSpaVec(1,p)    = CSP0;
            %
            clearvars('-except',Vars_4_Spat3DRegul{:});
            Vars_4_Spat3DRegul = who;  % AGAIN Listing all variables needed for spat3D.  
        end 
        % For debugging
        d = debug_point(fileID);
        fprintf(fileID, 'Outside p-loop');
    end
    
    %% Saving MWF MAP AS NII
    % Save analysis files if muS_Opt_Tag is not 'Yes'
% %     if ~strcmpi(muS_Opt_Tag, 'Yes')
% %         T2_dist_Spa_4D = AddMarginToProcessedData(T2_dist_Spa_4D, Sz3D0, Vert1, Vert2);
% %         % filename
% %         if (strcmpi(Spatial_ProcessMode, '2D_3D_Mixed'))
% %             Spatial_Filename_nii = ['MWF_Using_MixedIter', '.nii'];   
% %         elseif (strcmpi(Spatial_ProcessMode, '3D'))
% %             Spatial_Filename_nii = ['MWF_Using_3D_OP', '.nii'];   
% %         else % 2D
% %             Spatial_Filename_nii = ['MWF_Using_2D_OP', '.nii'];   
% %         end
% %         %
% %         MWF_Map = sum(T2_dist_Spa_4D(:,:,:,find(T2_Scale >= Myelin_Cutoff(1) &T2_Scale <= Myelin_Cutoff(2))),4)./(sum(T2_dist_Spa_4D,4)+eps);% Recalculate MWF map
% %         % If mask is nifti file, then "nii" already exists here; else create "nii"structure.
% %         if ~exist('nii')   % For simulated data, there is no header information.
% %             try 
% %                 % For exptl data, nifti files are created from dicom data 
% %                 % file; So, loading "nii" structure from there.
% %                 NiiFileName = findFile_with_GivenExt(AnalysisFolder, '', '*.nii');
% %                 nii =   load_untouch_nii([Convert_2_Forwardslash(AnalysisFolder), NiiFileName]);
% %             catch
% %                 % For simulation data, creating "nii" structure
% %                 nii = make_nii(MWF_Map);
% %             end 
% %         end          
% %         SaveMWFAsNifti(nii, MWF_Map, [Convert_2_Forwardslash(AnalysisFolder),Spatial_Filename_nii]);      
% %     end
end
%
%% 
function SaveMWFAsNifti(nii, MWF_Map, FilenameWtDir_MWF)
    % For 3D data
    nii.hdr.dime.bitpix=32; 
    nii.hdr.dime.datatype=16;
    nii.hdr.dime.dim(1)=3;
    nii.hdr.dime.dim(5)=1;     % nii.hdr.dime.dim(5) = image(4) = 3 or any suitable value
    nii.img = MWF_Map;
    save_untouch_nii(nii, FilenameWtDir_MWF);
end

%%
function [Vertex1, Vertex2] = Find_ROI_Coord_4_MaskedData(Mask3D)
    % If mask has NaN values, convert to zeros
    if any(isnan(Mask3D))
        Mask3D(isnan(Mask3D))=0.0;
    end
    %% THIS DOES NOT WORK FOR Z3
%     [Row3, Col3, z3] = find(abs(Mask3D-1)<eps);
%     figure(2)
%     plot(z3)

    %% THIS WORKS
    [SubV, SubH, SubZ] = ind2sub([size(Mask3D)],find(abs(Mask3D-1)<eps));
    Vertex1 = ones(1, 3);  Vertex2 = ones(1, 3); 
    Vertex1(1,1) = min(SubV);      Vertex1(1,2) = min(SubH);   Vertex1(1,3) = min(SubZ);
    Vertex2(1,1) = max(SubV);      Vertex2(1,2) = max(SubH);   Vertex2(1,3) = max(SubZ);
end

   %%
   function MultiDimVec = AddMarginToProcessedData(MultiDimVec, Sz3D0, Vert1, Vert2)
        Sz_MD =  size(MultiDimVec);
        Len = length(Sz_MD);   
       %
       if (Len==4)
           MultiDimVec0 = zeros(Sz3D0(1), Sz3D0(2), Sz3D0(3), Sz_MD(4));
           MultiDimVec0(Vert1(1):Vert2(1),Vert1(2):Vert2(2),Vert1(3):Vert2(3),:) = MultiDimVec;
       end
       %
       if (Len==3)
           MultiDimVec0 = zeros(Sz3D0(1), Sz3D0(2), Sz3D0(3));
           MultiDimVec0(Vert1(1):Vert2(1),Vert1(2):Vert2(2),Vert1(3):Vert2(3)) = MultiDimVec;
       end 
       %
       MultiDimVec = MultiDimVec0;
   end