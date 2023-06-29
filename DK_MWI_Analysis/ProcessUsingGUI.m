function varargout = ProcessUsingGUI(varargin)
    % ProcessUsingGUI MATLAB code for ProcessUsingGUI.fig
    %      ProcessUsingGUI, by itself, creates a new ProcessUsingGUI or raises the existing
    %      singleton*.
    %f
    %      H = ProcessUsingGUI returns the handle to a new ProcessUsingGUI or the handle to
    %      the existing singleton.
    %
    %      ProcessUsingGUI('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in ProcessUsingGUI.M with the given input arguments.
    %
    %      ProcessUsingGUI('Property','Value',...) creates a new ProcessUsingGUI or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before ProcessUsingGUI_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to ProcessUsingGUI_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help ProcessUsingGUI

    % Last Modified by GUIDE v2.5 18-Mar-2020 16:15:32

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @ProcessUsingGUI_OpeningFcn, ...
                       'gui_OutputFcn',  @ProcessUsingGUI_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
                   
    %%clc;      
                   
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT
end

% --- Executes just before ProcessUsingGUI is made visible.
function ProcessUsingGUI_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % varargin   command line arguments to ProcessUsingGUI (see VARARGIN)

%     % Choose default command line output for ProcessUsingGUI
%     handles.output = hObject;
%     
%     % Update handles structure
%     guidata(hObject, handles);
    
     global EchoNo_4_Mask  DicomFileExt;
     % global JavaPathStr1;
     
     %%% clc; 
     EchoNo_4_Mask = 1;  % 1st echo is used to make nifti mask using home-made code, similar to BET mask of FSL
     DicomFileExt    = '';       %'*.dcm'; 

    % DKU: Add code folder   
    [CodeFolder_Parent, CodeFilename, ext] = fileparts(mfilename('fullpath'));
    addpath(genpath(CodeFolder_Parent));   % Add all subfolders
    % JavaPathStr1 = [CodeFolder_Parent, '\Subcodes\ParforProgMonv2\java']; 
    handles.resetDK = false;
    guidata(hObject, handles); % Sharing data with other functions
    % Load handles is already saved. DKU
    cd(CodeFolder_Parent);    
    loadState(handles);  % DKU
end

% --- Outputs from this function are returned to the command line.
function varargout = ProcessUsingGUI_OutputFcn(hObject, eventdata, handles) 
%     % Get default command line output from handles structure
%     varargout{1} = handles.output;  % DKU commeted this out
    
%%        delete(hObject);  %% DKU not happy with this delete; else I get 2 GUIs.
end

%% Pushutton functions

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
    DialogTitle = 'Please Choose the Data Folder';
    
    StartFolder = pwd;   
    FolderName = uigetdir(StartFolder,DialogTitle);

    set(handles.edit1, 'string', FolderName); 
    cd(FolderName);
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
    DialogTitle = 'Please Choose the Data File';
    StartFolder = Try2StartWithRelevantFolder(handles.edit2);
    [FileName, FolderName, FilterIndex] = uigetfile({'*.nii';'*.mat'}, DialogTitle); 
    set(handles.edit2, 'string', [FolderName, FileName]); 
    cd(FolderName);
end

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
    DialogTitle = 'Please Choose the Data File';
    StartFolder = Try2StartWithRelevantFolder(handles.edit3);
    [FileName, FolderName, FilterIndex] = uigetfile({'*.nii';'*.mat'}, DialogTitle); 
    set(handles.edit3, 'string', [FolderName, FileName]); 
    cd(FolderName);
end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
    DialogTitle = 'Please Choose the Data File';
    StartFolder = Try2StartWithRelevantFolder(handles.edit4);
    [FileName, FolderName, FilterIndex] = uigetfile({'*.nii';'*.mat'}, DialogTitle); 
    set(handles.edit4, 'string', [FolderName, FileName]); 
    cd(FolderName);
end

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
    DialogTitle = 'Please Choose the Data File';
    StartFolder = Try2StartWithRelevantFolder(handles.edit5);
    [FileName, FolderName, FilterIndex] = uigetfile({'*.nii';'*.mat'}, DialogTitle); 
    set(handles.edit5, 'string', [FolderName, FileName]); 
    cd(FolderName);
end

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
    DialogTitle = 'Please Choose the Data File';
    StartFolder = Try2StartWithRelevantFolder(handles.edit6);
    [FileName, FolderName, FilterIndex] = uigetfile({'*.nii';'*.mat'}, DialogTitle); 
    set(handles.edit6, 'string', [FolderName, FileName]); 
    cd(FolderName);
end

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
%Read dicom data, average out and make a mat file
    global InputDataFilesCell;
    global EditBox7_Val EditBox8_Val EditBox9_Val ...
        EditBox11_Val EditBox12_Val EditBox13_Val EditBox14_Val EditBox15_Val ...
        EditBox16_Val EditBox17_Val EditBox18_Val EditBox42_Val EditBox56_Val EditBox57_Val EditBox58_Val
    global Popupmenu1_Str Popupmenu2_Str Popupmenu3_Str Popupmenu5_Str Popupmenu7_Str Popupmenu9_Str Popupmenu10_Str;
    global Popupmenu16_Str  Popupmenu18_Str Popupmenu20_Str;
% % % global Popupmenu6_Str 
% % %     global Popupmenu11_Str 
% % %     global Popupmenu12_Str;    
    global EditBox20_Val EditBox27_Val EditBox28_Val EditBox29_Val EditBox30_Val EditBox31_Val EditBox32_Val EditBox33_Val EditBox34_Val
    global EchoNo_4_Mask  DicomFileExt
    %global JavaPathStr1
    
    % Wait time before executing the code
    WaitHour_Val = str2num(get(handles.edit43, 'string'));    
    WaitMin_Val  = str2num(get(handles.edit44, 'string')); 
    TotalSeconds = WaitHour_Val*3600 + WaitMin_Val*60;
    pause(TotalSeconds); %in seconds  
    
   	[CodeFolder_Parent, CodeFilename, ext] = fileparts(mfilename('fullpath'));
    addpath(genpath(CodeFolder_Parent));   % Add all subfolders
   	JavaPathStr1 = Convert_2_Forwardslash([CodeFolder_Parent, '\Subcodes\ParforProgMonv2\java']); 

    % Choose solver
    Solver_Type = '2';    % '1' --> My dafult spatial solver, '2' --> fast spatial solver

    % These two function helps access all values in GUI
    GetAllPar_in_ProcessingParameterPanel(handles);
    GetAllDataFilenames_in_InputDataFilesPanel(handles);    

    Capsule1.FAE_in   = EditBox7_Val;
    Capsule1.FAE_Fin  = EditBox8_Val;
    Capsule1.FAE_step = EditBox9_Val;
    
    %assigning the values to proper fields
    SizeX1_DSW = EditBox11_Val;
    SizeY1_DSW = EditBox12_Val;
    SizeZ1_DSW = EditBox20_Val;   
    Size_DSW = [SizeX1_DSW, SizeY1_DSW, SizeZ1_DSW];

    ResX                    = EditBox27_Val;
    ResY                    = EditBox28_Val;
    ResZ                    = EditBox29_Val;
    Resol                   = [ResX, ResY, ResZ]; 

    WOR                     = EditBox13_Val;
    Capsule1.alphaEx_NominalRad      = (double(EditBox14_Val)/90)*(pi/2);
    Capsule1.alphaRef_NominalRad     = (double(EditBox15_Val)/180)*pi;
    Capsule1.nMatlabPoolWorkers_Conv = EditBox16_Val;
    Capsule1.nMatlabPoolWorkers_Spat = EditBox42_Val;    
    
    Capsule1.No_of_Iterations        = EditBox17_Val;
    Capsule1.No_of_MixedIterations        = EditBox56_Val;
    
    Capsule1.MixedIter.WinSz2D  = [EditBox57_Val, EditBox58_Val];
    
    Capsule1.alpha1                  = EditBox18_Val;    
    Capsule1.Myelin_Cutoff           = [EditBox30_Val; EditBox31_Val]*1e-3;   % in sec   
    Capsule1.T2_Initial              = EditBox32_Val;   % in ms   
    Capsule1.T2_Final                = EditBox33_Val;   % in ms   
    Capsule1.nT2                     = EditBox34_Val;   % in ms   
    

    FAECorrection_Tag           = Popupmenu1_Str; 
    UseMask                 = Popupmenu2_Str;
    RegTag                  = Popupmenu3_Str;  
    SelectEchoTag           = Popupmenu5_Str;
    Spatial_ProcessMode     = Popupmenu7_Str;
    InputDataMode           = Popupmenu9_Str;
    ImageTypeTag            = Popupmenu10_Str; 
%     Capsule1.TempoalAlgoTag = Popupmenu11_Str;
%    muS_Opt_TagTag   = Popupmenu12_Str;

    Capsule1.muS_OptRefine_Tag           = Popupmenu18_Str; 
    Capsule1.muFAI_OptRefine_Tag         = Popupmenu20_Str; 

    Capsule1.JavaPathStr1 = JavaPathStr1;
    Capsule1.SelectEchoTag = SelectEchoTag;
    Capsule1.FAECorrection_Tag = FAECorrection_Tag;    
    
    %UseMask has to be converted into 'Yes' and 'No' to be consistent with nonGUI code
    if strcmp(UseMask, 'Use Mask - automatic')
        UseMask     = 'Yes';
        MaskType    = 'nii';                % 'nii',  'roi', 'none'  
        MaskFileNamePrefix = '';
        MaskExt     = '*mask.nii*';             % could be .nii or nii.gz    
    elseif strcmp(UseMask, 'Use Mask - *.ROI file')
        UseMask = 'Yes';
        MaskType = 'roi';                % 'nii',  'roi', 'none'    
        MaskFileNamePrefix = 'Mask';
        MaskExt     = '*.roi';
    elseif strcmp(UseMask, 'Use Mask - *.NII file')
        UseMask     = 'Yes';
        MaskType    = 'nii';                % 'nii',  'roi', 'none'  
        MaskFileNamePrefix = '';
        MaskExt     = '*mask.nii*';             % could be .nii or nii.gz
    else
        UseMask     = 'No';
        MaskType    = 'none';                % 'nii',  'roi', 'none' 
    end
        
    %Find only non-empty entries in the cellf
    IND = find(~cellfun(@isempty,InputDataFilesCell));
    InputDataFilesCell0 = InputDataFilesCell(IND);
    
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

    % Changing some values in GUI fields for consistencies
    if strcmpi(Popupmenu16_Str,'Use set value')
        set(handles.edit18,'Visible','On')
    else
        set(handles.edit18,'Visible','Off')
    end
    %
    if strcmpi(Popupmenu1_Str, 'WITHOUT any FAE correction')
        set(handles.edit17, 'string',num2str(1));
    end

    for i = 1:length(InputDataFilesCell0)
        DicomFileDir = squeeze(InputDataFilesCell0{i});

        if strcmp(InputDataMode, 'Raw Dicom Folder')
            OutputFileDir = QT2R_Read_N_Avg_N_MakeMask_NII(DicomFileDir, DicomFileExt, EchoNo_4_Mask, ImageType_Desired);  % Making mask without any need for FSL
            Filename_wt_Folder = [Convert_2_Forwardslash(OutputFileDir), 'Avg_T2RelaxometryData.nii'];
            
        elseif strcmp(InputDataMode, 'Nifti-4D-Kumar-Format')
            OutputFileDir = DicomFileDir;
            Filename_wt_Folder = [Convert_2_Forwardslash(OutputFileDir), 'Avg_T2RelaxometryData.nii'];
        else
            OutputFileDir = DicomFileDir;
            Filename_wt_Folder = [Convert_2_Forwardslash(OutputFileDir), 'Avg_T2RelaxometryData.mat'];
        end

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
            if strcmp(SelectEchoTag, 'EvenEchoes')
                Temporal_Prefix = 'Temporal_EvenEchoes_';   
                Spatial_Prefix = 'Spatial_EvenEchoes_'; 
            elseif strcmp(SelectEchoTag, 'AllEchoes_Except_1st') 
                Temporal_Prefix = 'Temporal_Exc_1stEchoes_';   
                Spatial_Prefix = 'Spatial_Exc_1stEchoes_';                
            else  % strcmp(SelectEchoTag, 'AllEchoes') 
                Temporal_Prefix = 'Temporal_AllEchoes_';   
                Spatial_Prefix = 'Spatial_AllEchoes_'; 
            end

            AnalysisFolder = DataFolder;        
            MaskFolder = DataFolder;
            contents = cellstr(get(handles.popupmenu3,'String'));
            RegularizationTag = contents{get(handles.popupmenu3,'Value')}; 
            
            contents = cellstr(get(handles.popupmenu8,'String'));
            MatPriorityTag = contents{get(handles.popupmenu8,'Value')}; 
            
            % Some variables in GUI fields are found to be cleared when
            % multiple data are processed; so, reading everything back:
            GetAllPar_in_ProcessingParameterPanel(handles);
                
            
            % EvenEchoes tag not allowed for stimulated correction.
            if strcmp(SelectEchoTag, 'EvenEchoes')
                error('Error by DK. \n EvenEchoes tag not allowed for stimulated correction. Please correct your selection and then run')
            end

            % Open file for debug logging
            fileID = fopen(fullfile(DataFolder, 'Report_Debug_DK.txt'),'w')
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
            Capsule1.alphaS_DeltaVec = double(logspace(1,5,8));                                         
            % Initial values for determining muS (spatial regularizatio const for T2 dist. map) 
            %Capsule1.alphas_Info.alpha1	= [100, 200, 500, 1000, 2000, 4000, 7500, 10000];      
            Capsule1.alphas_Info.alpha1     = [10, 100, 200, 350, logspace(log10(5e2), log10(3.5e3), 10), 5e3, 1e4];    
           %%
           % If both temporal and spatial regularization ask for the
           % same number of Matpool worker, then open matpool only once
           if (Capsule1.nMatlabPoolWorkers_Conv == Capsule1.nMatlabPoolWorkers_Spat)
                % Open Matlab pool and set priority
                MatpoolPropSpat = OpenMatPool_N_SetPriority(MatpoolPropSpat);  % Conventional regularization is not that computational demanding; so, speed up is possble using virtual cores.  
           end

           Capsule1.NoRegulMode = 'Yes';               % 'Yes', 'No'
           % Performing temporal regularization ONLY
           debug_point(fileID)
           if (strcmp(RegTag, 'Both') || strcmp(RegTag, 'Temporal-Only'))  %'Temporal', 'Spatial', 'Both'
                % Open Matlab pool and set priority
                MatpoolPropConv = OpenMatPool_N_SetPriority(MatpoolPropConv);  % Conventional regularization is not that computational demanding; so, speed up is possble using virtual cores.  

                RegTag1 = 'Temporal-Only';  % Temporarily rewriting RegTag1
                % Capsule1.alphas_Info.alpha1 = 1e-6;    % This is a dummy value, not used for temporal regularizatio. But, the code needs it.
                muS_Opt_Tag = 'No'; 	  % 'Yes', 'No';   This is a dummy value, not used for temporal regularizatio. But, the code needs it.
                [ResVec, NormSpaVec] = QT2RProcess_wt_ParforMonitor_Stim3D(DataFolder, DataFileName, Size_DSW, WOR, Resol, MaskFileName, Spatial_ProcessMode, RegTag1, MatPriorityTag, Capsule1, muS_Opt_Tag);
           end
           
           %%
           %
           debug_point(fileID)
           % Find OUT if mus optimization is desired or not
           contents = cellstr(get(handles.popupmenu16,'String'));
           Popupmenu16_Str = contents{get(handles.popupmenu16,'Value')}; 
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
                RegTag1 = 'Spatial-Only'; % Temporarily rewriting RegTag1
                Capsule1.Popupmenu16_Str = Popupmenu16_Str; % Only to ensure that the variable stays even after being cleared by child 
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

               % Load residuals for unregularized case
                if strcmpi(FAECorrection_Tag, 'WITH FAE correction')
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

                %%
                if strcmpi(Popupmenu16_Str,'Use set value')
                    set(handles.edit18,'Visible','On')
                else
                    set(handles.edit18,'Visible','Off')
                end
                %
                Capsule1.NoRegulMode = 'No';               % 'Yes', 'No'
                debug_point(fileID)

                if strcmpi(Capsule1.Popupmenu16_Str,'Optimize muS')
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
                    debug_point(fileID)
                    Capsule1_temp.alphaS_DeltaVec = Capsule1.alphaS_DeltaVec;
                    Capsule1_temp.alphas_Info.alpha1 = Capsule1.alphas_Info.alpha1;

                    AlphaS_OPt                  = FindOptAlphaS(DataFolder, DataFileName, Size_DSW, WOR, Resol, MaskFileName, Spatial_ProcessModeTemp, RegTag1, MatPriorityTag, Capsule1_temp, muS_Opt_Tag, muS_OptRefine_Tag, fileID);
                    Capsule1.alphas_Info.alpha1 = AlphaS_OPt;
                    fprintf('Calculation ended  : normalizing factor for muS optimization \n')
                    debug_point(fileID)
                else % Pick up the supplied value of muS
                    Capsule1.alphas_Info.alpha1 = str2num(get(handles.edit18, 'string')); % AlphaS_OPt;
                end

                %%
                %
                % Now process entire data set
                % muS_Opt_Tag = 'Yes' ensures that there are no optimization over muS_Delta_Opt
                muS_Opt_Tag = 'No'; 	% 'Yes', 'No';
                %Capsule1.alphas_Info.alpha1 = 1000;
                d = debug_point(fileID)
                fprintf(fileID,'Optimization of muS parameters done;Still has to perform: muS_Delta optimization+ spatially regularized solution.');
                [ResVec, NormSpaVec] = QT2RProcess_wt_ParforMonitor_Stim3D(DataFolder, DataFileName, Size_DSW, WOR, Resol, MaskFileName, Spatial_ProcessMode, RegTag1, MatPriorityTag, Capsule1, muS_Opt_Tag);
                debug_point(fileID)
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
    
    %Priority set back to normal (Second way of priority setting.)
    if (~isunix)
        try % Even if it does not run, no big deal
            SetMatlabPriorityWindowDKU('normal'); 
        end
    end
    %clearvars;  % Clear all variables active in workspace
end

%% other GUIfunctions
% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

    global prevValIteration;
    % DK Changing some values in GUI fields for consistencies
    contents = cellstr(get(hObject,'String'));
    Popupmenu1_Str = contents{get(hObject,'Value')}; 
    if strcmpi(Popupmenu1_Str,'WITH FAE correction')
        % If value is 1, check for the previous value
        if (str2num(get(handles.edit17, 'string'))==1)
            if ~isempty(prevValIteration)
                set(handles.edit17,'string',prevValIteration);
            end
        end
    else
        prevValIteration = str2num(get(handles.edit17, 'string'));
        set(handles.edit17,'string',num2str(1));
    end
end

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    FAECorrectionOpt = {'WITH FAE correction', 'WITHOUT any FAE correction'};
    set(hObject, 'String', FAECorrectionOpt);
end

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
end

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    MaskType = {'Use Mask - automatic', 'Use Mask - *.NII file', 'Use Mask - *.ROI file', 'No Mask'};
    set(hObject, 'String', MaskType);
    
end

% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3

end

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    RegularizationTagOpt = {'Both', 'Temporal-Only', 'Spatial-Only'};
    set(hObject, 'String', RegularizationTagOpt);
    
end

% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    TempoEchoTag = {'AllEchoes', 'EvenEchoes', 'AllEchoes_Except_1st'};
    set(hObject, 'String', TempoEchoTag);
end


% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    FirstDiffOPTag = {'3D', '2D_3D_Mixed', '2D'};
    set(hObject, 'String', FirstDiffOPTag); 
end

% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    MatlabPriorityTag = {'normal', 'low', 'below normal', 'above normal', 'high priority', 'real time'};
    set(hObject, 'String', MatlabPriorityTag); 
end

% --- Executes during object creation, after setting all properties.
function popupmenu9_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    DataInputTag = {'Raw Dicom Folder', 'Nifti-4D-Kumar-Format', '.mat-4D-Kumar-Format'};
    set(hObject, 'String', DataInputTag); 
end

% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    ImageType_DesiredTag = {'Philips_Magn: ORIGINAL\PRIMARY\M_SE\M\SE', 'Siemens_Magn: ORIGINAL\PRIMARY\M\ND', 'Philips_Magn_2Slabs: ORIGINAL\PRIMARY\M_IR\M\IR'};
    set(hObject, 'String', ImageType_DesiredTag); 
end

%% three customized GUI fns
function GetAllPar_in_ProcessingParameterPanel(handles)
global EditBox7_Val EditBox8_Val EditBox9_Val ...
    EditBox11_Val EditBox12_Val EditBox13_Val EditBox14_Val EditBox15_Val ...
    EditBox16_Val EditBox17_Val EditBox18_Val ...
    EditBox20_Val EditBox27_Val EditBox28_Val EditBox29_Val EditBox30_Val EditBox31_Val ...
    EditBox32_Val EditBox33_Val EditBox34_Val EditBox42_Val EditBox56_Val EditBox57_Val EditBox58_Val

global Popupmenu1_Str Popupmenu2_Str Popupmenu3_Str Popupmenu5_Str;
global Popupmenu7_Str Popupmenu9_Str Popupmenu10_Str Popupmenu18_Str Popupmenu20_Str;
% % % global Popupmenu6_Str 
% % % global Popupmenu11_Str 
% % % global Popupmenu12_Str;

    EditBox7_Val = str2num(get(handles.edit7, 'string')); 
    EditBox8_Val = str2num(get(handles.edit8, 'string')); 
    EditBox9_Val = str2num(get(handles.edit9, 'string')); 
    %
    EditBox11_Val = str2num(get(handles.edit11, 'string')); 
    EditBox12_Val = str2num(get(handles.edit12, 'string')); 
    EditBox13_Val = str2num(get(handles.edit13, 'string')); 
    EditBox14_Val = str2num(get(handles.edit14, 'string')); 
    EditBox15_Val = str2num(get(handles.edit15, 'string'));
    EditBox16_Val = str2num(get(handles.edit16, 'string'));
   
    EditBox17_Val = str2num(get(handles.edit17, 'string')); 
    EditBox56_Val = str2num(get(handles.edit56, 'string')); 
    EditBox18_Val = str2num(get(handles.edit18, 'string'));         

    EditBox20_Val = str2num(get(handles.edit20, 'string')); 
    EditBox27_Val = str2num(get(handles.edit27, 'string'));
    EditBox28_Val = str2num(get(handles.edit28, 'string'));
    EditBox29_Val = str2num(get(handles.edit29, 'string'));
    EditBox30_Val = str2num(get(handles.edit30, 'string'));
    EditBox31_Val = str2num(get(handles.edit31, 'string'));
    
    EditBox32_Val = str2num(get(handles.edit32, 'string'));
    EditBox33_Val = str2num(get(handles.edit33, 'string'));
    EditBox34_Val = str2num(get(handles.edit34, 'string'));
    EditBox42_Val = str2num(get(handles.edit42, 'string')); 
    
    EditBox57_Val = str2num(get(handles.edit57, 'string'));  
    EditBox58_Val = str2num(get(handles.edit58, 'string')); 
    
    %Get the POPUPMENUs' values
    contents = cellstr(get(handles.popupmenu1,'String'));
    Popupmenu1_Str = contents{get(handles.popupmenu1,'Value')};
    
    contents = cellstr(get(handles.popupmenu2,'String'));
    Popupmenu2_Str = contents{get(handles.popupmenu2,'Value')};
    
    contents = cellstr(get(handles.popupmenu3,'String'));
    Popupmenu3_Str = contents{get(handles.popupmenu3,'Value')};    
    
    contents = cellstr(get(handles.popupmenu5,'String'));
    Popupmenu5_Str = contents{get(handles.popupmenu5,'Value')};
    
% % %     contents = cellstr(get(handles.popupmenu6,'String'));
% % %     Popupmenu6_Str = contents{get(handles.popupmenu6,'Value')};     
    
    contents = cellstr(get(handles.popupmenu7,'String'));
    Popupmenu7_Str = contents{get(handles.popupmenu7,'Value')};     

    contents = cellstr(get(handles.popupmenu9,'String'));
    Popupmenu9_Str = contents{get(handles.popupmenu9,'Value')};  
    
    contents = cellstr(get(handles.popupmenu10,'String'));
    Popupmenu10_Str = contents{get(handles.popupmenu10,'Value')};  
    
% % %     contents = cellstr(get(handles.popupmenu11,'String'));
% % %     Popupmenu11_Str = contents{get(handles.popupmenu11,'Value')};

    contents = cellstr(get(handles.popupmenu18,'String'));
    Popupmenu18_Str = contents{get(handles.popupmenu18,'Value')};  
    
    contents = cellstr(get(handles.popupmenu20,'String'));
    Popupmenu20_Str = contents{get(handles.popupmenu20,'Value')};  
end

function GetAllDataFilenames_in_InputDataFilesPanel(handles)
global InputDataFilesCell;
InputDataFilesCell = cell(6,1);
    InputDataFilesCell{1,1} = get(handles.edit1, 'string'); 
    InputDataFilesCell{2,1} = get(handles.edit2, 'string'); 
    InputDataFilesCell{3,1} = get(handles.edit3, 'string'); 
    InputDataFilesCell{4,1} = get(handles.edit4, 'string'); 
    InputDataFilesCell{5,1} = get(handles.edit5, 'string');
    InputDataFilesCell{6,1} = get(handles.edit6, 'string');
    %
    contents = cellstr(get(handles.popupmenu9,'String'));
    InputDataMode = contents{get(handles.popupmenu9,'Value')};
    %Start processing data
    if strcmp(InputDataMode, 'Raw Dicom Folder')       
        % If the last subfolder is not "PreprocessedData_N_Analysis"
        FolderName = EnsureTheCurrentDataFolder_Post_DataConv(InputDataFilesCell{1,1});
        set(handles.edit1, 'string', FolderName); 
        FolderName = EnsureTheCurrentDataFolder_Post_DataConv(InputDataFilesCell{2,1});
        set(handles.edit2, 'string', FolderName);
        FolderName = EnsureTheCurrentDataFolder_Post_DataConv(InputDataFilesCell{3,1});
        set(handles.edit3, 'string', FolderName);
        FolderName = EnsureTheCurrentDataFolder_Post_DataConv(InputDataFilesCell{4,1});
        set(handles.edit4, 'string', FolderName);
        FolderName = EnsureTheCurrentDataFolder_Post_DataConv(InputDataFilesCell{5,1});
        set(handles.edit5, 'string', FolderName);
        FolderName = EnsureTheCurrentDataFolder_Post_DataConv(InputDataFilesCell{6,1});
        set(handles.edit6, 'string', FolderName);
    end
end

function Folder = EnsureTheCurrentDataFolder_Post_DataConv(Folder)
    [~, SubFolder, ~] = fileparts(Folder);
    if ((~strcmpi(SubFolder,'PreprocessedData_N_Analysis')) && (~isempty(Folder)))
        if strcmpi(Folder(end), filesep)
            Folder = Folder(1:end-1);
        end
        Folder = [Folder, filesep, 'PreprocessedData_N_Analysis'];
    end
end

function StartFolder = Try2StartWithRelevantFolder(handles_2_Editbox)
    StartFolder = get(handles_2_Editbox, 'String')
    if isempty(StartFolder)
        StartFolder = pwd;
    else
        %check if StartFolder exists
        if ~(exist(StartFolder, 'file') == 7) % folder
            DialogTitle = 'Please Choose the MatlabUserLib folder';
            StartFolder = pwd;   
            FolderName = uigetdir(StartFolder,DialogTitle);
            set(handles_2_Editbox, 'string', FolderName); 
            cd(FolderName);
        else
            cd(StartFolder)
        end
    end
end

%% DK: some auxiliary non GUI codes
% function Folder = Convert_2_Forwardslash(Folder)
%     if ~strcmp(Folder(end), filesep)
%         Folder = [Folder, filesep];
%     end
%     %Convert to forward slash and also ends with forward slash    
%     Folder = strrep(Folder, '\', '/');  
% end

function MaskFileName = findFile_with_roi_OR_nii_Ext(MaskFolder,  MaskFileNamePrefix, MaskExt)
    %ensure that MaskFolder is ended with filesep
    MaskFolder = Convert_2_Forwardslash(MaskFolder)
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

function saveState(handles) 
    if ~isfield(handles, 'resetDK')
        handles.resetDK = false;
    end
    % Now save state
    if (handles.resetDK == false)
        CapsuleGUI.edit7_str = get(handles.edit7, 'string'); 
        CapsuleGUI.edit8_str = get(handles.edit8, 'string'); 
        CapsuleGUI.edit9_str = get(handles.edit9, 'string'); 
        %
        CapsuleGUI.edit11_str = get(handles.edit11, 'string'); 
        CapsuleGUI.edit12_str = get(handles.edit12, 'string'); 
        CapsuleGUI.edit13_str = get(handles.edit13, 'string'); 
        CapsuleGUI.edit14_str = get(handles.edit14, 'string'); 
        CapsuleGUI.edit15_str = get(handles.edit15, 'string');
        CapsuleGUI.edit16_str = get(handles.edit16, 'string');
        CapsuleGUI.edit17_str = get(handles.edit17, 'string');    
        CapsuleGUI.edit18_str = get(handles.edit18, 'string');         
    % %     CapsuleGUI.edit19_str = get(handles.edit19, 'string'); 
        CapsuleGUI.edit20_str = get(handles.edit20, 'string'); 
        CapsuleGUI.edit27_str = get(handles.edit27, 'string');
        CapsuleGUI.edit28_str = get(handles.edit28, 'string');
        CapsuleGUI.edit29_str = get(handles.edit29, 'string');
        CapsuleGUI.edit30_str = get(handles.edit30, 'string');
        CapsuleGUI.edit31_str = get(handles.edit31, 'string');
        CapsuleGUI.edit32_str = get(handles.edit32, 'string');
        CapsuleGUI.edit33_str = get(handles.edit33, 'string');
        CapsuleGUI.edit34_str = get(handles.edit34, 'string');  
        CapsuleGUI.edit42_str = get(handles.edit42, 'string'); 
        CapsuleGUI.edit56_str = get(handles.edit56, 'string');  
        CapsuleGUI.edit57_str = get(handles.edit57, 'string');  
        CapsuleGUI.edit58_str = get(handles.edit58, 'string');  

        CapsuleGUI.popupmenu1_Val = get(handles.popupmenu1,'Value');
        CapsuleGUI.popupmenu2_Val = get(handles.popupmenu2,'Value');
        CapsuleGUI.popupmenu3_Val = get(handles.popupmenu3,'Value');
        CapsuleGUI.popupmenu5_Val = get(handles.popupmenu5,'Value');
    % % %     CapsuleGUI.popupmenu6_Val = get(handles.popupmenu6,'Value');
        CapsuleGUI.popupmenu7_Val = get(handles.popupmenu7,'Value');
        CapsuleGUI.popupmenu8_Val = get(handles.popupmenu8,'Value');
        CapsuleGUI.popupmenu9_Val = get(handles.popupmenu9,'Value');
        CapsuleGUI.popupmenu10_Val = get(handles.popupmenu10,'Value');
    % % %     CapsuleGUI.popupmenu11_Val = get(handles.popupmenu11,'Value');
    % % %     CapsuleGUI.popupmenu12_Val = get(handles.popupmenu12,'Value');
        CapsuleGUI.popupmenu16_Val = get(handles.popupmenu16,'Value');
        CapsuleGUI.popupmenu18_Val = get(handles.popupmenu18,'Value');
        CapsuleGUI.popupmenu20_Val = get(handles.popupmenu20,'Value');
        % write the state to file
        [CodeFolder_Parent, CodeFilename, ext] = fileparts(mfilename('fullpath'));
        %
        save([Convert_2_Forwardslash(CodeFolder_Parent), 'state.mat'], 'CapsuleGUI');
    end
end

function loadState(handles)
    filename = 'state.mat';
    [CodeFolder_Parent, CodeFilename, ext] = fileparts(mfilename('fullpath'));
        
    if exist([Convert_2_Forwardslash(CodeFolder_Parent), 'state.mat'],'file')
        load(filename);
        set(handles.edit7, 'string', CapsuleGUI.edit7_str); 
        set(handles.edit8, 'string', CapsuleGUI.edit8_str); 
        set(handles.edit9, 'string', CapsuleGUI.edit9_str); 
        %
        set(handles.edit11, 'string', CapsuleGUI.edit11_str); 
        set(handles.edit12, 'string', CapsuleGUI.edit12_str); 
        set(handles.edit13, 'string', CapsuleGUI.edit13_str); 
        set(handles.edit14, 'string', CapsuleGUI.edit14_str); 
        set(handles.edit15, 'string', CapsuleGUI.edit15_str); 
        set(handles.edit16, 'string', CapsuleGUI.edit16_str); 
        set(handles.edit17, 'string', CapsuleGUI.edit17_str); 
        set(handles.edit18, 'string', CapsuleGUI.edit18_str); 
% %         set(handles.edit19, 'string', CapsuleGUI.edit19_str); 
        set(handles.edit20, 'string', CapsuleGUI.edit20_str); 
        set(handles.edit27, 'string', CapsuleGUI.edit27_str); 
        set(handles.edit28, 'string', CapsuleGUI.edit28_str); 
        set(handles.edit29, 'string', CapsuleGUI.edit29_str);
        set(handles.edit30, 'string', CapsuleGUI.edit30_str);
        set(handles.edit31, 'string', CapsuleGUI.edit31_str);
        set(handles.edit32, 'string', CapsuleGUI.edit32_str);
        set(handles.edit33, 'string', CapsuleGUI.edit33_str);
        set(handles.edit34, 'string', CapsuleGUI.edit34_str);
        set(handles.edit42, 'string', CapsuleGUI.edit42_str);
        set(handles.edit56, 'string', CapsuleGUI.edit56_str);
        set(handles.edit57, 'string', CapsuleGUI.edit57_str);
        set(handles.edit58, 'string', CapsuleGUI.edit58_str);
%%%%%%        
        set(handles.popupmenu1, 'Value', CapsuleGUI.popupmenu1_Val);
        set(handles.popupmenu2, 'Value', CapsuleGUI.popupmenu2_Val);
        set(handles.popupmenu3, 'Value', CapsuleGUI.popupmenu3_Val);
        set(handles.popupmenu5, 'Value', CapsuleGUI.popupmenu5_Val);
% % %         set(handles.popupmenu6, 'Value', CapsuleGUI.popupmenu6_Val);
        set(handles.popupmenu7, 'Value', CapsuleGUI.popupmenu7_Val);
        set(handles.popupmenu8, 'Value', CapsuleGUI.popupmenu8_Val);
        set(handles.popupmenu9, 'Value', CapsuleGUI.popupmenu9_Val);
        set(handles.popupmenu10, 'Value', CapsuleGUI.popupmenu10_Val);
% % %         set(handles.popupmenu11, 'Value', CapsuleGUI.popupmenu11_Val);
% % %         set(handles.popupmenu12, 'Value', CapsuleGUI.popupmenu12_Val);
        set(handles.popupmenu16, 'Value', CapsuleGUI.popupmenu16_Val);
        set(handles.popupmenu18, 'Value', CapsuleGUI.popupmenu18_Val);
        set(handles.popupmenu20, 'Value', CapsuleGUI.popupmenu20_Val);
        % Findout the contents of handles.popupmenu16
        contents = cellstr(get(handles.popupmenu16,'String'));
        Popupmenu16_Str = contents{get(handles.popupmenu16,'Value')}; 

        % Changing some values in GUI fields for consistencies
        if strcmpi(Popupmenu16_Str,'Use set value')
            set(handles.edit18,'Visible','On');
            set(handles.text12,'Visible','On');
        else
            set(handles.edit18,'Visible','Off');
            set(handles.text12,'Visible','Off');
        end
    end
end

% --- Executes when user attempts to close figure1 / GUI.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
    % Load handles is already saved DKumar
    saveState(handles);
    %Object is deleted
    delete(hObject);
end


% --- Executes on button press in pushbutton16.
% This only read the dicom data and convert to nifti.
% This will not participate in MWI analysis.
function pushbutton16_Callback(hObject, eventdata, handles)
    global InputDataFilesCell;
    global Popupmenu10_Str  Popupmenu9_Str;
    global EchoNo_4_Mask  DicomFileExt;
    global JavaPathStr1;
     
    % Get all data folders
    GetAllDataFilenames_in_InputDataFilesPanel(handles);  
    % Find important fields
    contents = cellstr(get(handles.popupmenu8,'String'));
    MatPriorityTag = contents{get(handles.popupmenu8,'Value')}; 
    contents = cellstr(get(handles.popupmenu9,'String'));
    Popupmenu9_Str = contents{get(handles.popupmenu9,'Value')};
    contents = cellstr(get(handles.popupmenu10,'String'));
    Popupmenu10_Str = contents{get(handles.popupmenu10,'Value')};
            
    MatpoolPropRead.nPoolMax        = FindMaxMatpool()-1;
    MatpoolPropRead.MatPriorityTag  = MatPriorityTag;     
    MatpoolPropRead = OpenMatPool_N_SetPriority(MatpoolPropRead);
    
    %Find only non-empty entries in the cell
    IND = find(~cellfun(@isempty,InputDataFilesCell));
    InputDataFilesCell0 = InputDataFilesCell(IND);
    InputDataMode           = Popupmenu9_Str;
        
    %   ImageType_Desired should be changed based on scanner and sequence type  
    ImageTypeTag            = Popupmenu10_Str; 
    if strcmp(ImageTypeTag, 'Philips_Magn: ORIGINAL\PRIMARY\M_SE\M\SE')
        ImageType_Desired = 'ORIGINAL\PRIMARY\M_SE\M\SE';       % Single slab expt:  'ORIGINAL\PRIMARY\M_SE\M\SE', 'ORIGINAL\PRIMARY\PHASE MAP\P\SE' for philips data
    elseif strcmp(ImageTypeTag, 'Siemens_Magn: ORIGINAL\PRIMARY\M\ND')
        ImageType_Desired = 'ORIGINAL\PRIMARY\M\ND';            % Single slab expt:  'ORIGINAL\PRIMARY\M\ND' for siemens data
    elseif strcmp(ImageTypeTag, 'Philips_Magn_2Slabs: ORIGINAL\PRIMARY\M_IR\M\IR')
        ImageType_Desired = 'ORIGINAL\PRIMARY\M_IR\M\IR';       % DOUBLE slab expt:  'ORIGINAL\PRIMARY\M_IR\M\IR', 'ORIGINAL\PRIMARY\PHASE MAP\P\IR' for philips data   
    else
        % This ImageType_Desired NOT possible as chosen from dropdown menu.
    end

    %Start processing data
    if strcmp(InputDataMode, 'Raw Dicom Folder')
        for i = 1:length(InputDataFilesCell0)
            DicomFileDir = squeeze(InputDataFilesCell0{i});
            % Reading data and making mask
            OutputFileDir = QT2R_Read_N_Avg_N_MakeMask_NII(DicomFileDir, DicomFileExt, EchoNo_4_Mask, ImageType_Desired);  % Making mask without any need for FSL
        end 
    else
        display('''Data Input Option'' should be ''Raw Dicom Folder''' );
        display('Please choose so');
        return
    end
end


% --- Executes during object creation, after setting all properties.
function popupmenu16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    OptimizatioTag = {'Optimize muS', 'Use set value'};
    set(hObject, 'String', OptimizatioTag);
end

% --- Executes on selection change in popupmenu3.
function popupmenu16_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3

        contents = cellstr(get(hObject,'String'));
        Popupmenu16_Str = contents{get(hObject,'Value')}; 

        % Changing some values in GUI fields for consistencies
        if strcmpi(Popupmenu16_Str,'Use set value')
            set(handles.edit18,'Visible','On');
            set(handles.text12,'Visible','On');
        else
            set(handles.edit18,'Visible','Off');
            set(handles.text12,'Visible','Off');
        end
end


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    %
    [CodeFolder_Parent, CodeFilename, ext] = fileparts(mfilename('fullpath'));
    % Delete file
    handles.resetDK = true;
    guidata(hObject, handles); % Sharing data with other functions
    delete(fullfile(CodeFolder_Parent,'state.mat'));  
end

% --- Executes during object creation, after setting all properties.
function popupmenu18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    OptimizatioTag = {'Sufficient', 'Refined'};
    set(hObject, 'String', OptimizatioTag);
end

% --- Executes during object creation, after setting all properties.
function popupmenu20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    OptimizatioTag = {'Sufficient', 'Refined'};
    set(hObject, 'String', OptimizatioTag);
end
