%%
function AlphaS_OPt = FindOptAlphaS(DataFolder, DataFileName, Size_DSW, WOR, Resol, MaskFileName, Spatial_ProcessMode, RegTag, MatPriorityTag, Capsule1, muS_Opt_Tag, muS_OptRefine_Tag, fileID)
    debug_point(fileID);
    fprintf(fileID, 'time: %s', datestr(now,'mm-dd-yyyy HH-MM'));
    
    ErrorSR_T2Dist = 25;         % extra percentage error allowed for spatial regularization of T2Dist
    PlotTag = 'Off';            % 'On', 'Off';   Whether it will plot L-curve or not
    
    JavaPathStr1                = Capsule1.JavaPathStr1;
    Capsule1.No_of_Iterations   = 1;
    
    % What should be my reference residual value
    Res_0       = Capsule1.alphas_Info.Res_0;     
    Res_Tempo   = Capsule1.alphas_Info.Res_Tempo;  
    if (Res_Tempo >= 1.5*Res_0)
        Res_0 = Res_Tempo;
    end
    %
    % Reconstruct for AlphaS_Vec
    [Res1D_SpatRegVec, NormSpa_SpatRegVec] = QT2RProcess_wt_ParforMonitor_Stim3D(DataFolder, DataFileName, Size_DSW, WOR, Resol, MaskFileName, Spatial_ProcessMode, RegTag, MatPriorityTag, Capsule1, muS_Opt_Tag);
    AlphaS = Capsule1.alphas_Info.alpha1; 
    debug_point(fileID);
%     %
%     % To avoid oversmoothing, we restrict to indices of normDataTermFAE_Vec < (1+ErrorSR_T2Dist/100)*resFAE0 and use only a subset of data
%     AlphaS              = AlphaS            (Res1D_SpatRegVec < (1+ErrorSR_T2Dist/100).*Res_0);
%     NormSpa_SpatRegVec	= NormSpa_SpatRegVec(Res1D_SpatRegVec < (1+ErrorSR_T2Dist/100).*Res_0);
%     Res1D_SpatRegVec    = Res1D_SpatRegVec  (Res1D_SpatRegVec < (1+ErrorSR_T2Dist/100).*Res_0);    
%     

    debug_point(fileID);
    display('First instance of L-curve')   
    min_nEdist = Find_Regularization_Par(Res1D_SpatRegVec, NormSpa_SpatRegVec, PlotTag);  % If only narrow range considered
    debug_point(fileID);
    %% Scan some more
    if strcmpi(muS_OptRefine_Tag, 'Refined')  
        % Plot L-curve and find three alphaS with lowest distance to origin
        %min3_nEdist = Find_Regularization_Par3(Res1D_SpatRegVec, NormSpa_SpatRegVec, PlotTag);
        min3_nEdist = Find_3Possible_LCorners(Res1D_SpatRegVec, NormSpa_SpatRegVec, AlphaS, PlotTag);
        AlphaS_Sep  = (AlphaS(min3_nEdist(1))-AlphaS(min3_nEdist(2)));           
        
        % Check: if the exit condition is satisfied. If yes, then exit; continue otherwise.
        DataTermImprovement = 100*abs(Res1D_SpatRegVec(min3_nEdist(1))-Res1D_SpatRegVec(min3_nEdist(2)))/min(Res1D_SpatRegVec(min3_nEdist(1)),Res1D_SpatRegVec(min3_nEdist(2)));
        PriorTermImprovement = 100*abs(NormSpa_SpatRegVec(min3_nEdist(1))-NormSpa_SpatRegVec(min3_nEdist(2)))/min(NormSpa_SpatRegVec(min3_nEdist(1)), NormSpa_SpatRegVec(min3_nEdist(2)));    
        [NormSpa_SpatRegVec(min3_nEdist(1)), NormSpa_SpatRegVec(min3_nEdist(2)), PriorTermImprovement]
    %
        % Check: if AlphaS for two lowest distances are less than 100 apart. If yes, then exit; continue otherwise.
        boolResponse = 1;
        PriorTermImprHistory = nan(1,100);
        DataTermImprHistory = nan(1,100);

        iCount = 1;
        PriorTermImprHistory(1,iCount) = PriorTermImprovement;
        DataTermImprHistory(1,iCount)   = DataTermImprovement;
        W1 = 0.5; W2 = 0.5;
        debug_point(fileID);
        while (boolResponse)    
            minVec = [AlphaS(min3_nEdist(1)), AlphaS(min3_nEdist(2)), AlphaS(min3_nEdist(3))];
            [minVec, minVecIND] = sort(minVec);
            AlphaS_New1 = round((minVec(2)-minVec(1))*W1+minVec(1));        
            AlphaS_New2 = round((minVec(3)-minVec(2))*W2+minVec(2)); 
            AlphaS_NewVec = [AlphaS_New1, AlphaS_New2];
            %
            Response1 = any(AlphaS==AlphaS_New1) || any(AlphaS==AlphaS_New2);
            if (Response1 == 1)
                a = 0.25; b = 0.75;
                W1 = (b-a)*rand+a;   % a random number between 0.25, 0.75
                W2 = (b-a)*rand+a;   % a random number between 0.25, 0.75
                % Recalculate
                AlphaS_New1 = round((minVec(2)-minVec(1))*W1+minVec(1));        
                AlphaS_New2 = round((minVec(3)-minVec(2))*W2+minVec(2)); 
            else
                W1 = 0.5;
                W2 = 0.5;
            end
            %
            [[AlphaS(min3_nEdist(1)), AlphaS(min3_nEdist(2)), AlphaS(min3_nEdist(3))]; [0, AlphaS_New1, AlphaS_New2]]
            display(['From inside mus optimization loop: -------------------------------->', num2str(iCount)]);

            %Reconstruct for two new values of alphaS and return two vectors
            %[Res1D_New1, Res1D_New2] and [NormSpa_New1, NormSpa_New2]
            Capsule1.alphas_Info.alpha1             = AlphaS_NewVec;
            [ResVec2, NormSpaVec2] = QT2RProcess_wt_ParforMonitor_Stim3D(DataFolder, DataFileName, Size_DSW, WOR, Resol, MaskFileName, Spatial_ProcessMode, RegTag, MatPriorityTag, Capsule1, muS_Opt_Tag);

            AlphaS = [AlphaS, AlphaS_New1, AlphaS_New2];
            [AlphaS, IND] = sort(AlphaS);

            Res1D_New1 = ResVec2(1);        Res1D_New2 = ResVec2(2);        %Res1D_New3 = ResVec2(3);
            NormSpa_New1 = NormSpaVec2(1);  NormSpa_New2 = NormSpaVec2(2);  %NormSpa_New3 = NormSpaVec2(3);

            Res1D_SpatRegVec    = [Res1D_SpatRegVec, Res1D_New1, Res1D_New2];
            Res1D_SpatRegVec    = Res1D_SpatRegVec(IND);
            NormSpa_SpatRegVec  =[NormSpa_SpatRegVec, NormSpa_New1, NormSpa_New2];
            NormSpa_SpatRegVec  = NormSpa_SpatRegVec(IND);
            %
            %min3_nEdist = Find_Regularization_Par3(Res1D_SpatRegVec , NormSpa_SpatRegVec, PlotTag);
            min3_nEdist = Find_3Possible_LCorners(Res1D_SpatRegVec, NormSpa_SpatRegVec, AlphaS, PlotTag);
            AlphaS_Sep = (AlphaS(min3_nEdist(1))-AlphaS(min3_nEdist(2)));  

            DataTermImprovement = 100*abs(Res1D_SpatRegVec(min3_nEdist(1))-Res1D_SpatRegVec(min3_nEdist(2)))/min(Res1D_SpatRegVec(min3_nEdist(1)),Res1D_SpatRegVec(min3_nEdist(2)));
            PriorTermImprovement = 100*abs(NormSpa_SpatRegVec(min3_nEdist(1))-NormSpa_SpatRegVec(min3_nEdist(2)))/min(NormSpa_SpatRegVec(min3_nEdist(1)), NormSpa_SpatRegVec(min3_nEdist(2)));    
            %[NormSpa_SpatRegVec(min3_nEdist(1)), NormSpa_SpatRegVec(min3_nEdist(2)), PriorTermImprovement]
            iCount = iCount + 1;
            PriorTermImprHistory(1,iCount) = PriorTermImprovement;
            DataTermImprHistory(1,iCount)   = DataTermImprovement;
            %
            % Termination condition:
            % Based on "PriorTermImprovement": notice that PriorTermImprovement is in percentage. So, 0.1% is really small, typical change in PriorTermImprovement is few %
            if ((abs(PriorTermImprovement)>0.1) && (abs(AlphaS_Sep) > 50)) 
                boolResponse = 1;       % Continue; do not exit 
            elseif (abs(AlphaS_Sep) <= 50)   % NEW CONDITION -- April 2019
                boolResponse = 0;       % stop further iteration
            else                        % (abs(PriorTermImprovement)<0.1)
                boolResponse = 0;       % Exit
            end
            debug_point(fileID);
        end
        
        % To avoid oversmoothing, we restrict to indices of normDataTermFAE_Vec < (1+ErrorSR_T2Dist/100)*resFAE0 and use only a subset
    % % %     AlphaS              = AlphaS            (Res1D_SpatRegVec < (1+ErrorSR_T2Dist/100).*Res_0);
    % % %     NormSpa_SpatRegVec	= NormSpa_SpatRegVec(Res1D_SpatRegVec < (1+ErrorSR_T2Dist/100).*Res_0);
    % % %     Res1D_SpatRegVec    = Res1D_SpatRegVec  (Res1D_SpatRegVec < (1+ErrorSR_T2Dist/100).*Res_0); 

        %
        % Once outside the loop, the closest point on L-curve is the optimal pt
        % min_nEdist = Find_Regularization_Par(Res1D_SpatRegVec , NormSpa_SpatRegVec);
        % AlphaS_OPt = AlphaS(min_nEdist);
        min3_nEdist = Find_3Possible_LCorners(Res1D_SpatRegVec, NormSpa_SpatRegVec, AlphaS, PlotTag);
        min_nEdist = min3_nEdist(1);
    end
    AlphaS_OPt = AlphaS(min_nEdist); 
    % Saving the optimization report
    save([Convert_2_Forwardslash(DataFolder), 'Report_AlphaS_Optimization_DK.mat'], 'AlphaS', 'Res1D_SpatRegVec','NormSpa_SpatRegVec', 'AlphaS_OPt');
    d = debug_point(fileID);
    fprintf(fileID,'Done completely with file %s, function="%s", line %i\n', d.file, d.name, d.line);
    fprintf(fileID, 'time: %s', datestr(now,'mm-dd-yyyy HH-MM'));
end

