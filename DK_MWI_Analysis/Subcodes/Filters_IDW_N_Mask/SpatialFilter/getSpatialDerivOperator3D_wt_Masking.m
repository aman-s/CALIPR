% This operator is only used for FAE smoothing and for the smoothing
% of T2 distributions

function  Dstr = getSpatialDerivOperator3D_wt_Masking(Filter, sz, Mask_ROI, Resol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: 3-D size of T2_Scale spectrum - nrows x ncols x nT2
% returns a matrix of size (nrows.ncols.T2_Scale x nrows.ncols.T2_Scale) whose
% operation on a stacked vector of size (nrows.ncols.T2_Scale x 1) 
%
% will give the same result as if you convolved the underlying image by the
% highpass filter. The highpassfilter should be a 2D filter only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    FilterSize  = Filter.FilterSize;
    FilterDim   = Filter.FilterDim;   % '2D', '3D'
    FilterTag   = Filter.FilterTag;   % 'IDW', 'ExpInv', 'Distance', 'Distance-Sq', 'Equal' 
    mm = sz(1);
    nn = sz(2);
    %Calculate spatial term
    if strcmp(FilterDim, '2D') 
        hipassfilt = Create_1stDiffOp([FilterSize(1), FilterSize(2), 1], [Resol(1), Resol(2)], FilterTag);
        hp_Org = hipassfilt;
        if (length(sz) ~= 3)
            disp('The dimension of parameter is not correct. It should be 3D array')
            stop
        else
            lt2 = sz(3);
            pp = 1;
            hp_Org(abs(hp_Org) < 0.02*max(abs(hp_Org(:)))) = 0;
            del = floor(size(hp_Org)/2);   
        end
    elseif strcmp(FilterDim, '3D')
        hipassfilt = Create_1stDiffOp([FilterSize(1), FilterSize(2), FilterSize(3)], Resol, FilterTag);
        hp_Org = hipassfilt;  
        if (length(sz) ~= 4)
            disp('The dimension of parameter is not correct. It should be 4D array');
            stop
        else
            pp = sz(3);
            lt2 = sz(4);
            hp_Org(abs(hp_Org) < 0.02*max(abs(hp_Org(:)))) = 0;
            del = floor(size(hp_Org)/2);            
        end
    else
        disp('High pass filter type not defined');
        stop
    end
    %%%   
    Dstr = sparse(mm*nn*pp, mm*nn*pp);
    for ind = 1:mm*nn*pp
        if strcmp(FilterDim, '2D')
            [i,j] = ind2sub([mm,nn], ind);  % pp = 1 for FilterDim--> '2D'
            D0 = zeros(mm+2*del(1), nn+2*del(2));
            %D0 = sparse(mm+2*del(1), nn+2*del(2));
            
            %Vertex1 = [max(i-1, 1), max(j-1, 1)]; 
            %Vertex2 = [min(i+1, size(Mask_ROI,1)), min((j+1), size(Mask_ROI,2))];
               
            if (Mask_ROI(i,j) == 1) %VOI is in WM
                Mask0 = Mask_ROI;
            else %VOI is in other tissue type
                Mask0 = 1-Mask_ROI;
            end
            
            D0(1:size(hp_Org,1), 1:size(hp_Org,2)) = hp_Org; 	% hp is at upper left corner.    
            
            D = circshift(D0, [i-1,j-1]);
            D1 = D(del(1)+1:del(1)+mm, del(2)+1:del(2)+nn);
            D1 = D1 .*Mask0;                                    % MODIFIED to account for TissueClassMask
            
            %Convert 2D to 1D properly
            d1 = reshape(D1,[size(D1, 2)*size(D1, 1), 1]);

        elseif strcmp(FilterDim, '3D')
            [i,j,k] = ind2sub([mm,nn, pp], ind);
            D0 = zeros(mm+2*del(1), nn+2*del(2), pp+2*del(3));
            %D0 = sparse(mm+2*del(1), nn+2*del(2), pp+2*del(3));  

            if (Mask_ROI(i,j,k) == 1) %VOI is in WM
                Mask0 = Mask_ROI;
            else %VOI is in other tissue type
                Mask0 = 1-Mask_ROI;
            end
                       
            D0(1:size(hp_Org,1), 1:size(hp_Org,2), 1:size(hp_Org,3))= hp_Org;               % hp is at upper left corner.          
            D = circshift(D0, [i-1,j-1,k-1]);
            D1 = D(del(1)+1:del(1)+mm, del(2)+1:del(2)+nn, del(3)+1:del(3)+pp);
            
            D1 = D1 .*Mask0;   %% MODIFIED to account for TissueClassMask    
            % Convert 3D to 1D: Spatial only
            d1 = ConvertArr_3D_to_1D_SpatialOnly(D1);
            
        else
            disp('High pass filter not defined')
            exit
        end 

        % The next block ensures that sum of coeff of high pass filter is zero as voxels @ the corner does not get the filter right
        if (abs(sum(d1))>= 1e-6)  % Ideally we would like to check if ~(sum(d1)==0)
            %Find position of of nonzero entires ~= -0.5 in that column  
            Sum0 = sum(d1(find(d1 ~= -0.5)))+eps; 
            d1(find((d1 == -0.5))) = -Sum0; 
        end
        
        % Dstr(ind, :) = d1; %    >---------------CHECK  trying OLD 
        Dstr(:, ind) = d1; %    >---------------CHECK  trying new
    end
end

function Weight = Create_1stDiffOp(FilterSize, Resol, FilterTag)
    theta0 = 0.1;   % NEED TO TWIQ EXP FUNCTION FURTHER
    %p = 2;         % Originally weighting was disance squared --till NEUROIMG SUBMISSION    
    p = 1;          % weighting distance .. correlation goes much far   
    p1 = p/2;    
    %
    R = Resol;          R_Min = min(Resol);         R_MinSq = R_Min^2;  
    Weight = zeros([FilterSize]);
    Center = ceil(FilterSize/2);
    Ic = Center(1); Jc = Center(2);    
    if ((length(FilterSize) == 3) && (FilterSize(3) ~= 1)) % 3D filter
        if strcmpi(FilterTag, 'Equal-Local')
            Weight = ones(3, 3, 3)./52;
            Weight(2, 2, 2) = -26/52;                        
        else
            Kc = Center(3);
            for k = 1:FilterSize(3)
                for j = 1:FilterSize(2)
                    for i = 1:FilterSize(1)
                        Dist_Sq = (R(1)*(i-Ic))^2 + (R(2)*(j-Jc))^2 + (R(3)*(k-Kc))^2;
                        if (Dist_Sq ~= 0.0)
                            if strcmpi(FilterTag, 'IDW-NonLoc')
                                Dist_pow_p = Dist_Sq^p1/(R_MinSq^p1);
                                Weight(i,j,k) = 1/Dist_pow_p; 
                            else %strcmpi(FilterTag, 'ExpDecay-NonLoc')
                                Weight(i,j,k) = exp(-theta0*((Dist_Sq/R_MinSq)^0.5));   
                            end
                        end
                    end
                end
            end
            SumW = sum(Weight(:));
            Weight(Center(1), Center(2), Center(3)) = -SumW;
            Weight = Weight./(2*SumW);
        end
    else
        if strcmpi(FilterTag, 'Equal-Local')
            Weight = [1 1 1; 1 -8 1; 1 1 1]/16;         
        else
            for j = 1:FilterSize(2)
                for i = 1:FilterSize(1)
                    Dist_Sq = (R(1)*(i-Ic))^2 + (R(2)*(j-Jc))^2;
                    if (Dist_Sq ~= 0.0)
                        if strcmpi(FilterTag, 'IDW-NonLoc') 
                            Dist_pow_p = Dist_Sq^p1/(R_MinSq^p1);
                            Weight(i,j) = 1/Dist_pow_p;
                        else %strcmpi(FilterTag, 'ExpDecay')
                            Weight(i,j) = exp(-theta0*((Dist_Sq/R_MinSq)^0.5));    
                        end
                    end
                end
            end
            SumW = sum(Weight(:));
            Weight(Center(1), Center(2)) =  -SumW;
            Weight = Weight./(2*SumW); 
        end        
    end

end
