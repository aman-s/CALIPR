
function TEST_INVERSE_DISTANCE()
    clc;
    %Basefolder = 'C:\Users\kumar\Dropbox\';
    Basefolder = '/home/dkumar/Dropbox/';

    cd([Basefolder, 'T2_Relaxometry_Codes_Aug2011/FAEProcess/Simulated_Data_SameFAE/FAs_90_180/Realistic_B1Map_MediumHigh@Corners/WithFAE/SNR_200_NoProfileError_repeat_test']);   
    S = load('Spatial_FAE_DK_SNR_200_Ex_n_Refoc_FA_Conven2_Slice1_alpha1_3000_Recur_1.mat');
    Img2Interpolate = S.FAErrorMap;

    figure(3);
    AX(1) = subplot(1,3,1); imagesc(Img2Interpolate); axis image;
    
    addpath('/home/dkumar/Dropbox/MatlabUserLib/IDW');

	[x, y]  = meshgrid(1:size(Img2Interpolate,1), 1:size(Img2Interpolate,2));
    xCol    = x(:);
    yCol    = y(:);
    ImgCol  = Img2Interpolate(:);

    ImgInterpolated = zeros(size(Img2Interpolate));
    for j = 1:size(Img2Interpolate,2)
        for i = 1: size(Img2Interpolate,1)
            C = logical((x == i) .* (y == j));
            [xInclude, yInclude] = find(C == false);
            Ind = find(C == false);
            ImgInclude = ImgCol(Ind);
            
            [xExclude, yExclude] = find(C == true);
             
            %Img2Interpolate(i,j) =IDW(xInclude, yInclude, ImgInclude, xExclude, yExclude,-2,'fr', 4);
            Img2Interpolate(i,j) =IDW(xInclude, yInclude, ImgInclude, xExclude, yExclude,-2,'ng',length(xInclude)/4);
        end
    end
    AX(2) = subplot(1,3,2); imagesc(Img2Interpolate); axis image;     
    
     ImgCol  = Img2Interpolate(:);

    ImgInterpolated = zeros(size(Img2Interpolate));
    
    for j = 1:size(Img2Interpolate,2)
        for i = 1: size(Img2Interpolate,1)
            C = logical((x == i) .* (y == j));
            [xInclude, yInclude] = find(C == false);
            Ind = find(C == false);
            ImgInclude = ImgCol(Ind);
            
            [xExclude, yExclude] = find(C == true);
            
            %ImgInterpolated(i,j) =IDW(xInclude, yInclude, ImgInclude, xExclude, yExclude,-2,'fr', 4);
            ImgInterpolated(i,j) =IDW(xInclude, yInclude, ImgInclude, xExclude, yExclude,-2,'ng',length(xInclude));
        end
    end   
    AX(3) = subplot(1,3,3); imagesc(ImgInterpolated); axis image; 
end