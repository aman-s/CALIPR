
function TEST_INVERSE_DISTANCE_3D()
    clc;
    %Basefolder = '/home/dkumar/Dropbox/T2_Relaxometry_Codes_Aug2011/FAEProcess/Simulated_Data_SameFAE/FAs_90_180/Realistic_B1Map_MediumHigh@Corners/WithFAE/SNR_200_NoProfileError_repeat_test');

    Basefolder = 'C:\Users\kumar\Dropbox\';
    Basefolder = '/home/dkumar/Dropbox/';

    addpath(['/home/dkumar/Dropbox/MatlabUserLib/3dRotation']);
    
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

    ImgInterpolated3d = zeros(50, 50, 11);
    ImgInterpolated3d(:,:,6) = Img2Interpolate;
    
    

    I = Img2Interpolate;
    [x y] = meshgrid(1:50);
    
    Image3D = zeros([size(Img2Interpolate), 11]);
    Image3D(:,:,6) = Img2Interpolate;
    for k =1:5
        [xi yi] = meshgrid(1:(1+k*.1):50);
        TempImage = 0.26*(1+0.2*rand(50, 50));
        Image2 = cast(interp2(x,y,double(I),xi,yi,'linear'), 'double');
        Sz = size(Image2);
        TempImage(25-Sz(1)/2:25+Sz(1)/2-1, 25-Sz(2)/2:25+Sz(2)/2-1) = Image2;
        
        AX(2) = subplot(1,3,2); imagesc(TempImage); axis image; 
        Image3D(:,:,6+k) = TempImage;
        Image3D(:,:,6-k) = TempImage;
    end
    
    % Let us try to implement 3D inverse 
    Image3D_Interpolated = Interpolate_3D_Based_on_IDW_Fn(Image3D);
    
    figure(10)
    AX(1) = subplot(3,4,1); imagesc(Image3D_Interpolated(:,:,1)); axis off; axis image
    AX(1) = subplot(3,4,1); imagesc(Image3D_Interpolated(:,:,1)); axis off; axis image
    AX(2) = subplot(3,4,2); imagesc(Image3D_Interpolated(:,:,2)); axis off; axis image
    AX(3) = subplot(3,4,3); imagesc(Image3D_Interpolated(:,:,3)); axis off; axis image
    AX(4) = subplot(3,4,4); imagesc(Image3D_Interpolated(:,:,4)); axis off; axis image
    AX(5) = subplot(3,4,5); imagesc(Image3D_Interpolated(:,:,5)); axis off; axis image
    AX(6) = subplot(3,4,6); imagesc(Image3D_Interpolated(:,:,6)); axis off; axis image
    AX(7) = subplot(3,4,7); imagesc(Image3D_Interpolated(:,:,7)); axis off; axis image
    AX(8) = subplot(3,4,8); imagesc(Image3D_Interpolated(:,:,8)); axis off; axis image
    AX(9) = subplot(3,4,9); imagesc(Image3D_Interpolated(:,:,9)); axis off; axis image
    AX(10) = subplot(3,4,10); imagesc(Image3D_Interpolated(:,:,10)); axis off; axis image
    AX(11) = subplot(3,4,11); imagesc(Image3D_Interpolated(:,:,11)); axis off; axis image
end

function Image3D_Interpolated = Interpolate_3D_Based_on_IDW_Fn(Image3D)
    [x, y, z]  = meshgrid(1:size(Image3D,1), 1:size(Image3D,2), 1:size(Image3D,3));
    ImgCol  = Image3D(:); 
    warning('off', 'MATLAB:colon:nonIntegerIndex');
    Image3D_Interpolated = zeros(size(Image3D));
    
    for k = 1:size(Image3D,3)
        for j = 1:size(Image3D,2)
            for i = 1: size(Image3D,1)
                [i,j,k]
                C = logical(((x == i) .* (y == j)) .*(z == k));
                [xInclude, yInclude, zInclude] = find(C == false);
                Ind = find(C == false);
                ImgInclude = ImgCol(Ind);

                [xExclude, yExclude, zExclude] = find(C == true);

                Image3D_Interpolated(i,j, k) =IDW_3D_DKU(xInclude, yInclude, zInclude, ImgInclude, xExclude, yExclude, zExclude, -2,'fr', 4);
            end
        end
    end
end