

clear; 
    cd('C:\Users\kumar\Dropbox\T2_Relaxometry_Codes_Aug2011\FAEProcess\Simulated_Data_SameFAE\Realistic_B1Map\SNR_100_NoProfileError')
%     S = load('FAEMap0rg.mat');
%     Img = S.FAEMapOrg;
    
    S = load('FAErrorMap1_NotAvged.mat');
    Img = S.FAE;
    Img = double(Img); 
    [nx ny] = size(Img); 

    u = Img; 
    Img = uint8(u); 
    imwrite(Img, 'micro.jpg'); 
    fftu = fft2(u,2*nx-1,2*ny-1); 
    fftu = fftshift(fftu); 
    subplot(1,2,1); mesh(log(1+(abs(fftu)))); 
    
    
    % Use Butterworth or Gaussian high pass filter. 
    filter = ones(2*nx-1,2*ny-1); 
    d0 = 12; 
    n = 10;         % Only used for Butterworth 
    for i = 1:2*nx-1 
         for j =1:2*ny-1 
             dist = ((i-(nx+1))^2 + (j-(ny+1))^2)^.5; 
             % Use Butterworth high pass filter. 
             %filter(i,j)= 1/(1 + (dist/d0)^(2*n)); 
             %filter(i,j)= 1.0 - filter(i,j); 
             % Use Gaussian high pass filter. 
             filter(i,j) = exp(-dist^2/(2*d0^2)); 
             %filter(i,j) = 1.0 - filter(i,j); 
         end 
    end 

    % Update image with high frequencies. 
    fil_micro = filter.*fftu; 
    %subplot(1,2,2);  mesh(log(1+abs(fil_micro-fftu))); 
    fil_micro = ifftshift(fil_micro); 
    fil_micro = ifft2(fil_micro,2*nx-1,2*ny-1); 
    fil_micro = real(fil_micro(1:nx,1:ny)); 
    fil_micro = uint8(fil_micro); 
    imwrite(fil_micro, 'micro_fil.jpg');
    
    figure(100);
    AX(1) = subplot(1,2,1); imagesc(Img); axis image; axis off; caxis([-10 10]);
    AX(1) = subplot(1,2,2); imagesc(fil_micro); axis image; axis off; caxis([-10 10]);