
function filtered_Image  = Filter_KSpace_3D(Img)
    Img = double(Img); 
    [nx ny nz] = size(Img); 
    % Use Butterworth or Gaussian high pass filter. 
    filter = ones(2*nx-1,2*ny-1); 
    d0 = 25;     % 10
    n = 10;         % Only used for Butterworth 
    
    for i = 1:2*nx-1 
         for j =1:2*ny-1
             for k =1:2*nz-1
                 dist = ((i-(nx+1))^2 + (j-(ny+1))^2 + (k-(nz+1))^2 )^.5; 
                 % Use Butterworth high pass filter. 
                 %filter(i,j)= 1/(1 + (dist/d0)^(2*n)); 
                 %filter(i,j)= 1.0 - filter(i,j); 
                 % Use Gaussian high pass filter.  
                 filter(i,j,k) = exp(-dist^2/(2*d0^2)); 
                 %filter(i,j) = 1.0 - filter(i,j); 
             end
         end 
    end 
    
    %Find Fouriern transform
    % fftu = fft2(Img,2*nx-1,2*ny-1); 
    fftu = fftn(Img, [2*nx-1,2*ny-1,2*nz-1]);                 
     
    fftu = fftshift(fftu); % rearranges the outputs of fft, fft2, and fftn by moving the zero-frequency component to the center of the array                                                                                                                 
    % Update image with high frequencies. 
    filtered_KSpace = filter.*fftu; 
    
    
    %subplot(1,2,2);  mesh(log(1+abs(filtered_KSpace-fftu))); 
    filtered_KSpace = ifftshift(filtered_KSpace); 
    filtered_Image   = ifftn(filtered_KSpace, [2*nx-1,2*ny-1,2*nz-1]);                  
    filtered_Image  = real(filtered_Image(1:nx,1:ny,1:nz)); 
end