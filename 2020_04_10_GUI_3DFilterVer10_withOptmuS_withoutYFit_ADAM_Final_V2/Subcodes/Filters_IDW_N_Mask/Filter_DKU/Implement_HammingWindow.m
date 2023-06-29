
function Implement_HammingWindow()
%     close all
%     clear all
    
    % Read the image and perform fft2
    %InFolder = '/home/dkumar/Dropbox/PAPER-REVISION/T2W/Vol2';
    InFolder = 'C:\Users\kumar\Dropbox\PAPER-REVISION\T2W\Vol4';
    InFile = 'tse.nii.gz';  %  'tse.nii.gz', 'BET_Brain.nii.gz'
    %InFile = 'BET_Brain.nii.gz';
    
    
    cd('C:\Users\kumar\Dropbox\PAPER-REVISION\T2W\Vol4');
    Img = double(dicomread('DKU_HEALTHYVOLUNTEER_4.MR.TEST_DKU_FABIAN_SCAN.0014.0002.2014.12.16.19.00.04.890020.155965566.ima'));
    
% Read the nifti file
%     [Img, nii] = Load_nii_Im_DK([Check_4_Filesep(InFolder),  InFile]);
    Img_INV = fftshift(fft2(Img));

    %% Create 2D Hamming window through rotation formulation of 1D ver
    L = 16;
    w1D = hamming(L); % Some 1D window
    M = (L-1)/2;
    xx = linspace(-M,M,L);
    [x,y] = meshgrid(xx);
    r = sqrt( x.^2 + y.^2 );
    w2D = zeros(L);
    w2D(r<=M) = interp1(xx,w1D,r(r<=M));

    HammingWin = zeros(size(Img_INV));
    HammingWin(size(Img_INV,1)/2-M: size(Img_INV,1)/2+M, size(Img_INV,2)/2-M: size(Img_INV,2)/2+M) = w2D;

    h = HammingWin;

    figure(1);
    AX(1) = subplot(1,3,1); imagesc(Img); axis image; axis off;
    AX(2) = subplot(1,3,2); imagesc(h); axis image; axis off;

    filtered = ifft2(Img_INV.*h); 
    %AX(3) = subplot(1,3,3); imagesc(Img-abs(filtered)); axis image; axis off;
    AX(3) = subplot(1,3,3); imagesc(Img./abs(filtered)); axis image; axis off;    
    
    figure(40)
    imshow(Img./abs(filtered),[0 2])
    %imshow(rot90(Img./abs(filtered),1),[0 2])
end


function Folder = Check_4_Filesep(Folder)
    if ~strcmp(Folder(end), filesep)
        Folder = [Folder, filesep];
    end 
end

function [Im, nii] = Load_nii_Im_DK(FileName_wtDir)
    Unzipped_FileName_wtDir    =  Unzip_if_Zipped(FileName_wtDir);
    nii=load_untouch_nii(Unzipped_FileName_wtDir);
    Im = double(nii.img);
    delete(Unzipped_FileName_wtDir);
end

function Filename_with_Dir = Unzip_if_Zipped(Filename_with_Dir)
        %see if the file is in zipped format .gz
        %If it is unzip it
        SearchIndex = uint8(strfind(Filename_with_Dir, '.gz'));
        if ~isempty(SearchIndex)
            gunzip(Filename_with_Dir)
            Filename_with_Dir = strrep(Filename_with_Dir, '.gz', '');
        end
end