%% Load the mask
%%%%%%%%%%%%% 

function Mask3D = LoadAppropriateMask(MaskType, MaskFileNameWTFolder, dims, RotAng_Tag, FlipTag)
    Mask3D = zeros([dims]);

    if strcmp(MaskType, 'nii')
        [Mask3D, nii] = Load_nii_Im_DK(MaskFileNameWTFolder);
        for SliceNo = 1:dims(3)
            MaskTemp = squeeze(Mask3D(:,:,SliceNo));
            MaskTemp = double(rot90(MaskTemp, RotAng_Tag));
            
            %Flip the mask
            if strcmp(FlipTag, 'Yes')
                Mask3D(:,:,SliceNo) = flipdim(MaskTemp,1);   
            else
                Mask3D(:,:,SliceNo) = MaskTemp;
            end
        end
    elseif strcmp(MaskType, 'roi')
        %dims = [SizeVol(2), SizeVol(1), SizeVol(3)];
        %dims = [SizeVol(1), SizeVol(2), TotSlices];
        roi3d =readroi_DK(MaskFileNameWTFolder, dims');
        for SliceNo = 1:dims(3)
            MaskTemp = roi3d(:, :, SliceNo);
            MaskTemp = double(rot90(MaskTemp, RotAng_Tag));
            
            %Flip the mask
            if strcmp(FlipTag, 'Yes')
                Mask3D(:,:,SliceNo) = flipdim(MaskTemp,2);   
            else
                Mask3D(:,:,SliceNo) = MaskTemp;
            end
        end
    else
        Mask3D = ones([dims]);  
    end
    
end

% function Folder = Check_4_FileSep(Folder)
%     if ~strcmp(Folder(end), filesep)
%         Folder = [Folder, filesep];
%     end
% end

function [Im, nii] = Load_nii_Im_DK(FileName_wtDir)
    Unzipped_FileName_wtDir    =  Unzip_if_Zipped(FileName_wtDir);
    nii=load_untouch_nii(Unzipped_FileName_wtDir);
    Im = double(nii.img);
    %delete(Unzipped_FileName_wtDir);
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

function [roi3d] = readroi_DK(filename, dims_DK)
%function [roi3d]=readroi(filename, dims)
%
% READROI loads the roi file created by MRIcro and expands it to an array
% with size given by the dims vector
% 
% Input: 
%   filename - the name of the roi-file
%   dims     - 3 element vector which gives the array size
%              it must match the array size of the Analyze data where the
%              roi was originated.
%   
% Output:
%   roi3d - array with size as given in the dims vector
%           the pixel in POI have values one and all other pixels are zero.
%
% Jan Sedlacik. 
%
% A minor modification regarding dims done by DK

    %DK modified dimension only
    dims = [dims_DK(2), dims_DK(1), dims_DK(3)]';

    % create 3d array for final roi
    roi3d=zeros(dims(1),dims(2),dims(3),'int8');
    
    % open and read roi file from MRIicro
    fid=fopen(filename,'r');
    roi=fread(fid,'uint16');
    fclose(fid);

    % evaluate roi file and fill 3D roi array
    l=1; % set counter through roi file
    while l <= numel(roi)
        roi2d=zeros(dims(1),dims(2),'int8'); % define single slice
        slice=roi(l); % grep slice number
        if slice > 2^15, slice = slice -2^15;, end % test and corrects slice number if exceeded
        entries=roi(l+1); % grep no of entries of single slice
        for m=2:2:entries-2 % steps through entries of single slice
            index=roi(l+m); % grep start index of ROI in slice
            pixels=roi(l+m+1); % grep number of following ROI pixels

            % some tests and corrections if index exceeds numerical range of uni16
            ofset=0;
            negof=0;
            for n=1:12
                if pixels > 4096*n
                       ofset=2^(16)*n;
                       negof=4096*n;
                end
            end 

            % fill 2d roi slice with ROI pixels
            roi2d(index+ofset:index+ofset+pixels-1-negof)=1;
        end

        roi3d(:,:,slice)=roi2d; % fill 3d roi array with filled 2d roi slice
         l=l+entries; % increase roi-file counter for next slice
    end
end
