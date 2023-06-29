% Written by Dushyant Kumar
%
% It corrects a 3D map based on inverse distance weighting filter.
% It DOES NOT use matlabpool.
%
% Use the other code (Interpolate_3D_Based_on_IDW_Fn_Parloop,m) if you want
% to expedite the processing using matlabpool
%
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
    warning('on', 'MATLAB:colon:nonIntegerIndex');
end