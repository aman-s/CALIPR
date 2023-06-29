% Written by Dushyant Kumar
%
% It corrects a 3D map based on inverse distance weighting filter.
% It also use matlabpool and parfor loop to expedite the processing.
%
function Image3D_Interpolated = Interpolate_3D_Based_on_IDW_Fn_Parloop(Image3D)
    [x, y, z]  = meshgrid(1:size(Image3D,1), 1:size(Image3D,2), 1:size(Image3D,3));
    ImgCol  = Image3D(:); 
    warning('off', 'MATLAB:colon:nonIntegerIndex');
    Sz = size(Image3D);
    
    Image3D_Interpolated = zeros(Sz);
    Image3D_InterpolatedVec = zeros(Sz(1)*Sz(2)*Sz(3),1);
    parfor iVox = 1:Sz(3)*Sz(2)*Sz(1)
        [i,j,k] = ind2sub([Sz(1), Sz(2), Sz(3)], iVox);
        C = logical(((x == i) .* (y == j)) .*(z == k));
        [xInclude, yInclude, zInclude] = find(C == false);
        Ind = find(C == false);
        ImgInclude = ImgCol(Ind);
        [xExclude, yExclude, zExclude] = find(C == true);
        Image3D_InterpolatedVec(iVox, 1) =IDW_3D_DKU(xInclude, yInclude, zInclude, ImgInclude, xExclude, yExclude, zExclude, -2,'fr', 4);
    end
    Image3D_Interpolated = reshape(Image3D_InterpolatedVec, [Sz]);
    warning('on', 'MATLAB:colon:nonIntegerIndex');
end