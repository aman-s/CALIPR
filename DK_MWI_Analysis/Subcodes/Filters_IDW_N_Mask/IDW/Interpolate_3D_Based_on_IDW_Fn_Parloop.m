% Iterpolation of 3D image is performed on 3D basis
function Image3D_Interpolated = Interpolate_3D_Based_on_IDW_Fn_Parloop(Image3D, Res,Mask)
    [x, y, z]  = meshgrid(1:size(Image3D,1), 1:size(Image3D,2), 1:size(Image3D,3));
    x = permute(x,[2 1 3]);   
    y = permute(y,[2 1 3]);   
    z = permute(z,[2 1 3]);     
    Sz = size(Image3D);
    Image3D_InterpolatedVec = zeros(Sz(1)*Sz(2)*Sz(3),1);    
    parfor iVox = 1:Sz(3)*Sz(2)*Sz(1)
        [i,j,k] = ind2sub([Sz(1), Sz(2), Sz(3)], iVox);
        if (Mask(i,j,k) ~= 0)
            C = logical(((x == i).* (y == j)).*(z == k));
            indexInclude = find(C == false);% JanS: find returny only index and not x,y,z coordinates
            [xInclude, yInclude, zInclude] = ind2sub([Sz(1), Sz(2), Sz(3)], indexInclude); % JanS: convert find index to x,y,z
            Image3D_InterpolatedVec(iVox, 1) =IDW_3D_mitRes_DKU(Res, xInclude, yInclude, zInclude, Image3D(~C), i, j, k, -2,'fr', 4);
        end
    end    
    Image3D_Interpolated = reshape(Image3D_InterpolatedVec, Sz(1), Sz(2), Sz(3));
end

