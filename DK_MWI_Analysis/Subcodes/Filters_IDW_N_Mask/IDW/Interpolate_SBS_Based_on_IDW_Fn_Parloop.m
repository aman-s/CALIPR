% Iterpolation of 3D image is performed on slice by slice 
function Image3D_Interpolated = Interpolate_SBS_Based_on_IDW_Fn_Parloop(Image3D, Res)
    [x, y]  = meshgrid(1:size(Image3D,1), 1:size(Image3D,2));
    x = permute(x,[2 1]);   
    y = permute(y,[2 1]);       
    Sz = size(Image3D);
    
    if (length(size(Sz))==2) % Single slice
       Sz = [Sz, 1]; % Just for consistency with other 3D processing, added third dimension as 1.
    end
    
    Image3D_Interpolated = Image3D .* 0.0;
    for k = 1:Sz(3)
        Image2D = squeeze(Image3D(:,:,k));
        Image2D_InterpolatedVec = zeros(Sz(1)*Sz(2),1); 
        parfor iVox = 1:Sz(2)*Sz(1)
            [i,j] = ind2sub([Sz(1), Sz(2)], iVox);
            C = logical(((x == i).* (y == j)));
            indexInclude = find(C == false);% JanS: find returny only index and not x,y,z coordinates
            [xInclude, yInclude] = ind2sub([Sz(1), Sz(2)], indexInclude); % JanS: convert find index to x,y,z
            Image2D_InterpolatedVec(iVox, 1) =IDW_2D_mitRes_DKU(Res, xInclude, yInclude, Image2D(~C), i, j, -2,'fr', 4);
        end
        Image3D_Interpolated(:,:,k) = reshape(Image2D_InterpolatedVec, Sz(1), Sz(2));
    end
end
