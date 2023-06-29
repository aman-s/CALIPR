function Test_Mask_Orient_WRT_Img(T2Data4D, TE, Mask3D, dims, SliceNo, ClickTag)    
    ImageSlice = squeeze(T2Data4D(:,:,SliceNo,1));
    MaskSlice  = squeeze(Mask3D(:,:,SliceNo));
    figure(100);
    colormap('jet')
    AX1     = subplot(2,2,1); imagesc(ImageSlice); axis image; axis off;
    MaskSlice   = bwperim(MaskSlice);
    AX2     = subplot(2,2,2); imagesc(MaskSlice ); axis image; axis off;
    
    %% Overlay mask-rim (black color)  on the top of image
    MaxVoxX = dims(1); MaxVoxY = dims(2);
    dtmp = ImageSlice(1:MaxVoxX, 1:MaxVoxY);
    dtmp=uint8(double(dtmp)./double(max(dtmp(:))).*255);
    rgbtmp(:,:,1)=dtmp+1000.*(uint8(MaskSlice(:,:)));
    AX3 = subplot(2,2,3); im1 = imshow(rgbtmp); axis image; axis off; colormap('jet') 
     
    if strcmp(ClickTag, 'Yes')
        for k=1:50
            disp('Button click')
            [y,x] = ginput(1);
            x = round(x); y = round(y);
            T2DataVox = squeeze(T2Data4D(x,y,SliceNo,:));
            AX3 = subplot(2,2,4); plot(TE, T2DataVox);
        end 
    end
    %close;
end

