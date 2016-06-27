
cropRad = 70;
nLinePxls =55;
maxFishPxls = 70;
maxDist = 10;
nHoodSize = 3;
headPos = [cropRad+1, cropRad+1];
I_crop = CropImgsAroundPxl(IM,fishPos, cropRad);
I_crop = max(IM(:))-I_crop;
I_crop = Standardize(I_crop);

%%
S = @(pxl,pxls)(sum((pxls-repmat(pxl,size(pxls,1),1)).^2,2).^0.5);
I_thresh = size(I_crop);
linePxls = zeros(nLinePxls,2,size(I_crop,3));
% indMat = reshape(1:size(I_crop,1)*size(I_crop,2),size(I_crop,1),size(I_crop,2));
[C,R] = meshgrid(1:size(I_crop,1),1:size(I_crop,2));
for jj = 1:size(I_crop,3)
    cla
    %# Get threshold for binarizing image
    try
        thresh = multithresh(I_crop(:,:,jj),3);        
    catch
        try
            thresh = multithresh(I_crop(:,:,jj),2);
        catch
            thresh = multithresh(I_crop(:,:,jj),1);
        end
    end
    thresh = 0.75*mean(thresh(1:min(2,length(thresh))));
% thresh = 0.8;
    seedPos = headPos;
    linePxls_temp = [];
    %# Loop through each image
    repeat = true;
    count =0;
    while repeat && count <=100
        img = Standardize(imquantize(I_crop(:,:,jj),thresh));
        supInds = find(img);
        img_bw = bwmorph(img,'thin','Inf');
        %         cc = bwconncomp(img_bw);
        %         dS = S(headPos,[x,y]);
        [y,x] = find(img_bw);
        count3= 0;
        while (size(x,1) < nLinePxls) && (count3 <=100)
            thresh = thresh*0.9;
            img = imquantize(I_crop(:,:,jj),thresh);
            img_bw = bwmorph(Standardize(img),'thin','Inf');
            [y,x] = find(img_bw);
        end
        %# Iteratively increase threshold until no more than a certain # of fish pxls are doung
        %         count2 = 0;
        %         while (size(x,1) > maxFishPxls) && (count2 < 100)
        %             thresh = thresh*1.1;
        %             img = imquantize(I_crop(:,:,jj),thresh);
        %             img_bw = bwmorph(Standardize(img),'thin','Inf');
        %             [y,x] = find(img_bw);
        %         end
        %# Starting from the head pxl find the next closest pxl
        if isempty(linePxls_temp)
            startInd = 1;
        else
            startInd = find(sum(linePxls_temp==0,2)==2,1,'first');
        end
        for kk = startInd:nLinePxls
            dS = S(seedPos,[x,y]);
            zerInds = find(dS==0);
            x(zerInds)=[];
            y(zerInds)=[];
            dS(zerInds)=[];
            [s, pxlInd] = min(dS);
            if s < maxDist
                pxl1 = [x(pxlInd),y(pxlInd)];
                
                %# Weighted sum of neighboring pxls to get a smooth line
                rInds = pxl1(1)-nHoodSize:pxl1(1)+nHoodSize;
                rInds(rInds<=0 | rInds > size(img,1))=[];
                cInds = pxl1(2)-nHoodSize:pxl1(2)+nHoodSize;
                cInds(cInds<=0 | cInds>size(img,2))=[];
                rNeighbors = R(rInds, cInds);
                cNeighbors = C(rInds, cInds);
                
                wts =  SparseIndex(I_crop(:,:,jj),rNeighbors(:),cNeighbors(:));%
                pxl1 =  ([rNeighbors(:),cNeighbors(:)]'*wts(:))'/sum(wts);
                pxl1 = round(pxl1*10)/10;
                seedPos =pxl1;
                x(pxlInd)=[]; y(pxlInd)=[];
                if kk <= nLinePxls
                    linePxls(kk,:,jj) = pxl1;
                end
            end
        end
        linePxls_temp = squeeze(linePxls(:,:,jj));
        if all(sum(linePxls_temp==0,1))
            repeat = true;
            thresh = 0.9*thresh;
            count = count + 1;
        else
            repeat = false;
        end
    end
    imagesc(img),axis image, colormap(gray)
%     imagesc(I_crop(:,:,jj)),axis image, colormap(gray)
%     cla
%     imagesc(img),colormap(gray)
%     cla
%     imagesc(img_bw),colormap(gray)
    title(num2str(jj))
    hold on
    plot(linePxls(:,1,jj),linePxls(:,2,jj),'r.','markersize',6)
    plot(cropRad+1,cropRad+1,'c+','markersize',10)
    shg
    drawnow
    if jj >= 52         
           pause()
    end
end