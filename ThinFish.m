headPos = [61, 61];
nLinePxls = 60;
maxFishPxls = 100;
maxDist = 20;
I_crop = CropImgsAroundPxl(IM,fishPos,70);
I_crop = max(IM(:))-I_crop;

%%
% img = I_crop(:,:,15);
% rImg = RadialFish(img,[61, 61],1,60);
S = @(pxl,pxls)(sum((pxls-repmat(pxl,size(pxls,1),1)).^2,2).^0.5);
I_thresh = size(I_crop);
linePxls = zeros(nLinePxls,2,size(I_crop,3));
for jj = 1:size(I_crop,3)
    cla
    try
        thresh = multithresh(I_crop(:,:,jj),3);
        
    catch
        try
            thresh = multithresh(I_crop(:,:,jj),2);
        catch
            thresh = multithresh(I_crop(:,:,jj),1);
        end
    end
    thresh = mean(thresh(1:min(2,length(thresh))));
    repeat = true;
    count =0;
    while repeat && count <=100
        img = imquantize(I_crop(:,:,jj),thresh);
        img_bw = bwmorph(Standardize(img),'thin','Inf');
        [y,x] = find(img_bw);
        count2 = 0;
        while (size(x,1) > maxFishPxls) && (count2 < 100)
            thresh = thresh*1.1;
            img = imquantize(I_crop(:,:,jj),thresh);
            img_bw = bwmorph(Standardize(img),'thin','Inf');
            [y,x] = find(img_bw);
        end
        seedPos = headPos;
        for kk = 1:size(x,1)
            dS = S(seedPos,[x,y]);
            [s, pxlInd] = min(dS);
            if s < maxDist
                pxl1 = [x(pxlInd),y(pxlInd)];
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
    imagesc(img),axis image
    title(num2str(jj))
    hold on
    plot(linePxls(:,1,jj),linePxls(:,2,jj),'g.','markersize',6)
    drawnow
%     pause(0.1)
end