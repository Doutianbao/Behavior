

%% Fish detection
% 
% imgDir = 'Z:\SPIM\Avinash\RS neurons\Ablations\mCellArray-11-4-2015\Behavior\11-4-2015_mCellArray\vibStim\Fish2_ctrl\Trial_00_20151105_011410_AM';
% % 
% IM = ReadImgSequence(imgDir);
% IM_proc = ProcessImages(IM, 110:210);
% fishPos = GetFishPos(IM_proc,50);
% motionFrames = GetMotionFrames(fishPos,2);

%%
% figure('Name','Fish Outline')
% % frames = motionFrames(:)';
% frames = 100:size(IM_proc,3);
% for jj = frames   
%     cla
%     imFish = GetFishPxls2(IM_proc(:,:,jj),fishPos(jj,:),[],80);
% %      imFish = GetFishPxls2(IM(:,:,jj),fishPos(jj,:),5,80);
%     imFish(imFish~=0)=1;
%     imFish = max(imFish(:))-imFish;
%     imFish_dt = bwdist(imFish);
%     blah = zeros(size(imFish_dt));
%     blah(imFish_dt > 3 & imFish_dt <6)=1;
%     imagesc(imFish), axis image, colormap(gray),drawnow
%     hold on
%     shg
%     title(['Fish Outline, Frame # ' num2str(jj)])
%     pause(0.15)
% end

% 
% %% Getting fish midlines for entire image stack
% figure
% % imgList = 1:size(IM,3);
% imgList = 80:300;
% for imgNum = imgList(:)'
% blah = IM(:,:,imgNum);
% blah  = max(blah(:))-blah;
% %blah(blah<intThr)=0;
% startPt = fishPos(imgNum,:);
% 
% lineLens = [18 16 14 10 8 8];
% lineInds = {};
% for jj = 1:length(lineLens)
%   if jj ==1
% %       lineInds{jj} = GetLineToEdge(blah,startPt,[],[],lineLens(jj));
%        lineInds{jj} = GetMidline(blah,startPt,[],[],lineLens(jj));
%   else
% %       lineInds{jj} = GetLineToEdge(blah,startPt,prevStartPt,[],lineLens(jj));
%       lineInds{jj} = GetMidline(blah,startPt,prevStartPt,[],lineLens(jj));
%   end
%     if jj ==3
%         a = 1;
%     end
%     si = lineInds{jj}(end);
%     [r,c] = ind2sub(size(blah),si);
%     x = c;
%     y = r;
%     prevStartPt = startPt;
%     startPt = [x,y];    
%     img = blah;
% end
% cla
% for kk = 1:length(lineInds)  
%     img(lineInds{kk}) = max(img(:)); 
%     imagesc(img),axis image        
% end
% hold on
% plot(fishPos(imgNum,1),fishPos(imgNum,2),'k*')
% title(num2str(imgNum))  
% shg
% % pause(0.2)
% end
% % lineInds = GetLineToEdge(blah,fishPos(560,:));
% 

%%
fishPos = tracexy;
midlineInds = GetMidline(IM,fishPos,18);


