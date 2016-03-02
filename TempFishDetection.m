

%% Fish detection
% 
% imgDir = 'Z:\SPIM\Avinash\RS neurons\Ablations\mCellArray-11-4-2015\Behavior\11-4-2015_mCellArray\vibStim\Fish2_ctrl\Trial_00_20151105_011410_AM';
% 
% IM = ReadImgSequence(imgDir);
% IM_proc = ProcessImages(IM, 110:210);
% fishPos = GetFishPos(IM_proc);
% motionFrames = GetMotionFrames(fishPos,2);

%%
figure('Name','Fish Outline')
% frames = motionFrames(:)';
frames = 100:size(IM_proc,3);
for jj = frames   
    cla
    imFish = GetFishPxls2(IM_proc(:,:,jj),fishPos(jj,:),[],80);
%      imFish = GetFishPxls2(IM(:,:,jj),fishPos(jj,:),5,80);
    imFish(imFish~=0)=1;
    imFish = max(imFish(:))-imFish;
    imFish_dt = bwdist(imFish);
    blah = zeros(size(imFish_dt));
    blah(imFish_dt > 3 & imFish_dt <6)=1;
    imagesc(imFish), axis image, colormap(gray),drawnow
    hold on
    shg
    title(['Fish Outline, Frame # ' num2str(jj)])
    pause(0.15)
end
