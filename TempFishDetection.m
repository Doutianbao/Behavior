

%% Fish detection

imgDir = 'Z:\SPIM\Avinash\RS neurons\Ablations\mCellArray-11-4-2015\Behavior\11-4-2015_mCellArray\vibStim\Fish2_ctrl\Trial_00_20151105_011410_AM';

IM = ReadImgSequence(imgDir);
IM_proc = ProcessImages(IM);
fishPos = GetFishPos(IM_proc);
motionFrames = GetMotionFrames(fishPos,1);

%%
figure('Name','Fish Outline')
frames = motionFrames(:)';
% frames = 1:size(IM_proc,3);
for jj = frames   
    cla
    imFish = GetFishPxls(IM_proc(:,:,jj),fishPos(jj,:),2,100);
    imFish(imFish~=0)=1;
    imFish = max(imFish(:))-imFish;
    imFish_dt = bwdist(imFish);
    blah = zeros(size(imFish_dt));
    blah(imFish_dt > 3 & imFish_dt <6)=1;
    imagesc(blah), axis image, colormap(gray),drawnow
    hold on
    shg
    title(['Fish Outline, Frame # ' num2str(jj)])
end
