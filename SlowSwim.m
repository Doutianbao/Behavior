
% switch  'LoadingNewFilm'  %'LoadingNewFilm' 'RerunAnalysis' 'LoadingCoordinates'
%     case 'LoadingNewFilm'
%

clear all
close all

cd 'S:\Avinash\Tracking demo'

readMode =  'fromImages';
% readMode = 'fromMishVid';

poolSize  = 10;
switch readMode
    case 'fromMishVid'
        fName_prefix = input('Enter fish name, e.g., Fish1: ','s');
        [IM, outDir] = ReadMishVid();
        imgInds = 1:size(IM,3);
    case 'fromImages'
        imgDir = input('Enter image dir path:  ', 's')
        imgExt = input('Enter image extension, e.g. jpg:  ','s')
        imgInds = input('Enter indices of images to read as a vector:  ');
        fName_prefix = input('Enter fish name, e.g., Fish1: ','s');
        bpMessg =['Filter bandpass for finding fish pos, e.g. [15 25]. ' ...
            'Press enter [] to skip (filtering is recommended if using collimated ' ...
            'light during behavior): '];
        fps = input('Enter frame rate (frames/sec): ');
        bp = input(bpMessg);
        outDir = fullfile(imgDir,'proc');
        IM = ReadImgSequence(imgDir,imgExt,imgInds);
end

if exist(outDir)~=7
    mkdir(outDir)
end


%% Processing images
tic
if matlabpool('size')==0
    matlabpool(poolSize)
end
disp('Processing images...')
[IM_proc, ref] = SubtractBackground(IM);
imgDims = size(IM_proc);
toc

%% Tracking the fish
tic
disp('Getting fish pos...')
if ~isempty(bp)
%     fishPos = GetFishPos(IM_proc, 30,'filter',bp,'process','parallel');
        fishPos = GetFishPos(IM_proc, 25,'filter',bp,'process','serial');
else
    fishPos = GetFishPos(IM_proc, 25,'process','parallel');
    %     fishPos = GetFishPos(IM_proc, 30,'process','serial');
end
toc

cropWid = 70;
disp('Cropping images...')
IM_proc_crop = CropImgsAroundPxl(IM_proc,fishPos,cropWid);
clear IM_proc
toc
disp('Saving fish position and ref image...')
if isempty(imgInds)
    imgInds = 1:size(IM,3);
end
fName_suffix = [num2str(round(imgInds(1)/60/30)) '-' num2str(round(imgInds(end)/60/30)) 'mins'];
fName = strcat(fName_prefix,'_',fName_suffix);
ts = datestr(now,30);
procData = matfile(fullfile(outDir,['procData_' ts '.mat']),'Writable',true);
procData.fishPos = fishPos;
procData.ref = ref;
toc

%% Fish Orientation
disp('Getting tail curvature...')
tic

% midlineInds = GetMidlines(IM_proc,fishPos,[20 15 15],'bmp','ref',ref);

%##### Not using this method at the moment.
% disp('Getting arena edge...')
% edgeInds  = GetArenaEdge(ref,'detThr',0.6,'nIter',100);
% disp('Getting extra arena pxls...')
% [~,outPxls] = GetInnerPxls(ref,fliplr(edgeInds));
% extraArenaInds = sub2ind(size(ref),outPxls(:,1),outPxls(:,2));
% midlineInds = GetMidlines(IM_proc,fishPos,[24 20 15],'bmp','extraArenaInds',extraArenaInds,'procType','parallel');
%###########

midlineInds = GetMidlines(IM_proc_crop,(fishPos./fishPos)*(size(IM_proc_crop,1)/2+1),...
    [15 10 10 10 10 5],'bmp','procType','parallel');
toc

tailCurv = SmoothenMidlines(midlineInds,IM_proc_crop,3,'plotBool',0,'pauseDur',0,'smoothFactor',5);

orientation = GetFishOrientationFromMidlineInds(midlineInds,imgDims(1:2),'s');
% orientation_backup = orientation;

disp('Saving midline inds, and tailCurv...')
% procData.orientation = orientation;
procData.midlineInds = midlineInds;
procData.tailCurv = tailCurv;
toc

%% Motion Info
motionThr = 1;
% [motionFrames, swimStartFrames] = GetMotionFrames(fishPos,motionThr);
motionInfo = GetMotionInfo(fishPos,orientation,imgDims(1),'motionThr',motionThr);

%% Saving processed images
saveOrNot = 'y';
% saveOrNot = input('Save cropped image stacks (y/n)?  ','s');
cropWid = input('Enter crop width in pxls: ');
if isempty(cropWid)
    cropWid = 70;
end
tic
if strcmpi('y',saveOrNot)   
    disp('Cropping images...')
    tic
    IM_crop = CropImgsAroundPxl(IM,fishPos,cropWid);
    clear IM   
    toc    
    disp('Saving cropped stacks...')
    tic
    procData.IM = IM_crop;
    toc
    procData.IM_proc = IM_proc_crop; % This is to save some time and space
    toc
else
    disp('Data not saved!')
end
toc

%% Getting and saving peak info
saveOrNot = input('Detect and save peak info (y/n)?  ','s');
tic
if strcmpi('y',saveOrNot)
    disp('Getting peak info...')
    GetTapDark1st2ndLeftRightPks;
    procData.pks = pks;
else
    disp('Pk info not saved!')
end
toc

break;
if matlabpool('size')>0
    matlabpool close
end

break

%% Plot max-int projected fish images for tap and dark flash trls
periPxls = 100;
nFramesInTrl = 750;
nTrls = size(IM_proc,3)/nFramesInTrl;
trls.tap = 1:2:nTrls;
trls.dark = 2:2:nTrls;
maxImgStacks = cell(2,1);
fldNames = fieldnames(trls);
for stim = 1:length(fldNames)
    N = length(trls.(fldNames{stim}));
    maxImgStacks{stim} = zeros(2*periPxls+1,2*periPxls+1,N);
    fh = figure('Name',fldNames{stim})
    nRows = 3;
    nCols = ceil(N/3);
    disp(['Stim: ' fldNames{stim}])
    count = 0;
    for trl = trls.(fldNames{stim})
        count = count + 1;
        disp(['Trl: ' num2str(trl)])
        inds = (nFramesInTrl*(trl-1)+1:nFramesInTrl*trl);
        IM_fix = PlayFixedFish(IM_proc(:,:,inds),fishPos(inds,:),orientation(1,inds)+180, ...
            'periPxls',periPxls,'dispBool',0);
        maxImg = max(IM_fix,[],3);
        maxImgStacks{stim}(:,:,count)= maxImg;
        figure(fh)
        subaxis(nRows,nCols,count,'SpacingHoriz',0, 'SpacingVert',0.05)
        imagesc(imNormalize999(maxImg)), axis image, colormap(gray), axis off
        set(gca,'clim',[0 0.9])
        title(num2str(trl))
        drawnow
        
    end
end
saveOrNot = input('Save max img stacks (y/n)?:', 's');
if strcmpi(saveOrNot,'y')
    save(fullfile(outDir,'Max img stacks for different trls.m'),'maxImgStacks');
end


%% Turn angles during swims and histogram
x = fishPos(:,1);
y = fishPos(:,2);
dx = diff(x);
dy = diff(y);
dS = sqrt(dx.^2 + dy.^2);
noMovInds  = find(dS < 10);
movInds = setdiff(1:length(dx),noMovInds);
x_mov = x(movInds);
y_mov = y(movInds);
dx_mov = zeros(size(movInds));
dy_mov = zeros(size(movInds));

orientation_mov = fish.orientation(movInds);
dO_mov = diff(orientation_mov);
dO_mov(dO_mov <-180) = dO_mov(dO_mov < -180) + 360;
dO_mov(dO_mov >180) = dO_mov(fish.dO_mov > 180) - 360;

figure('Name','Turn angle histogram')
hist(fish.dO_mov,75)
box off
set(gca,'tickdir','out','xtick',-180:45:180)
xlabel('Turn angle')
ylabel('Count')
title('Turn angle histogram')


%% Swim distances and histogram
dMovInds = diff(movInds);
movInds_jump = movInds;
movInds_jump(dMovInds==1)=[];
x_mov = x(movInds_jump);
y_mov = y(movInds_jump);
dS_mov = sqrt(x_mov.^2 + y_mov.^2);

figure('Name','Swim distance histogram')
hist(dS_mov,50)
box off
xlabel('Swim distances (px)')
title('Swim dist histogram')
ylabel('Count')
set(gca,'tickdir','out')
title('Swim distance histogram')
shg

break

