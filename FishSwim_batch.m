function FishSwim_batch(imgDir,fishName)


readMode =  'fromImages';
% readMode = 'fromMishVid';

poolSize  = 10;
switch readMode
    case 'fromMishVid'
        fName_prefix = input('Enter fish name, e.g., Fish1: ','s');
        [IM, outDir] = ReadMishVid();
        imgInds = 1:size(IM,3);
    case 'fromImages'
        %         imgDir = input('Enter fish dir path:  ', 's')
        %         imgExt = input('Enter image extension, e.g. jpg:  ','s')
        imgExt = 'bmp';
        %         imgInds = input('Enter indices of images to read as a vector:  ');
        imgInds = [];
        %         fName_prefix = input('Enter fish name, e.g., Fish1: ','s');
        fName_prefix  = fishName;
        %         bpMessg =['Filter bandpass for finding fish pos, e.g. [15 25]. ' ...
        %             'Press enter [] to skip (filtering is recommended if using collimated ' ...
        %             'light during behavior): '];
        %         fps = input('Enter frame rate (frames/sec): ');
        fps = 500;
        %         nFramesInTrl  = input('# of frames in each trl (default = 750): ');
        nFramesInTrl = 750;
        %         bp = input(bpMessg);
        bp = 30;
        outDir = fullfile(imgDir,'proc');
        IM = ReadImgSequence(imgDir,imgExt,imgInds);
end

if exist(outDir)~=7
    mkdir(outDir)
end


%% Processing images block-by-block
blockSize = 4;
cropWid = 90; %( For imgDims ~ [900,900], use 70)
lineLen = 15;

tic
if matlabpool('size')==0
    matlabpool(poolSize)
end

disp('Block-by-block processing...')
blockInds = round((1:blockSize)*(size(IM,3)/blockSize));
blockInds = [1,blockInds(:)'];
blockInds(end) = size(IM,3);
imgDims = size(IM);
IM_proc_crop = zeros(2*cropWid + 1, 2*cropWid + 1, size(IM,3));
fishPos = zeros(size(IM,3),2);
hOr = cell(size(IM,3),1);
ref = zeros(imgDims(1),imgDims(2),blockSize);
im_proc = [];
for block = 1:numel(blockInds)-1
    imgFrames = blockInds(block):blockInds(block+1);
    disp(['Block # ' num2str(block),...
        ', images: ' num2str(imgFrames(1)) ' - ' num2str(imgFrames(end))])
    disp('Subtracting background...')
    [im_proc,ref(:,:,block)] = SubtractBackground(IM(:,:,imgFrames));
    if ~isempty(bp)
        [fp,hOr_temp] = GetFishPos(im_proc, 25,'filter',bp,'process','parallel',...
            'lineLen',lineLen);
    else
        [fp,hOr_temp] = GetFishPos(im_proc, 25,'process','parallel','lineLen',lineLen);
    end
    fishPos(imgFrames,:) = fp;
    hOr(imgFrames) = hOr_temp;
    disp('Cropping images...')
    IM_proc_crop(:,:,imgFrames) = CropImgsAroundPxl(im_proc,fp,cropWid);
end
ref = mean(ref,3);
clear im_proc fp hOr_temp
toc

%% Background subtraction
% [IM_proc, ref] = SubtractBackground(IM);
% imgDims = size(IM_proc);


%% Tracking the fish
% tic
% disp('Getting fish pos...')
% if ~isempty(bp)
%     %     fishPos = GetFishPos(IM_proc, 30,'filter',bp,'process','parallel');
%     [fishPos,hOr] = GetFishPos(IM_proc, 25,'filter',bp,'process','parallel','lineLen',15);
% else
%     [fishPos,hOr] = GetFishPos(IM_proc, 25,'process','parallel','lineLen',15);
%     %     fishPos = GetFishPos(IM_proc, 30,'process','serial');
% end
% toc

%% Cropping images, adjusting head orientation vector for cropped images, and saving
tic

disp('Adjusting head orientation vector for cropped images...')
hOr_crop = hOr;
imgDims = size(IM);
imgDims_crop = size(IM_proc_crop);
dispChunk = round(size(IM_proc_crop,3)/10);
dim = size(hOr{1},2);
try
    for jj = 1:length(hOr)
        if dim ==2
            x = hOr{jj}(:,1);
            y = hOr{jj}(:,2);
        else
            [y,x] = ind2sub(imgDims(1:2),hOr{jj});
        end
        blah = [x,y] - repmat(fishPos(jj,:),length(x),1) + repmat([cropWid+1, cropWid+1],length(x),1);
        x = blah(:,1); y = blah(:,2);
        hOr_crop{jj} = [x,y];
        if mod(jj,dispChunk)==0
            disp(num2str(jj))
        end
    end
catch
    disp('Getting head orientation in cropped stack...')
    tic
    [~,hOr_crop] =  GetFishPos(IM_proc_crop, 25,'filter',bp,'process','parallel','lineLen',15);
    toc
end

clear IM_proc
toc
disp('Saving fish position, ref image, and cropped images....')
if isempty(imgInds)
    imgInds = 1:size(IM,3);
end
fName_suffix = [num2str(round(imgInds(1)/60/30)) '-' num2str(round(imgInds(end)/60/30)) 'mins'];
fName = strcat(fName_prefix,'_',fName_suffix);
ts = datestr(now,30);
procData = matfile(fullfile(outDir,['procData_' ts '.mat']),'Writable',true);
procData.fishPos = fishPos;
procData.ref = ref;
procData.fps = fps;
procData.imgDims = size(IM);
procData.IM_proc_crop = IM_proc_crop;
procData.imgDims_crop = size(IM_proc_crop);
toc

%% Fish Orientation
disp('Getting tail curvature...')
tic
fp = repmat(ceil([size(IM_proc_crop,1), size(IM_proc_crop,2)]/2),size(fishPos,1),1);
[midlineInds,dsVecs,failedInds] = GetMidlinesByThinning(IM_proc_crop,...
    'fishPos',fp,'process','parallel','plotBool',1,'kerSize',9);

toc

if ~exist('nFramesInTrl')
    nFramesInTrl = 750; %(Default when collecting images at 500 fps)
elseif isempty(nFramesInTrl)
    nFramesInTrl = 750;
end

nTrls = size(IM_proc_crop,3)/nFramesInTrl;
trlStartFrames =  (0:nTrls-1)*750 + 1;
[tailCurv, tailCurv_uncorrected] = SmoothenMidlines(midlineInds,IM_proc_crop,3,'plotBool',0,...
    'pauseDur',0,'smoothFactor',8,'dsVecs',dsVecs,'trlStartFrames',trlStartFrames);

% orientation = GetFishOrientationFromMidlineInds(midlineInds,imgDims(1:2),'s');
% orientation_backup = orientation;

disp('Saving midline inds, and tailCurv...')
% procData.orientation = orientation;
procData.hOr_crop = hOr_crop;
procData.midlineInds = midlineInds;
procData.dsVecs = dsVecs;
procData.tailCurv = tailCurv;
procData.nFramesInTrl = nFramesInTrl;
toc

%% Plot trialized tail bends
TrializeTailBendsAndPlot

%% Motion Info
% motionThr = 1;
% % [motionFrames, swimStartFrames] = GetMotionFrames(fishPos,motionThr);
% motionInfo = GetMotionInfo(fishPos,orientation,imgDims(1),'motionThr',motionThr);

%% Saving processed images
saveOrNot = 'y';

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
    procData.IM_crop = IM_crop;
    toc
else
    disp('Data not saved!')
end
toc



end
