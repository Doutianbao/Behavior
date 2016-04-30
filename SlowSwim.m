

% switch  'LoadingNewFilm'  %'LoadingNewFilm' 'RerunAnalysis' 'LoadingCoordinates'
%     case 'LoadingNewFilm'
%

clear all
close all

cd 'Z:\Avinash\Ablations & Behavior'

readMode =  'fromImages'; 
% readMode = 'fromMishVid';

poolSize  = 12;
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
%         IM = ReadImgSequence_beta(imgDir,imgExt,imgInds);
         IM = ReadImgSequence(imgDir,imgExt,imgInds);
        outDir = fullfile(imgDir,'spont');
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
IM_proc = ProcessImages(IM);

%% Tracking the fish
disp('Tracking fish...')
fishPos = GetFishPos(IM_proc, 60);
toc
disp('Creating mean reference frame...')
ref = mean(IM,3);
toc

%% Fish Orientation
disp('Getting fish orientation...')
tic

% midlineInds = GetMidline_template_parallel(IM_orient,fishPos,[30]);
midlineInds = GetMidline_beta(IM_proc,fishPos,[32 20]);
% midlineInds = GetMidline_beta(IM_proc,fishPos,[32 20]);

%    orientation_corr = CorrectOrientation(orientation, 90);
imgDims = size(IM_proc);
orientation = GetFishOrientationFromMidlineInds(midlineInds,imgDims(1:2),'s');
orientation = orientation';
orientation_backup = orientation;
toc

%% Motion Info
% motionThr = 5;
% [motionFrames, swimStartFrames] = GetMotionFrames(fishPos,motionThr);
% motionInfo = GetMotionInfo(fishPos,orientation,imgDims(1));

%% Save timeseries
% fName = input('Enter fish name (e.g. Fish7): ','s');
fName_suffix = [num2str(round(imgInds(1)/60/30)) '-' num2str(round(imgInds(end)/60/30)) 'mins'];
fName = strcat(fName_prefix,'_',fName_suffix);
ts = datestr(now,30);
save(fullfile(outDir,[fName, '_orientation_' ts '.mat']),'orientation');
save(fullfile(outDir,[fName, '_imgDims_'  ts '.mat']),'imgDims')
save(fullfile(outDir,[fName, '_midlineInds_' ts '.mat']),'midlineInds');
save(fullfile(outDir,[fName, '_ref_' ts '.mat']),'ref');
save(fullfile(outDir,[fName, '_tracexy_' ts '.mat']),'fishPos');
% save(fullfile(outDir,[fName '_motionInfo_' ts '.mat']),'motionInfo');
disp(['Saved orientation, imgDims ,midlineInds, ref, tracexy at ' outDir])

break;

%% Saving processed images
saveOrNot = 'y';
%         saveOrNot = input('Save the variables (y/n)?  ','s');
tic
if strcmpi('y',saveOrNot)
    disp('Saving relevant variables...') 
    savefast(fullfile(outDir,[fName, '_IM_proc.mat']),'IM_proc');    
else
    disp('Data not saved!')
end
toc
break;
if matlabpool('size')>0
    matlabpool close
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

