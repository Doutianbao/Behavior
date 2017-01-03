function varargout = ProcessFishImages(varargin)
% ProcessFishImages(pathList);
% procData = ProcessFishImages(pathList);
% procData = ProcessFishImages(pathList,'readMode',readMode,'fps',fps,'imgExt',imgExt,'nFramesInTrl',nFramesInTrl,...
%   'nHeadPxls',nHeadPxls,'lineLen',lineLen,'spatialFilt',spatialFilt,'imgInds',imgInds,'blockSize',blockSize,'cropWid',cropWid);
% Inputs:
% If no inputs are specified, allows for interactive selection of image
%   directory and use default parameters for processing.
% pathList - Cell array of paths, where each path points to the location of
%   a set of images such as all vibration or dark flash trials from an
%   experiment.
% readMode - 'fromImages' or 'fromMishVid'. The former results in reading
%   of an image sequence, whereas the latter reads .mishVid created by
%   Misha Ahrens.
% imgExt - Image extension; reads only images in specified folder with this
%   extension ['bmp']
% fps - Frames per second [500].
% nFramesInTrl = Number of frames in a single trial [750];
% spatialFilt - Spatial filter used to smooth images before finding fish's
%   head centroid
% imgInds - Vector of image indices to read. If empty, reads all images.
%   By default reads all images.
% nHeadPxls - Number of head pixels; This is the number of pixels to use
%   for finding head centroid of fish [25].
% blockSize - Loading all images can take up too much memory, so this
%   allows processing by splitting total number of images in this many
%   blocks and processing one block at a time. If empty, or by default
%   blockSize = 4.
% cropWid - Crop width. After finding fish, crops image by this width
%   around the fish to save storage space.
% lineLen - Length of line in pixels to use for estimating head orienation
%   [15].
% Outputs:
% procData - Mat file saving all relevant info after processing.
% 
% Avinash Pujala, Koyama lab/HHMI, 2016

readMode =  'fromImages';
imgExt = 'bmp';
imgInds = [];
fps = 500;
nFramesInTrl = 750;
spatialFilt = 30;
nHeadPxls = 25;
lineLen = 15;
blockSize = 4;
cropWid = 90; %( For imgDims ~ [900,900])

for jj  = 1:nargin
    if ischar(varargin{jj})        
        switch lower(varargin{jj})
            case 'readmode'
                readMode = varargin{jj+1};
            case 'imgext'
                imgExt = varargin{jj+1};
            case 'imginds'
                imgInds = varargin{jj+1};
            case 'fps'
                fps = val;
            case 'nframesintrl'
                nFramesInTral = varargin{jj+1};
            case 'spatialfilt'
                spatialFilt = varargin{jj+1};
            case 'nheadpxls'
                nHeadPxls = varargin{jj+1};
            case 'linelen'
                lineLen = varargin{jj+1};
            case 'blocksize'
                blockSize = varargin{jj+1};
            case 'cropwid'
                cropWid = varargin{jj+1};
        end
    end
end

pathList = varargin{1};
if ~iscell(pathList)
    pathList = {pathList}; % Encapsulate single path string in a cell.
    
end
procData = cell(length(pathList),1);

tic
for pp = 1:length(pathList)
    [currPath,currFile] = fileparts(pathList{pp});
   if ~strcmpi(currFile,'proc')
      currPath = fullfile(currPath,currFile);
   end
    fprintf('Processing ... \n')
    disp(currPath)
    procData{pp} = ProcessFishImages_subset(currPath,'readMode',readMode,'imgInds',imgInds,'fps',fps,...
        'imgExt',imgExt,'nFramesInTrl',nFramesInTrl,'nHeadPxls',nHeadPxls,'lineLen',lineLen,...
        'spatialFilt',spatialFilt,'blockSize',blockSize,'cropWid',cropWid);
    toc
end
toc

varargout{1}= procData;

end


function varargout = ProcessFishImages_subset(varargin)
% ProcessFishImages();
% procData = ProcessFishImages();
% procData = ProcessFishImages(imgDir);
% procData = ProcessFishImages(imgDir,'readMode',readMode,'fps',fps,'imgExt',imgExt,'nFramesInTrl',nFramesInTrl,...
%   'nHeadPxls',nHeadPxls,'lineLen',lineLen,'spatialFilt',spatialFilt,'imgInds',imgInds,'blockSize',blockSize,'cropWid',cropWid);
% Inputs:
% If no inputs are specified, allows for interactive selection of image
%   directory and use default parameters for processing.
% imgDir - Path to image directory. If empty, allows for interactive
%   selection of path by selecting any file in the directory.
% readMode - 'fromImages' or 'fromMishVid'. The former results in reading
%   of an image sequence, whereas the latter reads .mishVid created by
%   Misha Ahrens.
% imgExt - Image extension; reads only images in specified folder with this
%   extension ['bmp']
% fps - Frames per second [500].
% nFramesInTrl = Number of frames in a single trial [750];
% spatialFilt - Spatial filter used to smooth images before finding fish's
%   head centroid
% imgInds - Vector of image indices to read. If empty, reads all images.
%   By default reads all images
% nHeadPxls - Number of head pixels; This is the number of pixels to use
%   for finding head centroid of fish [25]
% blockSize - Loading all images can take up too much memory, so this
%   allows processing by splitting total number of images in this many
%   blocks and processing one block at a time. If empty, or by default
%   blockSize = 4.
% cropWid - Crop width. After finding fish, crops image by this width after
%   aroung the fish to save space
% lineLen - Length of line in pixels to use for estimating head orienation
%   [15].
% Outputs:
% procData - Mat file saving all relevant info after processing
% 
% Avinash Pujala, Koyama lab/HHMI, 2016

imgDir = [];
readMode =  'fromImages';
imgExt = 'bmp';
imgInds = [];
fps = 500;
nFramesInTrl = 750;
spatialFilt = 30;
nHeadPxls = 25;
lineLen = 15;
blockSize = 4;
cropWid = 90; %( For imgDims ~ [900,900])
poolSize  = 10;

for jj  = 2:nargin
    if ischar(varargin{jj})
        val = varargin{jj+1};
        switch lower(varargin{jj})
            case 'readmode'
                readMode = val;
            case 'imgext'
                imgExt = val;
            case 'imginds'
                imgInds = val;
            case 'fps'
                fps = val;
            case 'nframesintrl'
                nFramesInTral = val;
            case 'spatialfilt'
                spatialFilt = val;
            case 'nheadpxls'
                nHeadPxls = val;
            case 'linelen'
                lineLen = val;
            case 'blocksize'
                blockSize = val;
            case 'cropwid'
                cropWid = val;
        end
    end
end

imgExt(strfind(imgExt,'.'))=[];
if nargin ==0
    [~, imgDir] = uigetfile(['*.' imgExt]);
else
    imgDir = varargin{1};
    if isempty(imgDir)
        [~,imgDir] = uigetfile(['*.' imgExt]);
    end
end

switch readMode
    case 'fromMishVid'       
        [IM, outDir] = ReadMishVid();
        imgInds = 1:size(IM,3);
    case 'fromImages'           
        IM = ReadImgSequence(imgDir,imgExt,imgInds);
end

outDir = fullfile(imgDir,'proc');
if ~exist(outDir,'dir')
    mkdir(outDir)
end


%% Processing images block-by-block
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
    if ~isempty(spatialFilt)
        [fp,hOr_temp] = GetFishPos(im_proc, nHeadPxls,'filter',spatialFilt,'process','parallel',...
            'lineLen',lineLen);
    else
        [fp,hOr_temp] = GetFishPos(im_proc, nHeadPxls,'process','parallel','lineLen',lineLen);
    end
    fishPos(imgFrames,:) = fp;
    hOr(imgFrames) = hOr_temp;
    disp('Cropping images...')
    IM_proc_crop(:,:,imgFrames) = CropImgsAroundPxl(im_proc,fp,cropWid);
end
ref = mean(ref,3);
clear im_proc fp hOr_temp
toc

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
    [~,hOr_crop] =  GetFishPos(IM_proc_crop, nHeadPxls,'filter',spatialFilt,'process','parallel','lineLen',lineLen);
    toc
end

toc
disp('Saving fish position, ref image, and cropped images....')
if isempty(imgInds)
    imgInds = 1:size(IM,3);
end
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

nTrls = size(IM_proc_crop,3)/nFramesInTrl;
trlStartFrames =  (0:nTrls-1)*750 + 1;
[tailCurv, tailCurv_uncorrected] = SmoothenMidlines(midlineInds,IM_proc_crop,3,'plotBool',0,...
    'pauseDur',0,'smoothFactor',8,'dsVecs',dsVecs,'trlStartFrames',trlStartFrames);

disp('Saving midline inds, and tailCurv...')
procData.hOr_crop = hOr_crop;
procData.midlineInds = midlineInds;
procData.dsVecs = dsVecs;
procData.tailCurv = tailCurv;
procData.nFramesInTrl = nFramesInTrl;
toc

%% Plot trialized tail bends
close all
TrializeTailBendsAndPlot
disp('Saving trialized tail curvature figures figures...')
suffix = ['_Trialized tail curvatures_' ts];
h = get(0,'children');
for fig = 1:length(h)
    prefix = ['Fig_' sprintf('%0.2d', fig)];
    saveas(h(fig), fullfile(outDir,[prefix, suffix]))
end


%% Saving processed images
saveOrNot = 'y';
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

varargout{1} = procData;

end
