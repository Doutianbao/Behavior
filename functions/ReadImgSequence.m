function IM = ReadImgSequence(imgDir,varargin)
%ReadImgSequence Reads a sequences of images with a specified extension in
%   a specfied dir
% IM = ReadImgSequence(imgDir)
% IM = ReadImgSequence(imgDir,imgExt)
% IM = ReadImgSequence(imgDir,imgExt,imgNums)
% Inputs:
% imgDir - Directory where images are stored
% imgExt - Extension of the img, e.g., {'bmp'}['jpg']['tif'] ['png']
% imgNums - Vector specifying the subset of images to read from all the
%   images in the directory. For e.g., imgNuums = [1 3 5], will only cause
%   the images 1, 3, and 5 to be read from the image sequence
% imgInds - A vector of indices indicating which images in the entire
%   sequence to read


imgInds = [];
imgExt = 'bmp';
poolSize = 10;

if nargin == 2
    imgExt = varargin{1};
elseif nargin ==3
    imgExt = varargin{1};
    imgInds = varargin{2};
end
disp('Scanning all image files in the dir...')
tic
searchToken = ['*.' imgExt];
files = dir(fullfile(imgDir,searchToken));
fNames = {files.name};
if isempty(fNames)
    error('No files found in directory, please check path!')
end

if ~isempty(imgInds)
    fNames = fNames(imgInds);
end

imgInfo = imfinfo(fullfile(imgDir,fNames{1}));
imSize = [imgInfo.Height imgInfo.Width];
IM = zeros(imSize(1),imSize(2),length(fNames));
% IM = MappedTensor(imSize(1),imSize(2),length(fNames));
disp(['Reading all .' imgExt ' images from dir...'])
if matlabpool('size')==0
    matlabpool(poolSize)
end
imgNums = 1:length(fNames);
parfor jj = imgNums    
    img = imread(fullfile(imgDir,fNames{jj}));
    if length(size(img))==3
        if jj ==1
            disp('Great! Images are rgb and are being converted to gray...')
        end
        img = rgb2Gray(img);
    end
    IM(:,:,jj) = img;
end 
toc

end

function img_gray = rgb2Gray(img_rgb)
img_rgb(:,:,1) = img_rgb(:,:,1)*0.299;
img_rgb(:,:,2) = img_rgb(:,:,2)*0.587;
img_rgb(:,:,3) = img_rgb(:,:,3)*0.114;
img_gray = sum(img_rgb,3);
end
