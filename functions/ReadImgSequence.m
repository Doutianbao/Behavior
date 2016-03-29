function IM = ReadImgSequence(imgDir,varargin)
%ReadImgSequence Reads a sequences of images with a specified extension in
%   a specfied dir
% IM = ReadImgSequence(imgDir,imgExt)
% Inputs:
% imgDir - Directory where images are stored
% imgExt - Extension of the img, e.g., {'jpg'} ['tif'] ['png']
% imgInds - A vector of indices indicating which images in the entire
%   sequence to read

imgInds = [];
imgExt = 'jpg';
if nargin == 2
    imgExt = varargin{1}
elseif nargin == 3
    imgExt = varargin{1};
    imgInds = varargin{2};
end
disp('Scanning all image files in the dir...')
tic
files = dir(imgDir);
fNames = {files.name};
if isempty(fNames)
    error('No files found in directory, please check path!')
end
remInds = [];
dispChunk = round(length(fNames)/30);
for jj = 1:length(fNames)
    if isempty(findstr(lower(fNames{jj}),lower(imgExt)))
        remInds = [remInds;jj];
    end
    if mod(jj,dispChunk)==0
        disp([ num2str(jj) '  images...'])
    end
end
disp([num2str(length(fNames)) ' images found'])
toc
disp('Subsampling image sequence based on input indices...')
fNames(remInds)=[];
if ~isempty(imgInds)
    fNames = fNames(imgInds);
end

imgInfo = imfinfo(fullfile(imgDir,fNames{1}));
imSize = [imgInfo.Height imgInfo.Width];
IM = zeros(imSize(1),imSize(2),length(fNames));
disp(['Reading all .' imgExt ' images from dir...'])
tic
for jj = 1 :length(fNames)
    img = imread(fullfile(imgDir,fNames{jj}));
    if length(size(img)==3)
        if jj ==1
            disp('Great! Images are rgb and are being converted to gray...')
        end
        img = rgb2Gray(img);
    end
    IM(:,:,jj) = img;
    if mod(jj,50)==0
        disp(['Img # ' num2str(jj)])
    end
end 
toc

end

function img_gray = rgb2Gray(img_rgb)
img_rgb(:,:,1) = img_rgb(:,:,1)*0.299;
img_rgb(:,:,2) = img_rgb(:,:,2)*0.587;
img_rgb(:,:,3) = img_rgb(:,:,3)*0.114;
img_gray = sum(img_rgb,3);
end
