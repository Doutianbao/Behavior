function IM = ReadImgSequence(imgDir,varargin)
%ReadImgSequence Reads a sequences of images with a specified extension in
%   a specfied dir
% IM = ReadImgSequence(imgDir,imgExt)
% Inputs:
% imgDir - Directory where images are stored
% imgExt - Extension of the img, e.g., {'jpg'} ['tif'] ['png']

imgExt = 'jpg';
if nargin == 2
    imgExt = varargin{1}
end

files = dir(imgDir);
fNames = {files.name};
if isempty(fNames)
    error('No files found in directory, please check path!')
end
remInds = [];
for jj = 1:length(fNames)
    if isempty(findstr(lower(fNames{jj}),lower(imgExt)))
        remInds = [remInds;jj];
    end
end
fNames(remInds)=[];

imgInfo = imfinfo(fullfile(imgDir,fNames{1}));
imSize = [imgInfo.Height imgInfo.Width];
IM = zeros(imSize(1),imSize(2),length(fNames));
disp(['Reading all .' imgExt ' images from dir...'])
tic
for jj = 1 :length(fNames)
    IM(:,:,jj) = imread(fullfile(imgDir,fNames{jj}));
    if mod(jj,50)==0
        disp(['Img # ' num2str(jj)])
    end
end 
toc

end

