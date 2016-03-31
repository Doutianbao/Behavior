function IM_proc = ProcessImages_parallel(varargin)
%ProcessImages_parallel Smooth and background subtract images
% IM_proc = ProcessImages(IM)
% IM_proc = ProcessImages(IM,refFrames); 
% Inputs:
% IM - Image stack where 3rd time is time
% refFrames - A vector specifying the frames to average for background
%   subtraction. This is useful for considering only frames in which fish
%   moves

poolSize = 10;
IM = varargin{1};
disp('Computing smoothed mean frame...')
tic
imKer = ones(5)/5^2;
refFrames = 1:size(IM,3);
if (nargin == 2)
    refFrames = varargin{2};
end

im = conv2(squeeze(mean(IM(:,:,refFrames),3)),imKer,'same');
toc


disp('Processing images...')
tic
IM_proc = zeros(size(IM));
imgFrames= 1:size(IM,3);
poolOpened = 0;
if matlabpool('size')==0
 matlabpool(poolSize)
 poolOpened = 1;
end
dispChunk = round(numel(imgFrames)/30);
parfor jj=imgFrames
    IM_proc(:,:,jj) = conv2(-squeeze(IM(:,:,jj))+im,ones(5)/25,'same');
   if mod(jj, dispChunk)==0
       disp(['Img# ' num2str(jj)])
   end
end
if poolOpened
matlabpool close
end
% im = repmat(mean(IM_proc,3),[1 1 size(IM_proc,3)]);
% IM_proc = IM_proc-im;
toc
end
