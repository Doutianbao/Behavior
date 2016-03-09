function IM_proc = ProcessImages(varargin)
%ProcessImages Smooth and background subtract images
% IM_proc = ProcessImages(IM)
% IM_proc = ProcessImages(IM,refFrames); 
% Inputs:
% IM - Image stack where 3rd time is time
% refFrames - A vector specifying the frames to average for background
%   subtraction. This is useful for considering only frames in which fish
%   moves

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

% nWorkers= 10;
disp('Processing images...')
tic
IM_proc = zeros(size(IM));
imgFrames= 1:size(IM,3);
% matlabpool(nWorkers)
for jj=imgFrames
    ii= conv2(-squeeze(IM(:,:,jj))+im,ones(5)/25,'same'); 
    IM_proc(:,:,jj) = ii;
  
    if mod(jj,500)==0
        disp(['Img # ' num2str(jj)])
    end
end
% matlabpool close
% im = repmat(mean(IM_proc,3),[1 1 size(IM_proc,3)]);
% IM_proc = IM_proc-im;
toc
end
