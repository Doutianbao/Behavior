function PlayFixedFish(IM,fishPos,orientation,varargin)
% PlayFixedFish - Plays the video of the fish after adjustment for position
%   and orientation
% PlayFixedFish(IM,fishPos,orientation)
% PlayFixedFish(...,frameInds,pauseDur)
% Inputs:
% IM - Image stack of size M x N x T, where M = image height, N = image
%   width, T = # of time points
% fishPos - T x 2 vec where 1st & 2nd cols contain row & col positions.
%   fishPos = [], results in playing the frames without position tracking
%   of the fish in the images
% frameInds - Indices of frames to play; frameInds  = [] results in playing all frames
% pauseDur - Duration of pause between subsequent frames

frameInds = 1:size(IM,3);
pauseDur = 0.1;
if nargin < 3
    error('Minimum 3 inputs not entered!')    
elseif nargin == 4;
    frameInds = varargin{1};   
elseif nargin == 5;
    frameInds = varargin{1};
    pauseDur = varargin{2};
elseif nargin > 5
    error('Too many inputs!')
end
if isempty(frameInds);
    frameInds = 1:size(IM,3);
end

rMid = floor(size(IM,1)/2);
cMid = floor(size(IM,2)/2);
% I_adj = zeros(size(IM,1),size(IM,2),numel(frameInds));
figure('Name','Fixed Fish Movie')
for imgNum = frameInds(:)'
    cla
    im = IM(:,:,imgNum);
%     im = GetFishPxls(im,fishPos(imgNum,:));
%     [~,inds] = GetFishPxls(im,fishPos(imgNum,:));
%     blah = zeros(size(im));
%     blah(inds) = im(inds);
%     im = blah;
    rShift = rMid-fishPos(imgNum,2);
    cShift = cMid-fishPos(imgNum,1);
    im_rot = imrotate(im,orientation(imgNum),'nearest','crop');
    im_trans = circshift(im,[rShift,cShift]);
    im_trans_rot = imrotate(im_trans,mod(orientation(imgNum)-90,360),'nearest','crop');
    im_trans_rot(im_trans_rot<-0.5) = nan;    
%     I_adj(:,:,count) = im_trans_rot;    
    imagesc(im_trans_rot), axis image, drawnow, colormap(gray)
    shg
    title(num2str(imgNum))
    pause(pauseDur)    
end






