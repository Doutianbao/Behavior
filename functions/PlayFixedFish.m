function varargout= PlayFixedFish(IM,fishPos,varargin)
% PlayFixedFish - Plays the video of the fish after adjustment for position
%   and orientation
% PlayFixedFish(IM,fishPos)
% PlayFixedFish(IM,fishPos,orientation)
% PlayFixedFish(...,frameInds,pauseDur)
% PlayFixedFish(...,'frameInds',frameInds,'pauseDur',pauseDur,'zoom',zoom)
% IM_adj = PlayFixedFish(...)
% Inputs:
% IM - Image stack of size M x N x T, where M = image height, N = image
%   width, T = # of time points
% fishPos - T x 2 vec where 1st & 2nd cols contain row & col positions.
%   fishPos = [], does not center the fish in the center of the image
%   of the fish in the images
% orientation - Orientation of the fish's head; If orientation = [], does
%   not rotation correct img
% frameInds - Indices of frames to play; frameInds  = [] results in playing all frames
% pauseDur - Duration of pause between subsequent frames
% periPxls - # of pixels row and column pixels to keep on either side of
%   the fish when cropping image. periPxls = [], results in displaying of
%   full image
% Outputs:
% IM_adj = Adjusted IM, where each frame is transformed to account for
%   position and orientation of fish

frameInds = 1:size(IM,3);
pauseDur = 0.1;
orientation = [];
periPxls = [];
if nargin < 2
    error('Minimum 2 inputs not entered!')
elseif nargin > 2
    orientation = varargin{1};
end
if isempty(frameInds);
    frameInds = 1:size(IM,3);
end

for jj = 1:length(varargin)
    if strcmpi(varargin{jj},'frameInds')
        frameInds = varargin{jj+1};
    end
    if strcmpi(varargin{jj},'pauseDur')
        pauseDur = varargin{jj+1};
    end
    if strcmpi(varargin{jj},'periPxls')
        periPxls = round(varargin{jj+1});
    end
end

rMid = floor(size(IM,1)/2);
cMid = floor(size(IM,2)/2);
if ~isempty(periPxls)
    rInds = [max(rMid - periPxls,1), min(rMid + periPxls,size(IM,1))];
    cInds = [max(cMid - periPxls,1), cMid + min(periPxls,size(IM,3))];
else
    rInds = [1,size(IM,1)];
    cInds = [1,size(IM,2)];
end

if nargout > 0
    I_adj = zeros(size(IM,1),size(IM,2),numel(frameInds));
end
figure('Name','Fixed Fish Movie')
count = 0;
for imgNum = frameInds(:)'
    cla
    im = IM(:,:,imgNum);
    if isempty(fishPos)
        im_trans = im;
    else
        rShift = rMid-fishPos(imgNum,2);
        cShift = cMid-fishPos(imgNum,1);
        im_trans = circshift(im,[round(rShift),round(cShift)]);
    end
    if isempty(orientation)
        im_trans_rot = im_trans;
    else
        im_trans_rot = imrotate(im_trans,mod(orientation(imgNum)-90,360),'nearest','crop');
    end
    %     im_trans_rot(im_trans_rot<-0.5) = nan;
    count = count + 1;
    im_trans_rot = im_trans_rot(rInds(1):rInds(2),cInds(1):cInds(2));
    imagesc(im_trans_rot), axis image, drawnow, colormap(gray)
    shg
    title(num2str(imgNum))
    pause(pauseDur)
    if nargout > 0
        I_adj(:,:,count) = im_trans_rot;
    end
end

varargout{1} = I_adj;



