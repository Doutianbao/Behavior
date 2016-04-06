function varargout= PlayFixedFish(IM,fishPos,orientation,varargin)
% PlayFixedFish - Plays the video of the fish after adjustment for position
%   and orientation
% PlayFixedFish(IM,fishPos,orientation)
% PlayFixedFish(...,frameInds,pauseDur)
% IM_adj = PlayFixedFish(...)
% Inputs:
% IM - Image stack of size M x N x T, where M = image height, N = image
%   width, T = # of time points
% fishPos - T x 2 vec where 1st & 2nd cols contain row & col positions.
%   fishPos = [], results in playing the frames without position tracking
%   of the fish in the images
% frameInds - Indices of frames to play; frameInds  = [] results in playing all frames
% pauseDur - Duration of pause between subsequent frames
% Outputs:
% IM_adj = Adjusted IM, where each frame is transformed to account for
%   position and orientation of fish

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
if nargout > 0
    I_adj = zeros(size(IM,1),size(IM,2),numel(frameInds));
end
figure('Name','Fixed Fish Movie')
count = 0;
for imgNum = frameInds(:)'
    cla
    im = IM(:,:,imgNum);
    rShift = rMid-fishPos(imgNum,2);
    cShift = cMid-fishPos(imgNum,1);
    im_rot = imrotate(im,orientation(imgNum),'nearest','crop');
    im_trans = circshift(im,[round(rShift),round(cShift)]);
    im_trans_rot = imrotate(im_trans,mod(orientation(imgNum)-90,360),'nearest','crop');
    im_trans_rot(im_trans_rot<-0.5) = nan; 
    count = count + 1;
    if nargout > 0
    I_adj(:,:,count) = im_trans_rot;
    end
    imagesc(im_trans_rot), axis image, drawnow, colormap(gray)
    shg
    title(num2str(imgNum))
    pause(pauseDur)    
end

varargout{1} = I_adj;





