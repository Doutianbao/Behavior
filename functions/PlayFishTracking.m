function PlayFishTracking(IM,fishPos,varargin)
% PlayFishTracking - Plays the video of the fish being tracked
% PlayFishTracking(IM, fishPos,startFrame,endFrame,pauseDur)
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
if nargin == 3;
    frameInds = varargin{1};   
elseif nargin == 4;
    frameInds = varargin{1};
    pauseDur = varargin{2};
elseif nargin > 4
    error('Too many inputs!')
end
if frameInds == [];
    frameInds = 1:size(IM,3);
end
figure('Name', 'Fish Tracking')
tic
for imgNum = frameInds(:)'
    cla
    imagesc(IM(:,:,imgNum)),axis image, axis off, colormap(gray), drawnow
    hold on
    if ~isempty(fishPos)
    plot(fishPos(imgNum,1),fishPos(imgNum,2),'ro'), drawnow
    end
    title(['Frame: ' num2str(imgNum) ', Frame Rate: ' num2str((round(10*((imgNum+1-startFrame)/toc)))/10) ' fps'])
    shg
    pause(pauseDur)   
end
