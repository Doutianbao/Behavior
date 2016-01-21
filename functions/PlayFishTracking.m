function PlayFishTracking(IM,fishPosVec,varargin)
% PlayFishTracking - Plays the video of the fish being tracked
% PlayFishTracking(IM, fishPos,startFrame,endFrame,pauseDur)
% Inputs:
% IM - Image stack of size M x N x T, where M = image height, N = image
%   width, T = # of time points
% fishPosVec - T x 2 vec where 1st & 2nd cols contain row & col positions
%   of the fish in the images
% startFrame - Frame from which to display video
% endFrame - Last Frame of video display
% pauseDur - Duration of pause between subsequent frames

if nargin == 2
    startFrame = 1;
    endFrame = 1000;
    pauseDur = 0.1;
elseif nargin == 3;
    startFrame = varargin{1};
    endFrame = startFrame + 999;
    pauseDur = 0.1;
elseif nargin == 4;
    startFrame = varargin{1};
    endFrame = varargin{2};
    pauseDur = 0.1;
elseif nargin == 5
    startFrame = varargin{1};
    endFrame = varargin{2};
    pauseDur = varargin{3};
else
    error('Too many inputs!')
end


figure('Name', 'Fish Tracking')
tic
for imgNum = startFrame:endFrame
    imagesc(IM(:,:,imgNum)),axis image, axis off, colormap(gray), drawnow
    hold on
    plot(fishPosVec(1,imgNum),fishPosVec(2,imgNum),'ro'), drawnow
    title(['Frame: ' num2str(imgNum) ', Frame Rate: ' num2str((round(10*((imgNum+1-startFrame)/toc)))/10) ' fps'])
    shg
    pause(pauseDur)
    cla
end
