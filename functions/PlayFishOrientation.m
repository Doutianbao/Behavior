function PlayFishOriention(IM,fishPos,orientation,varargin)
% PlayFishOrientation - Plays the video of the fish along with orientation
%   line
% PlayFishOrientation(IM, fishPos,orientation, startFrame,endFrame,pauseDur)
% PlayFishOrientation(IM,fishPos,midlineInds,startFrame,endFrame,pauseDur)
% Inputs:
% IM - Image stack of size M x N x T, where M = image height, N = image
%   width, T = # of time points
% fishPos - T x 2 vec where 1st & 2nd cols contain row & col positions
%   of the fish in the images
% orientation - Orientation of the fish in degrees
% midlineInds - Indices of the fish midline in the image
% startFrame - Frame from which to display video
% endFrame - Last Frame of video display
% pauseDur - Duration of pause between subsequent frames

lineLength = 25;
skipFrames = 1;

if nargin == 3
    startFrame = 1;
    endFrame = size(IM,3);
    pauseDur = 0.1;
    
elseif nargin == 4;
    startFrame = varargin{1};
    endFrame = startFrame + size(IM,3);
    pauseDur = 0.1;
    
elseif nargin == 5;
    startFrame = varargin{1};
    endFrame = varargin{2};
    pauseDur = 0.1;
    
elseif nargin == 6
    startFrame = varargin{1};
    endFrame = varargin{2};
    pauseDur = varargin{3};
    
elseif nargin >6
    error('Too many inputs!')
end

if iscell(orientation) % Implies this input is midlineInds, rather than orientation
    [orientation, origins] = GetFishOrientationFromMidlineInds(orientation);
end

figure('Name', 'Fish Tracking')
tic
rho = 0:lineLength-1;
for imgNum = startFrame:skipFrames:endFrame
    cla
    imagesc(IM(:,:,imgNum)),axis image, axis on, colormap(gray), drawnow
    hold on
    %     plot(fishPos(imgNum,1),fishPos(imgNum,2),'ro'), drawnow
    theta = orientation(imgNum);
    theta= mod(theta+180,360)*pi/180;
    
    
    thetas = repmat(theta,1,lineLength);
    [blahX,blahY] = pol2cart(thetas,rho);
    img = IM(:,:,imgNum);
    try
    blahX = blahX + fishPos(imgNum,1);
    blahY = blahY + fishPos(imgNum,2);
    catch
        a = 1;
    end
    
    
    plot(blahX,blahY,'color','r','linewidth',2),drawnow
    title(['Frame: ' num2str(imgNum) ', Angle: ' num2str(round(theta*180/pi)) '^o ' ...
        ', Frame Rate: ' num2str((round(10*((imgNum+1-startFrame)/toc)))/10) ' fps'])
    shg
    pause(pauseDur)
    % pause()
    
end
