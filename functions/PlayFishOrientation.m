function PlayFishOriention(IM,fishPos,orientation,varargin)
% PlayFishOrientation - Plays the video of the fish along with orientation
%   line
% PlayFishOrientation(IM, fishPos,orientation, 'startFrame',startFrame,...
%           'endFrame',endFrame,'pauseDur',pauseDur,'curv',curv, 'periPxls',periPxls)
%
% Inputs:
% IM - Image stack of size M x N x T, where M = image height, N = image
%   width, T = # of time points
% fishPos - T x 2 vec where 1st & 2nd cols contain row & col positions
%   of the fish in the images
% orientation - Orientation of the fish in degrees
% midlineInds - Indices of the fish midline in the image
% frameInds - Indices of frames to display
% pauseDur - Duration of pause between subsequent frames
% plotCurv - 0 or 1, 0 results in plotting of curvatures
% saveDir - Directory to save images to

lineLength = 25;
skipFrames = 1;
plotCurv = 0;
frameInds  = [];
midlineInds = [];
saveDir = [];
pauseDur =0;
if nargin < 3
    error('3 inputs required!')
end
for jj = 1:numel(varargin)
    if strcmpi(varargin{jj},'frameInds')
        frameInds = varargin{jj+1};
    end
    if strcmpi(varargin{jj},'pauseDur')
        pauseDur = varargin{jj+1};
    end
    if strcmpi(varargin{jj},'midlineInds')
        midlineInds = varargin{jj+1};
    end
    if strcmpi(varargin{jj},'saveDir')
        saveDir = varargin{jj+1};
    end
    if strcmpi(varargin{jj},'plotCurv')
        plotCurv = varargin{jj+1};
        curv = GetCurvInfo(orientation);
        curv = fix(curv/8)*8;
        ker = gausswin(6); ker = ker/sum(ker);
        for kk = 1:size(curv,2)
            curv(:,kk) = conv2(curv(:,kk),ker(:),'same');
        end
   
    end
end


if iscell(orientation) % Implies this input is midlineInds, rather than orientation
    imgDims = size(IM);
    [orientation, origins] = GetFishOrientationFromMidlineInds(orientation,imgDims(1:2),'s');
end
orientation = orientation';

figure('Name', 'Fish Tracking')
tic
if isempty(frameInds)
    frameInds = 1:size(IM,3);
end
imgDims =[size(IM,1), size(IM,2)];
if ~isempty(midlineInds)
    clrs = jet(length(midlineInds{1}));
end
for imgNum = frameInds(:)'
    [x,y] = Or2LnInds(orientation(imgNum,1),lineLength);
    x = x + fishPos(imgNum,1);
    y = y + fishPos(imgNum,2);
    if plotCurv
        subplot(2,1,1)
        cla
        imagesc(IM(:,:,imgNum)),axis image, axis on, colormap(gray)
        hold on
        plot(x,y,'color','r','linewidth',2),drawnow
        eTime = toc;
        title(['Frame: ' num2str(imgNum) ', Angle: ' num2str(round(orientation(imgNum))) '^o ' ...
            ', Frame rate = ' num2str(round(1/eTime)) ' fps'])
        subplot(2,1,2)
        if imgNum >= 100
            cla
            set(gca,'color','k'), hold on
            plot(imgNum-99:imgNum,curv(imgNum-99:imgNum,1),'r')
            plot(imgNum-99:imgNum,curv(imgNum-99:imgNum,2),'g')
            plot(imgNum-99:imgNum,ones(1,100)*10,'w:')
            plot(imgNum-99:imgNum,ones(1,100)*-10,'w:')
            xlim([imgNum-99 imgNum])          
            %             plot([imgNum, imgNum ],[-20,20],'c--')
        else
            cla
            set(gca,'color','k'), hold on
            plot(curv(1:imgNum,1),'r')
            plot(curv(1:imgNum,2),'g')
            plot(ones(imgNum,1)*10,'w:')
            plot(ones(imgNum,1)*-10,'w:')
            xlim([0 imgNum])            
        end
        ylim([-120 120])
        lh = legend('head','tail');
        set(lh,'color','w','location','NorthWest')
        shg
        pause(pauseDur)
    else
        cla
        imagesc(IM(:,:,imgNum)),axis image, axis on, colormap(gray)
        hold on
        if ~isempty(midlineInds)
            for line = 1:length(midlineInds{imgNum})
            [y,x] = ind2sub(imgDims,midlineInds{imgNum}{line});
            plot(x,y,'color',clrs(line,:),'linewidth',2)          
            end            
        else
            plot(x,y,'color','r','linewidth',2)
        end   
        axis off
        if ~isempty(saveDir)            
            suffix = sprintf('%0.5d',imgNum);
            fName = ['Img_' suffix];
            print(fullfile(saveDir,fName),'-dpng')
        end
        drawnow
        eTime = toc;
        title(['Frame: ' num2str(imgNum) ', Angle: ' num2str(round(orientation(imgNum))) '^o ' ...
            ', Frame rate = ' num2str(round(1/eTime)) ' fps'])
        shg
        pause(pauseDur)
    end
    tic
end

end

function [x,y] = Or2LnInds(orientation,lineLength)
rho = 0:lineLength-1;
theta = orientation;
theta= mod(theta+180,360)*pi/180;
thetas = repmat(theta,1,lineLength);
[x,y] = pol2cart(thetas,rho);
end

function curv = GetCurvInfo(orientation)
if (size(orientation,1)<size(orientation,2))
    orientation = orientation';
end
curv = nan(size(orientation,1),size(orientation,2)-1);
toVec = @(x)pol2cart(x*pi/180,ones(size(x)));
for seg = 1:size(curv,2)
    [x1,y1] = toVec(orientation(:,seg));
    [x2,y2] = toVec(orientation(:,seg+1));
    curv(:,seg)= angle((x1(:) + y1(:)*1i).* conj(x2(:) + y2(:)*1i));
end
curv = -curv*180/pi;
end