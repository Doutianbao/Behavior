function PlayFishTracking(IM,varargin)
% PlayFishOrientation - Plays the video of the fish along with orientation
%   line
% PlayFishOrientation(IM, fishPos,orientation, 'startFrame',startFrame,...
%           'endFrame',endFrame,'pauseDur',pauseDur,'curv',curv, 'periPxls',periPxls)
%
% Inputs:
% IM - Image stack of size M x N x T, where M = image height, N = image
%   width, T = # of time points
% 'fishPos' - T x 2 vec where 1st & 2nd cols contain row & col positions
%   of the fish head centroid in the images
% 'midlineInds' - A cell array of length T, where each array contains indices
%   of the fish midline 
% 'frameInds' - Indices of frames to display
% 'pauseDur' - Duration of pause between subsequent frames
% 'plotCurv' - 0 or 1; 1 results in plotting of curvatures in a subplot
%   below the fish images
% saveDir - Directory to save images to

plotCurv = 0;
frameInds  = [];
saveDir = [];
pauseDur = 0;
tailAngles = [];

for jj = 1:numel(varargin)
    if isstr(varargin{jj})
        switch lower(varargin{jj})
            case 'frameinds'
                frameInds = varargin{jj+1};
            case 'fishpos'
                fishPos = varargin{jj+1};
                if numel(fishPos)==2
                    fishPos = fishPos(:)';
                    fishPos = repmat(fishPos,size(IM,3),1);
                end
            case 'pausedur'
                pauseDur = varargin{jj+1};
            case 'tailcurv'
                tailCurv = varargin{jj+1};
            case 'savedir'
                saveDir = varargin{jj+1};
            case 'plotcurv'
                plotCurv = varargin{jj+1};
            case 'tailangles'
                tailAngles = varargin{jj+1};
        end
    end
end

if isempty(frameInds)
    frameInds = 1:size(IM,3);
end
if isempty(tailAngles)
    tailAngles = GetTailTangents(tailCurv);
end

IM = IM(:,:,frameInds);
fishPos = fishPos(frameInds,:);
tailCurv = tailCurv(:,:,frameInds);
tailAngles = tailAngles(:,frameInds);
imgDims =[size(IM,1), size(IM,2)];
cLim = [min(IM(:)), max(IM(:))*0.9];

figure('Name', 'Fish Tracking')
count = 0;
blah = zeros(size(tailAngles,1),100);
for imgNum = frameInds(:)'
    count = count + 1;
   x = tailCurv(:,1,count);
   y = tailCurv(:,2,count);  
    if plotCurv
        ax1 = subplot(2,1,1);
        cla
        imagesc(IM(:,:,count),'parent',ax1); axis image, axis on, colormap(ax1,jet),axis off
        hold on
        set(gca,'clim',cLim)
        plot(x,y,'color','r','linewidth',2),drawnow
        title(['Frame: ' num2str(imgNum) ', Angle: ' num2str(round(tailAngles(end,count)*10)/10) '^o '])  
        if imgNum >= 100            
           blah = tailAngles(:,count-99:count);           
        else            
            blah = zeros(size(tailAngles,1),100);
            blah(:,end-count+1:end) = tailAngles(:,1:count);          
        end    
        ax2 = subplot(2,1,2);
        cla
        imagesc(blah,'parent',ax2); colormap(ax2,jet),set(gca,'clim',[-150 150]);
        box off, axis off
        pause(pauseDur)
    else
        cla
        imagesc(IM(:,:,imgNum)),axis image, axis on, colormap(gray)
        hold on
        set(gca,'clim',cLim)
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
        title(['Frame: ' num2str(imgNum) ', Angle: ' num2str(round(orientation(imgNum))) '^o ' ...
            ', Frame rate = ' num2str(round(1/eTime)) ' fps'])
        shg
        pause(pauseDur)
    end
   
end

end


