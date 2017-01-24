function PlayFishTrials(procData,trlList, varargin)
% PlayFishTrials - Plays the videos of the fish for specified trials with superimposed
%   midline and bend angle timeseries(all bend angles along the body or total bend angle)
% PlayFishOrientation(procData,trlList,'frameInds','frameInds','pauseDur',pauseDur,
%   'bendAngleDispMode',bendAngleDispMode,'writeVideo',writeVideo,'saveDir',saveDir);
% PlayFishTrials([],trlList,....);
%
% Inputs:
% procData - Processed data with all the relevant info created by
%   ProcessFishImages. If [], then selects for interactive selection of
%   procData;
% trlList - List of trials to play for a given fish. e.g., [1 3  10];
% frameInds - Indices of frames to play, which are a subset of the
%   specified trials. For instance, frameInds = 1:100, plays only the first
%   100 frames of each of the trials.
% pauseDur - The duration of pause between displaying successive frames.
% bendAngleDispMode - 'all','total','none';
%   'all' results in the displaying in a subplot under the images the of all the angles
%       along the length of the fish as colored matrix
%   'total' - displays the total bend angle from head to tail as a single
%       trace
%   'none' - does not display bend angle
% nDispPts - The maximum number of  points of the bend angle to display at
%   any given time; determines the scaling of the bend angle timeseries plotted below the
%   image frames. If thsi exceeds the total number of available frames,
%   defaults to the number of available frames.
% 'writeVideo' - 0 or 1; the latter causes the fish video to be saved in
%   the current directory, or another directory if specified by saveDir
% 'saveDir' - Directory in which to save written movies, if [], then writes
%   in current directory
%
% Avinash Pujala, Koyama lab/ JRC (HHMI), 2017

frameInds  = [];
pauseDur = 0;
bendAngleDispMode = 'all';
nDispPts = 150;
writeVideo = 0;
saveDir = [];
lineClr = 'r';



if nargin < 2
    error('Minimum 2 inputs required!')
end

if isempty(procData)
    procData = OpenMatFile();
end

for jj = 1:numel(varargin)
    if ischar(varargin{jj})
        switch lower(varargin{jj})
            case 'frameinds'
                frameInds = varargin{jj+1};
            case 'pausedur'
                pauseDur = varargin{jj+1};
            case lower('bendAngleDispMode')
                bendAngleDispMode = varargin{jj+1};
            case 'ndisppts'
                nDispPts = varargin{jj+1};
            case 'writevideo'
                writeVideo = varargin{jj+1};
            case 'savedir'
                saveDir = varargin{jj+1};
        end
    end
end

nFramesInTrl  = procData.nFramesInTrl;
if isempty(frameInds)
    frameInds = 1:nFramesInTrl;
end

nDispPts = min(max(1,nDispPts),numel(frameInds));

disp('Reading relevant info from procData...')
I = procData.IM_proc_crop;
tC = procData.tailCurv;
tA = GetTailTangents(tC);
curv = tA(end,:);
F= @(trlNum)(trlNum-1)*nFramesInTrl + 1: (trlNum-1)*nFramesInTrl + nFramesInTrl;
cLim = [0 0.9];
for trl = trlList(:)'
    disp(['Trl # ' num2str(trl)])
    if writeVideo
        if isempty(saveDir)
            saveDir = cd;
        end
        ts = datestr(now,30);
        fName = ['FishTrackingVideo_trl #' num2str(trl) '_' ts '.avi'];
        vidObj = VideoWriter(fullfile(saveDir,fName));
        open(vidObj);
        disp('Writing to video...')
    end
    inds_trl = F(trl);
    I_this = Standardize(I(:,:,inds_trl(frameInds)));
    tC_this = tC(:,:,inds_trl(frameInds));
    tA_this = tA(:,inds_trl(frameInds));
    curv_this = curv(inds_trl(frameInds));
    
    fh = figure('Name','Fish Tracking');
    %     set(fh,'units','normalized','pos',[0.2 0.2 0.5 0.5]);
    ax = cell(2,1);
    ax{1} = [0.6 0.6 0.2 0.35];
    ax{2}=[1 0.3 0 0];
    axH = CreateSubaxes(fh, ax{1}, ax{2});
    cMap = jet(64*4);
    
    count = 0;
    blah = zeros(size(tA,1),100);
    for imgNum = 1:length(frameInds)
        count = count + 1;
        x = tC_this(:,1,count);
        y = tC_this(:,2,count);
        if strcmpi(bendAngleDispMode,'all')
            figure(fh)
            axes(axH(1))
            imagesc(I_this(:,:,count),'parent',axH(1)),axis image, axis off
            colormap(jet)
            hold on
            set(gca,'clim',cLim)
            plot(x,y,'color',lineClr,'linewidth',2)
%             plot(x(1),y(1),'bo','markersize',10)
            %             plot(fishPos(count,1),fishPos(count,2),'g*','markersize',10)
            drawnow
            title(['Abs Frame: ' num2str(imgNum) ', Frame # ' num2str(count)...
                ', Angle: ' num2str(round(tC_this(end,count)*10)/10) '^o '])
            hold off
            if count >= nDispPts
                blah = tA_this(:,count-(nDispPts-1):count);
            else
                blah = zeros(size(tA_this,1),nDispPts);
                blah(:,end-count+1:end) = tA_this(:,1:count);
            end
            figure(fh)
            axes(axH(2))
            imagesc(1:size(blah,2),1:size(blah,1),blah,'parent',axH(2));
            colormap(axH(2),cMap)
            set(gca,'clim',[-150 150],'ytick',[5 size(blah,1)-5],'yticklabel',{'Head','Tail'},'xtick',[]);
            colorbar
            drawnow
            if isempty(pauseDur)
                pause()
            else
                pause(pauseDur)
            end
        elseif strcmpi(bendAngleDispMode,'total')
            figure(fh)
            axes(axH(1))
            imagesc(I_this(:,:,count),'parent',axH(1)),axis image, axis off
            colormap(jet)
            hold on
            set(gca,'clim',cLim)
            plot(x,y,'color',lineClr,'linewidth',2)
%             plot(x(1),y(1),'bo','markersize',10)
            %             plot(fishPos(count,1),fishPos(count,2),'g*','markersize',10)
            drawnow
            title(['Frame # ' num2str(count)])
            hold off
            
            axes(axH(2))
            box off
            inds = (count-nDispPts:count);
            inds(inds<=0)=[];
            plot(inds,curv_this(inds),'g')
            box off
            xlim([count-nDispPts count]);
            yl = [min(min(curv_this),-nDispPts) max(max(curv_this),nDispPts)];
            ylim(yl)
            set(gca,'color','k','tick','out','ytick',...
                [round((yl(1)*0.75)/10)*10 0 round((yl(2)*0.75)/10)*10],'xtick',[]);
            ylabel('Body bend angle (deg)')
            drawnow
            if isempty(pauseDur)
                pause()
            else
                pause(pauseDur)
            end
        else
            figure(fh)
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
            if isempty(pauseDur)
                pause()
            else
                pause(pauseDur)
            end
        end
        if writeVideo
            currFrame = getframe(fh);
            vidObj.writeVideo(currFrame);
        end
    end
    if writeVideo
        close(vidObj)
    end
end

end


