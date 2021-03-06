%% Inputs for fast
fps = 500; % (30 for slow, 500 for fast)
nFramesInTrl = 750; %(1800 for slow, 750 for fast)
preStimPer = 0.1; % In seconds (1.5 for slow)
mult = 1; % Determines whether preStimPer gets added or subtracted while plotting rasters (-1 for slow)
nTrls = size(fishPos,1)/nFramesInTrl;
tapTrls = 1:2:nTrls;
flashTrls = 2:2:nTrls;
yOff = [220 200];
seg = 3; % [1 - head, 2 = tail, 3 - combined]
motionThr = 4;
peakDetThr = 0.2;

%## Sort all frame indices based on whether they are tap or dark trials
tapInds = [];
for trl = tapTrls
    tapInds = [tapInds, (trl-1)*nFramesInTrl + 1 : trl*nFramesInTrl];
end
darkInds = setdiff(1:length(time),tapInds);

% %% Inputs for slow
% fps = 30; % (30 for slow, 500 for fast)
% nFramesInTrl = 1800; %(1800/3600 for slow, 750 for fast)
% preStimPer = 0.6; % In seconds (1.5 for slow)
% mult = -1; % Determines whether preStimPer gets added or subtracted while plotting rasters (-1 for slow)
% nTrls = size(fishPos,1)/nFramesInTrl;
% tapTrls = 1:1:nTrls;
% flashTrls = 1:nTrls;
% yOff = [220 200];
% seg = 3; % [1 - head, 2 = tail, 3 - combined]
% motionThr = 4;
% peakDetThr = 0.1;

%% Getting motion info
imgDims = [size(ref), size(fishPos,1)];
motionInfo = GetMotionInfo(fishPos,orientation,imgDims(1));
time = (0:size(fishPos,1)-1)/fps;
dS_all = motionInfo.dS_all;
trlFrames = (0:nTrls-1)*nFramesInTrl + 1;
dS_all(trlFrames)=0; % Zeroing artifactual distances resulting from trial transitions

%% Correct fish pos, midlineInds, orientation
dS_thr = 50;
S = @(x)(sum(diff(x,[],1).^2,2)).^0.5;
badInds = find(dS_all>dS_thr);
badInds(badInds==1)=[];
badInds(badInds >= (length(dS_all)-1))=[];
iter = 0;
while numel(badInds)>0
    iter = iter + 1;
    disp(['Iter # ' num2str(iter)])
    if iter > 10
        break
    end
    fishPos(badInds,:) = 0.5*(fishPos(badInds-1,:) + fishPos(badInds+1,:));
    dS_all = [0; S(fishPos)];
    dS_all(trlFrames)=0;
    badInds = find(dS_all>dS_thr);
    badInds(badInds==1)=[];
    badInds(badInds >= (length(dS_all)-1))=[];    
end

for indNum = 1:length(badInds)
    ind = badInds(indNum);
    for seg = 1:length(midlineInds{ind})
        midlineInds{ind}{seg} = 0.5*(midlineInds{ind-1}{seg} + midlineInds{ind+1}{seg});
    end
end
orientation = GetFishOrientationFromMidlineInds(midlineInds,imgDims(1:2),'s');
orientation = orientation';
motionInfo = GetMotionInfo(fishPos,orientation,imgDims(1));
curv = motionInfo.curv;

%%  Swim episodes
stimFrames = (0:nTrls-1)*nFramesInTrl + (preStimPer*fps);
stimVec = zeros(size(fishPos,1),1);
stimVec(stimFrames)=1;
curv = motionInfo.curv;
ker = gausswin(round(100e-3*fps)); ker = ker/sum(ker);
try
    swim.smooth = Standardize(conv2(abs(chebfilt(curv(:,end),1/fps,20,'high')),ker(:),'same'));
catch
    swim.smooth = Standardize(conv2(abs(curv(:,end)),ker(:),'same'));
end
% swim.smooth = Standardize(conv2(chebfilt(curv(:,end),1/fps,20,'high').^2,ker(:),'same'));
[maxtab,~] = peakdet(swim.smooth,peakDetThr);
swim.inds = maxtab(:,1);
swim.times = swim.inds/fps;
swim.bool = zeros(size(swim.smooth));
swim.bool(swim.inds)= 1;
disp([num2str(sum(swim.bool)) ' swim episodes detected'])


%% Trializing some variables
if mod(nTrls,1) ~= 0
    error('Check value of nFramesInTrl, does not evenly divide into total # of frames')
else
    disp([num2str(nTrls) ' trials detected'])
end
SortTrls = @(var,nFramesInTrl,nTrls) permute(reshape(var,[nFramesInTrl,nTrls,size(var,2)]),[1 3 2]);
SepTapAndFlashTrls = @(var,tapTrls,flashTrls) cat(4, var(:,:,tapTrls), var(:,:,flashTrls));
curv_trl = SepTapAndFlashTrls(SortTrls(curv,nFramesInTrl,nTrls),tapTrls,flashTrls);
if size(orientation,1) < size(orientation,2)
    orientation = orientation';
    disp('Transposed orientation!')
end
or_trl =  SepTapAndFlashTrls(SortTrls(orientation,nFramesInTrl,nTrls),tapTrls,flashTrls);
time_trl = ((0:nFramesInTrl-1)/fps);


%% Separating tap and dark flash trls
% tap = curv_trl(:,:,tapTrls);
% dark = curv_trl(:,:,flashTrls);
% curv_trl = cat(4,tap,dark);

%% Finding left and right curve peaks
if fps > 100
    blah = zscore(chebfilt(curv(:,3),1/fps,15,'high'));
    [lPks,rPks] = peakdet(blah,1);
end

%% Plotting tap and dark trl curvatures
postStimPer = (nFramesInTrl*(1/500))-0.1;
for stim = 1:size(curv_trl,4)
    figure('Name',['Stim ' num2str(stim)])
    count = 0;
    for trl = 1:size(curv_trl,3)
        yShift = (trl-1)*yOff(stim);
        if mod(count,2)==0
            clr = 'r';
        else
            clr = 'g';
        end
        plot(time_trl*1000,curv_trl(:,seg,trl,stim)-yShift, 'color',clr)
        hold on
        drawnow
        count = count + 1;
    end
    box off
    set(gca,'tickdir','out','color','k')
    xlabel('Time (ms)')
    ylim([-yShift-yOff(stim) yOff(stim)])
    xlim([-preStimPer*1000 postStimPer*1000])
end

%% Raster plot of trialized swim episodes
spikeAmp = 0.9;
swim.bool_trl = reshape(swim.bool,[nFramesInTrl,nTrls])';
% stimVec_trl  = reshape(stimVec,[nFramesInTrl,nTrls])';
figure('Name','Raster plot of trialized swim episodes')
[r,c] = find(swim.bool_trl);
plot(([c c]'*(1/fps)-(mult*preStimPer)), [r-spikeAmp/2 r+spikeAmp/2]','g-');
set(gca,'color','k','tickdir','out','ydir','reverse','ytick',1:2:nTrls)
box off
ylim([-inf inf])
xlim([-inf inf])
ylabel('Trl #')
xlabel('Time (sec)')

%% Traj plots
fishPos_trl = permute(reshape(fishPos,[nFramesInTrl, nTrls, size(fishPos,2)]),[1 3 2]);
tap = fishPos_trl(:,:,tapTrls);
df = fishPos_trl(:,:,flashTrls);
fishPos_trl = cat(4,tap,df);
clrMap = hsv(nFramesInTrl*2);
clrMap = clrMap(nFramesInTrl+1:end,:);
% clrMap = hsv(nFramesInTrl);
% edgeInds = GetArenaEdge(ref);
% close
for stim = 1:size(fishPos_trl,4)
    figure('Name',['Stim ' num2str(stim)])
    %     imagesc(ref),axis image
    %     plot(edgeInds(:,2),edgeInds(:,1),'w--'),
    %     axis image
    set(gca,'color','k')
    hold on
    count = 0;
    for trl = 1:size(fishPos_trl,3)
        for pt = 1:2:nFramesInTrl
            plot(fishPos_trl(pt,1,trl,stim),fishPos_trl(pt,2,trl,stim),'.','color',clrMap(pt,:),'markersize',5)
        end
        drawnow
        count = count + 1;
    end
    box off
    set(gca,'tickdir','out','color','k')
    axis image
    set(gca,'xtick',[],'ytick',[])
    ylim([1 imgDims(1)])
    xlim([-inf imgDims(2)])
end

%% Traj plots - adjusted
if size(orientation,2) == size(fishPos,1)
    orientation = orientation';
end
stimFrameInTrl = round(preStimPer*fps);

for stim = 1:size(fishPos_trl,4)
    disp(['Plotting for stim # ' num2str(stim)])
    figure('Name',['Stim ' num2str(stim)])
    %     imagesc(ref),axis image
    %     plot(edgeInds(:,2),edgeInds(:,1),'w--'),
    %     axis image
    set(gca,'color','k')
    hold on
    count = 0;
    plot([0 0], [-100 100],'g--')
    plot([-100 100],[0 0],'g--')
    for trl = 1:size(fishPos_trl,3)
        xy_start = [fishPos_trl(stimFrameInTrl,1,trl,stim),fishPos_trl(stimFrameInTrl,2,trl,stim)];
        or = or_trl(stimFrameInTrl,1,trl,stim);
        %         or = mod((360-or)+90,360);
        or = (360-or) + 90;
        %         cla
        for pt = 1:nFramesInTrl
            %             xy = RotateVecIn2D([fishPos_trl(pt,1,trl,stim)-xy_start(1),fishPos_trl(pt,2,trl,stim)-xy_start(2)],or);
            xy = RotateTraj([fishPos_trl(pt,1,trl,stim)-xy_start(1),fishPos_trl(pt,2,trl,stim)-xy_start(2)],or);
            plot(xy(1),xy(2),'.','color',clrMap(pt,:),'markersize',5)
        end
        count = count + 1;
    end
    drawnow
    shg
    box off
    set(gca,'tickdir','out','color','k')
    axis image
    set(gca,'xtick',[],'ytick',[])
    ylim([-100 100])
    xlim([-100 100])
end

%% Plot left, right curvature amplitudes
if fps >= 100
    ker = gausswin(round(0.02*fps)); ker = ker/sum(ker);
    curv_flt = conv2(chebfilt(curv(:,seg),1/fps,15,'high'),ker(:),'same');
    [lPks,rPks] = peakdet(zscore(curv_flt),1);
    lPks(:,2) = curv_flt(lPks(:,1));
    rPks(:,2) = curv_flt(rPks(:,1));
    binC = (0:0.5:11)*std(curv_flt);
    [lCount,vals] = hist(lPks(:,2),binC);
    [rCount,~] = hist(-rPks(:,2),binC);
    vals_interp = linspace(vals(1),vals(end),length(vals)*5);
    lCount_interp = interp1(vals,lCount,vals_interp,'spline');
    rCount_interp = interp1(vals,rCount,vals_interp,'spline');
    N = sum(lCount_interp) + sum(rCount_interp);
    figure('Name','Left and right curv amp hist')
    %     plot(vals_interp,lCount_interp/N,'.-')
    %     plot(vals_interp,lCount_interp/sum(lCount_interp),'.-')
    plot(vals_interp,lCount_interp/max(lCount_interp),'.-')
    hold on
    %     plot(vals_interp,rCount_interp/N,'r.-')
    %      plot(vals_interp,rCount_interp/sum(rCount_interp),'r.-')
    plot(vals_interp,rCount_interp/max(rCount_interp),'r.-')
    legend('Left','Right');
    xlim([-inf inf])
    ylim([-inf inf])
    xlabel('Curv amp')
    ylabel('Prob')
    box off
    set(gca,'tickdir','out')
    shg
end

%% Plot left, right curv velocities
if fps >=100
    time_interp = linspace(time(1),time(end),length(time)*2);
    curv_flt_interp = interp1(time,curv_flt,time_interp,'spline');
    dCurv = diff(curv_flt_interp)./diff(time_interp);
    [lPks,rPks] = peakdet(zscore(dCurv),1);
    lPks(:,2) = dCurv(lPks(:,1));
    rPks(:,2) = dCurv(rPks(:,1));
    binC = (0:0.5:11)*std(dCurv);
    [lCount,vals] = hist(lPks(:,2),binC);
    [rCount,~] = hist(-rPks(:,2),binC);
    vals_interp = linspace(vals(1),vals(end),length(vals)*5);
    lCount_interp = interp1(vals,lCount,vals_interp,'spline');
    rCount_interp = interp1(vals,rCount,vals_interp,'spline');
    N = sum(lCount_interp) + sum(rCount_interp);
    figure('Name','Left and right angular velocity hist')
    %     plot(vals_interp,lCount_interp/N,'.-')
    plot(vals_interp,lCount_interp/max(lCount_interp),'.-')
    hold on
    %     plot(vals_interp,rCount_interp/N,'r.-')
    plot(vals_interp,rCount_interp/max(rCount_interp),'r.-')
    legend('Left','Right');
    xlim([-inf inf])
    ylim([-inf inf])
    xlabel('Curv velocities')
    ylabel('Prob')
    box off
    set(gca,'tickdir','out')
    shg
end

%% Plot swim activity based on displacement
time = (0:size(fishPos,1)-1)/fps;
blah = dS_all;
blah(trlFrames)=0;
blah = chebfilt(blah,1/fps,10,'low');
blah(abs(blah)<=1.5)=0;
blah(blah<0)=0;
plot(time,blah,'g')
hold on
stem(time(trlFrames),ones(size(trlFrames))*60,'r','marker','none')
set(gca,'color','k')
shg


