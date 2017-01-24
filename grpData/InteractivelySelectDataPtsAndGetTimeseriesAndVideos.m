

%% Read in grpData from here 'S:\Avinash\Ablations and behavior\GrpData'
% currDir = cd;
% cd('S:\Avinash\Ablations and behavior\GrpData')
% disp('Reading grpData...')
% grpData = OpenMatFile();
% grpData = grpData.grpData;
%
if isfield(grpData,'data_xls')
    data_xls = grpData.data_xls;
else
    data_xls = UnrolledDataFrmGrpData(grpData);
end
data_xls.bendPer = data_xls.bendPer + rand(size(data_xls.bendPer))*2-1;
% cd(currDir)

%% Select trials to play from interactive data plots
grps = unique(data_xls.GrpName);
disp('Select Group Index: ');
disp(grps)
grpIdx = input('Select group index, e.g., 1 : ');

stims = unique(data_xls.StimType);
disp('Select Stim Type:')
disp(stims)
stimIdx = input('Select stim index, e.g., 1 : ');

grpInds = strcmpi(data_xls.GrpName,grps{grpIdx});
stimInds = strcmpi(data_xls.StimType,stims{stimIdx});
bendOneInds = data_xls.BendNum == 1;
nonNanInds =  1-isnan(data_xls.onset);
ctrlInds = strcmpi(data_xls.Treatment,'ctrl');
ablInds = strcmpi(data_xls.Treatment,'abl');

ctrlInds = find(grpInds(:).*stimInds(:).*bendOneInds(:).*nonNanInds(:).*ctrlInds(:));
ablInds = find(grpInds(:).*stimInds(:).*bendOneInds(:).*nonNanInds(:).*ablInds(:));


figure('Name','Select data points...')
subplot(2,1,1)
x = data_xls.bendPer(ctrlInds);
y = abs(data_xls.bendAmp(ctrlInds));
plot(x,y,'.')
hold on
x = data_xls.bendPer(ablInds);
y = abs(data_xls.bendAmp(ablInds));
plot(x,y,'ro')
legend('Ctrl','Abl');
box off
ylabel('Bend amp')
xlabel('Bend per')

subplot(2,1,2)
x = data_xls.onset(ctrlInds)+ data_xls.bendPer(ctrlInds)/3;
y = abs(data_xls.bendAmp(ctrlInds));
plot(x,y,'.')
hold on
x = data_xls.onset(ablInds)+ data_xls.bendPer(ablInds)/3;
y = abs(data_xls.bendAmp(ablInds));
plot(x,y,'ro')
legend('Ctrl','Abl');
box off
ylabel('Bend amp')
xlabel('Onset')


% Matching interactively chosen data points to corresponding indices in data_xls
selInds = cell(2,1);
% --- Control
subplot(2,1,1)
title('Choose control data points from one of the 2 figures')
[xx,yy] = ginput_plot();

chosenVar = input('Did you chose from Period [1] or Onset [2]? ');
if chosenVar ==1
    x = data_xls.bendPer(ctrlInds);
else
    x = data_xls.onset(ctrlInds);
end
y = abs(data_xls.bendAmp(ctrlInds));
V = [x(:), y(:)];

inds = nan(length(xx),1);
for jj = 1:length(xx)
    d = sqrt(sum((V - repmat([xx(jj),yy(jj)],length(V),1)).^2,2));
    [~, ind] = min(abs(d));
    inds(jj) = ctrlInds(ind);
end

selInds{1} = inds;

% --- Ablated
subplot(2,1,1)
title('Choose ablated data points from one of the 2 figures')
[xx,yy] = ginput_plot();

chosenVar = input('Did you chose from Period [1] or Onset [2]? ');
if chosenVar ==1
    x = data_xls.bendPer(ablInds);
else
    x = data_xls.onset(ablInds);
end
y = abs(data_xls.bendAmp(ablInds));
V = [x(:), y(:)];

inds = nan(length(xx),1);
for jj = 1:length(xx)
    d = sqrt(sum((V - repmat([xx(jj),yy(jj)],length(V),1)).^2,2));
    [~, ind] = min(abs(d));
    inds(jj) = ablInds(ind);
end

selInds{2} = inds;


%% Getting the tail curves, the image frames, and other relevant data
tic
nFramesInTrl = 750;
fps = 500;
data = cell(size(selInds));
for jj = 1:2
    if jj ==1
        disp('Ctrl trls...')
    else
        disp('Abl trls...')
    end
    temp = cell(length(selInds{jj}),1);
    for kk = 1:length(selInds{jj})
        disp(['Trl # ' num2str(kk)])
        blah = struct;
        ind = selInds{jj}(kk);
        grp = data_xls.GrpName{ind};
        stim = data_xls.StimType{ind};
        trtmnt = data_xls.Treatment{ind};
        fishNum = data_xls.FishNum(ind);
        trlNum = data_xls.TrlNum(ind);
        onset = data_xls.onset(ind);
        amp = data_xls.bendAmp(ind);
        per = data_xls.bendAmp(ind);
        trlInds = ...
            (trlNum-1)*nFramesInTrl+1:(trlNum-1)*nFramesInTrl+nFramesInTrl;
        pd = grpData.(grp).(trtmnt).(stim).procData{fishNum};
        blah.path = pd.Properties.Source;
        blah.im = pd.IM_crop(:,:,trlInds);
        curv = GetTailTangents(pd.tailCurv(:,:,trlInds),5);
        curv = curv(end,:);
        blah.curv_unflipped = curv;
        if amp < 0
            curv = -curv; % For aligning
        end
        blah.curv = curv;
        blah.time = (0:length(trlInds)-1)*(1/fps)*1000 - onset;       
        blah.trlInds = trlInds;
        temp{kk} = blah;
    end
    data{jj} = temp;
end
toc


%% Plot traces for chosen data points, onset aligned
nFramesInTrl = 750;
fps = 500;
xl = [50 400];
alphaVal = 0.5;
% yShift = 300;
figure('Name','Example traces')
% clrs = {'b','r'};
clrs = {[0 0 1], [1 0 0]};
m = nan(length(data),1);
for jj = 1:length(data)
    m(jj) = length(data{jj});
end
m = max(m);
for jj = 1:length(data)
    for kk = 1:length(data{jj})
        x = data{jj}{kk}.time;
        y = data{jj}{kk}.curv;
        %         y = y-mean(y)-((kk-1)*yShift);
        y = y-mean(y);
        %         plot(x,y,'color',clrs{jj})
        subplot(m,1,kk)
        patchline(x,y,'edgecolor',clrs{jj},'edgealpha',alphaVal,'linewidth',1.5)
        hold on
        alpha(alphaVal)
        box off
        set(gca,'tickdir','out','ytick',[0 200])
        xlim(xl)
        ylim([-inf inf])
        if kk == 2
            ylabel('Total body curvature (deg)')
        end
    end
end


%% Plot traces for chosen data points, stim aligned
xl = [-50 600];
alphaVal = 0.5;
% yShift = 300;
figure('Name','Example traces')
% clrs = {'b','r'};
clrs = {[0 0 1], [1 0 0]};
m = nan(length(data),1);
for jj = 1:length(data)
    m(jj) = length(data{jj});
end
m = max(m);
for jj = 1:length(data)
    for kk = 1:length(data{jj})
        %         x = data{jj}{kk}.time;
        x = (0:length(data{jj}{kk}.time)-1)*(1/fps)*1000 -100;
        y = data{jj}{kk}.curv;
        %         y = y-mean(y)-((kk-1)*yShift);
        y = y-mean(y);
        %         plot(x,y,'color',clrs{jj})
        subplot(m,1,kk)
        patchline(x,y,'edgecolor',clrs{jj},'edgealpha',alphaVal,'linewidth',1.5)
        hold on
        alpha(alphaVal)
        box off
        set(gca,'tickdir','out','ytick',[0 200])
        xlim(xl)
        ylim([-inf inf])
        if kk == 2
            ylabel('Total body curvature (deg)')
        end
    end
end
xlabel('Time (ms)')

%% Videos
fh = figure('Name','Example trial for ctrl and abl');
pauseDur = 0;
trls = [9 1]; % [trl for ctrl fish, trl for ablated fish] - these are played side-by-side
frameInds = [1:250];
dispTraces = 1; % Logical for displaying curvature timeseries underneath video frames
writeVideo = 0; % For writing videos
saveDir = 'S:\Avinash\Ablations and behavior\GrpData\Figs_20170113\20170119\Example traces and videos\M-homolog';

if isempty(frameInds)
    frameInds = 1:nFramesInTrl;
end
if numel(trls)==1
    trls = repmat(trls,1,2);
end

if writeVideo
    vidObj = VideoWriter(fullfile(saveDir,'ExampleTrlVideos.avi'));
    open(vidObj)
    disp('Writing to video...')
end

if dispTraces
    ax = cell(4,1);
    ax{1} = [0.5 0.7 0 0.3];
    ax{2} = [0.5 0.7 0.5 0.3];
    ax{3} = [0.5 0.3 0 0];
    ax{4} = [0.5 0.3 0.5 0];
    axH = CreateSubaxes(fh,ax{1},ax{2},ax{3},ax{4});
    colormap(gray)
    for ff = frameInds(:)'
        im1 = data{1}{trls(1)}.im(:,:,ff);
        im2 = data{2}{trls(2)}.im(:,:,ff);
        y1 = data{1}{trls(1)}.curv;
        y2 = data{2}{trls(2)}.curv;
        
        axes(axH(1))
        imagesc(im1)
        axis off
        title(['Trls: ' num2str(trls) ', frame: ' num2str(ff) '/' num2str(nFramesInTrl)])
        
        axes(axH(2))
        imagesc(im2)
        axis off
        
        inds = ff-100:ff;
        inds(inds<=0)=[];
        axes(axH(3))
        plot(inds, y1(inds),'b-')
        xlim([ff-100 ff])
        ylim([-100 250])
        box off
        set(gca,'tickdir','out','ytick',[0 200]);
        
        axes(axH(4))
        plot(inds,y2(inds),'r-')
        xlim([ff-100 ff])
        ylim([-100 250])
        box off
        set(gca,'tickdir','out','ytick',[]);
        drawnow
        
        if writeVideo
            currFrame = getframe(fh);
            vidObj.writeVideo(currFrame)
        end
        
        if isempty(pauseDur)
            pause()
        else
            pause(pauseDur)
        end
    end
else
    for ff = frameInds(:)'
        im1 = data{1}{trls(1)}.im(:,:,ff);
        im2 = data{2}{trls(2)}.im(:,:,ff);
        fh;
        imshowpair(im1,im2,'montage','Scaling','joint')
        axis off
        title(['Trls: ' num2str(trls) ', frame: ' num2str(ff) '/' num2str(nFramesInTrl)])
        drawnow
        
        if writeVideo
            currFrame = getframe(fh);
            vidObj.writeVideo(currFrame)
        end
        
        if isempty(pauseDur)
            pause()
        else
            pause(pauseDur)
        end
    end
end

if writeVideo
    close(vidObj)
end
