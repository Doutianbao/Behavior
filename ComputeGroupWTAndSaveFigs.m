
%% Setting inputs
xlsPath = 'S:\Avinash\Ablations and behavior\Ablation data summary.xlsx';
outDir = 'S:\Avinash\Ablations and behavior\GrpData';
xLim_tap = [-50 650];
xLim_dark = [-100 300];
onsetAlign = 1;
plotOrNot = 0;
sigmaXY = 100; %(default = NaN);
cLim = [0 80];
freqRange_dark = [5 35];

%% Ctrl, tap
disp('Getting ctrl, tap data...')
[paths,paramVals, filtInds] = GetFilteredPathsFromXLS(xlsPath);
for pp  = 1:length(paths)
    [path_curr, file_curr] = fileparts(paths{pp});
    if ~strcmpi(file_curr,'proc')
        paths{pp} = fullfile(paths{pp},'proc');
    end
end
% out = AnalyzeFreeSwims_nCycles_batch(paths, 'xLim',[0 750],'paramList',...
%     {'bodyAmp','headAmp'});
[~, procData] = GetFishWaves_group(paths, 'saveToProc',1,'xLim',xLim_tap,...
    'onsetAlign',onsetAlign,'sigmaXY',sigmaXY,'plotOrNot',plotOrNot,'cLim',cLim,'traceType','curv');
data.ctrl.vib.procData = procData;
data.ablationType = paramVals.AblationType;

for fish = 1:length(procData);
    [outDir, ~] = fileparts(procData{fish}.Properties.Source);
    subDirName = 'WT';
    outDir = fullfile(outDir,subDirName);
    if ~exist(outDir,'dir')
        mkdir(outDir)
    end
    close all
    PlotWTs(procData{fish}.W);
    h = get(0,'children');
    disp(['Saving figures for fish # ' num2str(fish)]);
    for fig = 1:length(h)
        saveas(h(fig), fullfile(outDir,['Fig_' sprintf('%.2d',h(fig))]))
    end
    dispNext = input('Display next set of figures?(y/n): ','s');
    if strcmpi(dispNext,'n')
        break
    end
end


%% Abl tap
disp('Getting abl, tap data...')
[paths,~] = GetFilteredPathsFromXLS(xlsPath,filtInds);
for pp  = 1:length(paths)
    [path_curr, file_curr] = fileparts(paths{pp});
    if ~strcmpi(file_curr,'proc')
        paths{pp} = fullfile(paths{pp},'proc');
    end
end
% out = AnalyzeFreeSwims_nCycles_batch(paths, 'xLim',[0 750],'paramList',...
%     {'bodyAmp','headAmp'});
[~, procData] = GetFishWaves_group(paths(1:end-1), 'saveToProc',1,'xLim',xLim_tap,...
    'onsetAlign',onsetAlign,'sigmaXY',sigmaXY,'plotOrNot',plotOrNot,'cLim',cLim,'traceType','curv');
data.abl.vib.procData = procData;
for fish = 1:length(procData);
    [outDir, ~] = fileparts(procData{fish}.Properties.Source);
    subDirName = 'WT';
    outDir = fullfile(outDir,subDirName);
    if ~exist(outDir,'dir')
        mkdir(outDir)
    end
    close all
    PlotWTs(procData{fish}.W);
    h = get(0,'children');
    disp(['Saving figures for fish # ' num2str(fish)]);
    for fig = 1:length(h)
        saveas(h(fig), fullfile(outDir,['Fig_' sprintf('%.2d',h(fig))]))
    end
    dispNext = input('Display next set of figures?(y/n): ','s');
    if strcmpi(dispNext,'n')
        break
    end
end

%% Ctrl, dark flash
disp('Getting ctrl, dark flash data...')
[paths,~] = GetFilteredPathsFromXLS(xlsPath,filtInds);
for pp  = 1:length(paths)
    [path_curr, file_curr] = fileparts(paths{pp});
    if ~strcmpi(file_curr,'proc')
        paths{pp} = fullfile(paths{pp},'proc');
    end
end
% out = AnalyzeFreeSwims_nCycles_batch(paths, 'xLim',[0 1500],'paramList',...
%     {'bodyAmp'});

[~, procData] = GetFishWaves_group(paths, 'saveToProc',1,'xLim',xLim_dark,...
    'onsetAlign',onsetAlign,'sigmaXY',sigmaXY,'plotOrNot',plotOrNot,'cLim',[0 50],'freqRange',freqRange_dark,'traceType','curv');

data.ctrl.dark.procData = procData;
for fish = 1:length(procData);
    [outDir, ~] = fileparts(procData{fish}.Properties.Source);
    subDirName = 'WT';
    outDir = fullfile(outDir,subDirName);
    if ~exist(outDir,'dir')
        mkdir(outDir)
    end
    close all
    PlotWTs(procData{fish}.W);
    h = get(0,'children');
    disp(['Saving figures for fish # ' num2str(fish)]);
    for fig = 1:length(h)
        saveas(h(fig), fullfile(outDir,['Fig_' sprintf('%.2d',h(fig))]))
    end
    dispNext = input('Display next set of figures?(y/n): ','s');
    if strcmpi(dispNext,'n')
        break
    end
end

%% Abl, dark flash
disp('Getting abl, dark flash data...')
[paths,~] = GetFilteredPathsFromXLS(xlsPath,filtInds);
for pp  = 1:length(paths)
    [path_curr, file_curr] = fileparts(paths{pp});
    if ~strcmpi(file_curr,'proc')
        paths{pp} = fullfile(paths{pp},'proc');
    end
end

% out = AnalyzeFreeSwims_nCycles_batch(paths(2:end), 'xLim',[0 1500],'paramList',...
%     {'bodyAmp'});

% out = AnalyzeFreeSwims_nCycles_batch(paths, 'xLim',[0 1500],'paramList',...
%     {'bodyAmp','headAmp'});
[~, procData] = GetFishWaves_group(paths([1:6 8:end]), 'saveToProc',1,'xLim',xLim_dark,...
    'onsetAlign',onsetAlign,'sigmaXY',sigmaXY,'plotOrNot',plotOrNot,'cLim',[0 50],'freqRange',freqRange_dark,'traceType','curv');

data.abl.dark.procData = procData;
for fish = 1:length(procData);
    [outDir, ~] = fileparts(procData{fish}.Properties.Source);
    subDirName = 'WT';
    outDir = fullfile(outDir,subDirName);
    if ~exist(outDir,'dir')
        mkdir(outDir)
    end
    close all
    PlotWTs(procData{fish}.W);
    h = get(0,'children');
    disp(['Saving figures for fish # ' num2str(fish)]);
    for fig = 1:length(h)
        saveas(h(fig), fullfile(outDir,['Fig_' sprintf('%.2d',h(fig))]))
    end
    dispNext = input('Display next set of figures?(y/n): ','s');
    if strcmpi(dispNext,'n')
        break
    end
end

%% WT averages across fish - Head & Tail Orientations
tic
W_all = cell(2,1);
%--- Ctrl ---
disp('Avg WTs from all control fish...')
procData = data.ctrl.vib.procData;
nFish = length(procData);
for fishNum = 1:nFish
    disp(['Fish # ' num2str(fishNum)])
    W = procData{fishNum}.W;
    W = cat(1,W.head.avg,W.tail.avg);
    if fishNum ==1
        W_all{1} = W;
    else
        W_all{1} = cat(3,W_all{1},W);
    end
end
data.ctrl.vib.mean = mean(W_all{1},3);
data.ctrl.vib.std = std(W_all{1},[],3);
data.ctrl.vib.cv = data.ctrl.vib.std./(data.ctrl.vib.mean + 1);
corrVec = nan(1,nFish);
for fishNum = 1:nFish
    corrVec(fishNum) = corr2(data.ctrl.vib.mean,W_all{1}(:,:,fishNum));
end
data.ctrl.vib.corrVec = corrVec;

%--- Abl ---
disp('Avg WTs from all ablated fish...')
procData = data.abl.vib;
nFish = length(procData);
for fishNum = 1:nFish
    disp(['Fish # ' num2str(fishNum)])
    W = procData{fishNum}.W;
    W = cat(1,W.head.avg,W.tail.avg);
    if fishNum ==1
        W_all{2} = W;
    else
        W_all{2} = cat(3,W_all{2},W);
    end
end

data.abl.vib.mean = mean(W_all{2},3);
data.abl.vib.std = std(W_all{2},[],3);
data.abl.vib.cv = data.abl.vib.std./(data.abl.vib.mean + 1);
corrVec = nan(1,nFish);
for fishNum = 1:nFish
    corrVec(fishNum) = corr2(data.abl.vib.mean,W_all{2}(:,:,fishNum));
end
data.abl.vib.corrVec = corrVec;
toc


%% WT averages across fish - Vib stim and Curvatures
tic
W_all = cell(2,1);

%--- Ctrl ---
disp('Avg WTs from all control fish...')
procData = data.ctrl.vib;
nFish = length(procData);
for fishNum = 1:nFish
    disp(['Fish # ' num2str(fishNum)])
    W = procData{fishNum}.W;
    W = W.curv.avg;
    if fishNum ==1
        W_all{1} = W;
    else
        W_all{1} = cat(3,W_all{1},W);
    end
end
data.ctrl.vib.mean = mean(W_all{1},3);
data.ctrl.vib.std = std(W_all{1},[],3);
data.ctrl.vib.cv = data.ctrl.vib.std./(data.ctrl.vib.mean + 1);
corrVec = nan(1,nFish);
for fishNum = 1:nFish
    corrVec(fishNum) = corr2(data.ctrl.vib.mean,W_all{1}(:,:,fishNum));
end
data.ctrl.vib.corrVec = corrVec;

%--- Abl ---
disp('Avg WTs from all ablated fish...')
procData = data.abl.vib;
nFish = length(procData);
for fishNum = 1:nFish
    disp(['Fish # ' num2str(fishNum)])
    W = procData{fishNum}.W;
    W = W.curv.avg;
    if fishNum ==1
        W_all{2} = W;
    else
        W_all{2} = cat(3,W_all{2},W);
    end
end

data.abl.vib.mean = mean(W_all{2},3);
data.abl.vib.std = std(W_all{2},[],3);
data.abl.vib.cv = data.abl.vib.std./(data.abl.vib.mean + 1);
corrVec = nan(1,nFish);
for fishNum = 1:nFish
    corrVec(fishNum) = corr2(data.abl.vib.mean,W_all{2}(:,:,fishNum));
end
data.abl.vib.corrVec = corrVec;
toc

%% The average WT across fish from average WTs for fish across trials - for dark flash stim and curvatures
tic
W_all = cell(2,1);

%--- Ctrl ---
disp('Avg WTs from all control fish...')
procData = data.ctrl.dark;
nFish = length(procData);
for fishNum = 1:nFish
    disp(['Fish # ' num2str(fishNum)])
    W = procData{fishNum}.W;
    W = W.curv.avg;
    if fishNum ==1
        W_all{1} = W;
    else
        W_all{1} = cat(3,W_all{1},W);
    end
end
data.ctrl.dark.mean = mean(W_all{1},3);
data.ctrl.dark.std = std(W_all{1},[],3);
data.ctrl.dark.std = data.ctrl.dark.std./(data.ctrl.dark.mean + 1);
corrVec = nan(1,nFish);
for fishNum = 1:nFish
    corrVec(fishNum) = corr2(data.ctrl.dark.mean,W_all{1}(:,:,fishNum));
end
data.ctrl.dark.corrVec = corrVec;

%--- Abl ---
disp('Avg WTs from all ablated fish...')
procData = data.abl.dark;
nFish = length(procData);
for fishNum = 1:nFish
    disp(['Fish # ' num2str(fishNum)])
    W = procData{fishNum}.W;
    W = W.curv.avg;
    if fishNum ==1
        W_all{2} = W;
    else
        W_all{2} = cat(3,W_all{2},W);
    end
end

data.abl.dark.mean = mean(W_all{2},3);
data.abl.dark.std = std(W_all{2},[],3);
data.abl.dark.std = data.abl.dark.std./(data.abl.dark.mean + 1);
corrVec = nan(1,nFish);
for fishNum = 1:nFish
    corrVec(fishNum) = corr2(data.abl.dark.mean,W_all{2}(:,:,fishNum));
end
data.abl.dark.corrVec = corrVec;
toc


%% Saving data
saveOrNot = input('Save group data?(y/n): ', 's');
if strcmpi(saveOrNot,'y')
    disp('Saving data...')
    tic
    timeStamp = datestr(now,30);
    fName = ['procData_grp_' data.ablationType timeStamp '.mat'];
    save(fullfile(outDir,fName),'data');
    toc
end

%% Plot Figs - for head and tail orientations

% -- Fig 1: Head and Tail WT for ctrl and abl ---
ax = {};
ax{1} = [0.5 0.5 0 0.5];
ax{2} = [0.5 0.5 0 0];
ax{3} = [0.5 0.5 0.5 0.5];
ax{4} = [0.5 0.5 0.5 0];
fh = figure('Name','Avg WT for ctrl and abl');
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3},ax{4});
W = data.ctrl.vib{1}.W;
freq = W.freq;
time = W.time;

cl = W.cLim;
cl(2) = cl(2);
cl_cv = [0 10];
midRow = size(data.ctrl.vib.mean,1)/2;

%--- Mean head for ctrl --
axes(axH(1))
imagesc(time,freq,abs(data.ctrl.vib.mean(1:midRow,:)).^1)
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
title('CONTROL')
ylabel('Freq (Hz)')

%--- Mean tail for ctrl --
axes(axH(2))
imagesc(time,freq,abs(data.ctrl.vib.mean(midRow+1:end,:)).^1)
set(gca,'clim',cl,'ydir','normal','tickdir','out')
box off
xlabel('Time (ms)')
ylabel('Freq (Hz)')

%--- Mean head for abl --
axes(axH(3))
imagesc(time,freq,abs(data.abl.vib.mean(1:midRow,:)).^1)
set(gca,'clim',cl,'ydir','normal','tickdir','out','ytick',[],'xtick',[],'yaxislocation','right')
box off
ylabel('HEAD')
title('ABLATED')

%--- CV tail --
axes(axH(4))
imagesc(time,freq,abs(data.abl.vib.mean(midRow+1:end,:)).^1)
set(gca,'clim',cl,'ydir','normal','tickdir','out','ytick',[],'yaxislocation','right')
box off
xlabel('Time (ms)')
ylabel('TAIL')


%--- Fig 2: Difference WT maps ---
fh = figure('Name','Difference WT maps')
ax = {};
ax{1} = [1 0.5 0 0.5];
ax{2} = [1 0.5 0 0];
axH = CreateSubaxes(fh,ax{1},ax{2});

axes(axH(1))
W_diff = data.abl.vib.mean(1:midRow,:) - data.ctrl.vib.mean(1:midRow,:);
% W_diff(abs(W_diff)<5) = 0;
imagesc(time,freq,W_diff);
set(gca,'ydir','normal','tickdir','out','clim',[-50 50])
box off
title('W_{abl} - W_{ctrl}')
ylabel({'HEAD';'Freq (Hz)'})
ch = colorbar;
set(ch,'location','North')

axes(axH(2))
W_diff = data.abl.vib.mean(midRow+1:end,:) - data.ctrl.vib.mean(midRow+1:end,:);
% W_diff(abs(W_diff)<5) = 0;
imagesc(time,freq,W_diff);
set(gca,'ydir','normal','tickdir','out','clim',[-50 50])
box off
xlabel('Time (ms)')
ylabel({'TAIL';'Freq (Hz)'})


% --- Fig 3: Variablity plot for ctrl and abl ---
figure('Name','Reponse variability across fish')
plot(sort(data.ctrl.vib.corrVec),'.-')
hold on
plot(sort(data.abl.vib.corrVec),'ro-')
set(gca,'tickdir','out'), box off
xlabel('Fish #')
ylabel('Response correlation to average response')
lh = legend('Control','Ablated');
set(lh,'location','best')
title('Response variability across fish for control and ablated groups')


%% Plot Figs - for body curvatures

% -- Fig 1: Curvature WT for ctrl and abl ---
ax = {};
ax{1} = [0.5 0.5 0 0.5];
ax{2} = [0.5 0.5 0 0];
ax{3} = [0.5 0.5 0.5 0.5];
ax{4} = [0.5 0.5 0.5 0];
fh = figure('Name','Avg WT for ctrl and abl');
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3},ax{4});
W = data.ctrl.vib{1}.W;
freq = W.freq;
time = W.time;

cl = W.cLim;
cl(2) = cl(2);
cl_cv = [0 10];
midRow = size(data.ctrl.vib.mean,1)/2;

%--- Mean head for ctrl --
axes(axH(1))
imagesc(time,freq,abs(data.ctrl.vib.mean(1:midRow,:)).^1)
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
title('CONTROL')
ylabel('Freq (Hz)')

%--- Mean tail for ctrl --
axes(axH(2))
imagesc(time,freq,abs(data.ctrl.vib.mean(midRow+1:end,:)).^1)
set(gca,'clim',cl,'ydir','normal','tickdir','out')
box off
xlabel('Time (ms)')
ylabel('Freq (Hz)')

%--- Mean head for abl --
axes(axH(3))
imagesc(time,freq,abs(data.abl.vib.mean(1:midRow,:)).^1)
set(gca,'clim',cl,'ydir','normal','tickdir','out','ytick',[],'xtick',[],'yaxislocation','right')
box off
ylabel('HEAD')
title('ABLATED')

%--- CV tail --
axes(axH(4))
imagesc(time,freq,abs(data.abl.vib.mean(midRow+1:end,:)).^1)
set(gca,'clim',cl,'ydir','normal','tickdir','out','ytick',[],'yaxislocation','right')
box off
xlabel('Time (ms)')
ylabel('TAIL')


%--- Fig 2: Difference WT maps ---
fh = figure('Name','Difference WT maps')
ax = {};
ax{1} = [1 0.5 0 0.5];
ax{2} = [1 0.5 0 0];
axH = CreateSubaxes(fh,ax{1},ax{2});

axes(axH(1))
W_diff = data.abl.vib.mean(1:midRow,:) - data.ctrl.vib.mean(1:midRow,:);
% W_diff(abs(W_diff)<5) = 0;
imagesc(time,freq,W_diff);
set(gca,'ydir','normal','tickdir','out','clim',[-50 50])
box off
title('W_{abl} - W_{ctrl}')
ylabel({'HEAD';'Freq (Hz)'})
ch = colorbar;
set(ch,'location','North')

axes(axH(2))
W_diff = data.abl.vib.mean(midRow+1:end,:) - data.ctrl.vib.mean(midRow+1:end,:);
% W_diff(abs(W_diff)<5) = 0;
imagesc(time,freq,W_diff);
set(gca,'ydir','normal','tickdir','out','clim',[-50 50])
box off
xlabel('Time (ms)')
ylabel({'TAIL';'Freq (Hz)'})


% --- Fig 3: Variablity plot for ctrl and abl ---
figure('Name','Reponse variability across fish')
plot(sort(data.ctrl.vib.corrVec),'.-')
hold on
plot(sort(data.abl.vib.corrVec),'ro-')
set(gca,'tickdir','out'), box off
xlabel('Fish #')
ylabel('Response correlation to average response')
lh = legend('Control','Ablated');
set(lh,'location','best')
title('Response variability across fish for control and ablated groups')



