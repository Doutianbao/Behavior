
%% Setting inputs
outDir = 'S:\Avinash\Ablations and behavior\GrpData';
xLim_tap = [-50 650];
xLim_dark = [-100 300];
onsetAlign = 1;
plotOrNot = 0;
sigmaXY = 100; %(default = NaN);
cLim = [0 400];

%% Ctrl, tap
disp('Getting ctrl, tap data...')
[paths,paramVals] = GetFilteredPathsFromXLS();
for pp  = 1:length(paths)
    [path_curr, file_curr] = fileparts(paths{pp});
    if ~strcmpi(file_curr,'proc')
        paths{pp} = fullfile(paths{pp},'proc');
    end
end
% out = AnalyzeFreeSwims_nCycles_batch(paths, 'xLim',[0 750],'paramList',...
%     {'bodyAmp','headAmp'});
[~, procData] = GetFishWaves_group(paths, 'saveToProc',1,'xLim',xLim_tap,...
    'onsetAlign',onsetAlign,'sigmaXY',sigmaXY,'plotOrNot',plotOrNot,'cLim',cLim);
data.ctrl.vib = procData;
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
[paths,~] = GetFilteredPathsFromXLS();
for pp  = 1:length(paths)
    [path_curr, file_curr] = fileparts(paths{pp});
    if ~strcmpi(file_curr,'proc')
        paths{pp} = fullfile(paths{pp},'proc');
    end
end
% out = AnalyzeFreeSwims_nCycles_batch(paths, 'xLim',[0 750],'paramList',...
%     {'bodyAmp','headAmp'});
[~, procData] = GetFishWaves_group(paths, 'saveToProc',1,'xLim',xLim_tap,...
    'onsetAlign',onsetAlign,'sigmaXY',sigmaXY,'plotOrNot',plotOrNot,'cLim',cLim);
data.abl.vib = procData;
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
[paths,~] = GetFilteredPathsFromXLS();
for pp  = 1:length(paths)
    [path_curr, file_curr] = fileparts(paths{pp});
    if ~strcmpi(file_curr,'proc')
        paths{pp} = fullfile(paths{pp},'proc');
    end
end
out = AnalyzeFreeSwims_nCycles_batch(paths(1), 'xLim',[0 1500],'paramList',...
    {'bodyAmp','headAmp'});

[~, procData] = GetFishWaves_group(paths(1), 'saveToProc',1,'xLim',xLim_dark,...
    'onsetAlign',onsetAlign,'sigmaXY',sigmaXY,'plotOrNot',plotOrNot,'cLim',cLim);
data.ctrl.dark = procData;
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
[paths,~] = GetFilteredPathsFromXLS();
for pp  = 1:length(paths)
    [path_curr, file_curr] = fileparts(paths{pp});
    if ~strcmpi(file_curr,'proc')
        paths{pp} = fullfile(paths{pp},'proc');
    end
end
% out = AnalyzeFreeSwims_nCycles_batch(paths, 'xLim',[0 1500],'paramList',...
%     {'bodyAmp','headAmp'});
[~, procData] = GetFishWaves_group(paths, 'saveToProc',1,'xLim',xLim_dark,...
    'onsetAlign',onsetAlign,'sigmaXY',sigmaXY,'plotOrNot',plotOrNot,'cLim',cLim);
data.abl.dark = procData;
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

%% The average WT across fish from average WTs for fish across trials
tic
W_all = cell(2,1);

%--- Ctrl ---
disp('Avg WTs from all control fish...')
var = data.ctrl.vib;
nFish = length(var);
for fishNum = 1:nFish
    disp(['Fish # ' num2str(fishNum)])
    W = var{fishNum}.W;
    W = cat(1,W.head.avg,W.tail.avg);
    if fishNum ==1
        W_all{1} = W;
    else
        W_all{1} = cat(3,W_all{1},W);
    end
end
data.ctrl.vib_mean = mean(W_all{1},3);
data.ctrl.vib_std = std(W_all{1},[],3);
data.ctrl.vib_cv = data.ctrl.vib_std./(data.ctrl.vib_mean + 1);
corrVec = nan(1,nFish);
for fishNum = 1:nFish
    corrVec(fishNum) = corr2(data.ctrl.vib_mean,W_all{1}(:,:,fishNum));
end
data.ctrl.vib_corrVec = corrVec;

%--- Abl ---
disp('Avg WTs from all ablated fish...')
var = data.abl.vib;
nFish = length(var);
for fishNum = 1:nFish
    disp(['Fish # ' num2str(fishNum)])
    W = var{fishNum}.W;
    W = cat(1,W.head.avg,W.tail.avg);
    if fishNum ==1
        W_all{2} = W;
    else
        W_all{2} = cat(3,W_all{2},W);
    end
end

data.abl.vib_mean = mean(W_all{2},3);
data.abl.vib_std = std(W_all{2},[],3);
data.abl.vib_cv = data.abl.vib_std./(data.abl.vib_mean + 1);
corrVec = nan(1,nFish);
for fishNum = 1:nFish
    corrVec(fishNum) = corr2(data.abl.vib_mean,W_all{2}(:,:,fishNum));
end
data.abl.vib_corrVec = corrVec;
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

%% Plot Figs

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
midRow = size(data.ctrl.vib_mean,1)/2;

%--- Mean head for ctrl --
axes(axH(1))
imagesc(time,freq,abs(data.ctrl.vib_mean(1:midRow,:)).^1)
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
title('CONTROL')
ylabel('Freq (Hz)')

%--- Mean tail for ctrl --
axes(axH(2))
imagesc(time,freq,abs(data.ctrl.vib_mean(midRow+1:end,:)).^1)
set(gca,'clim',cl,'ydir','normal','tickdir','out')
box off
xlabel('Time (ms)')
ylabel('Freq (Hz)')

%--- Mean head for abl --
axes(axH(3))
imagesc(time,freq,abs(data.abl.vib_mean(1:midRow,:)).^1)
set(gca,'clim',cl,'ydir','normal','tickdir','out','ytick',[],'xtick',[],'yaxislocation','right')
box off
ylabel('HEAD')
title('ABLATED')

%--- CV tail --
axes(axH(4))
imagesc(time,freq,abs(data.abl.vib_mean(midRow+1:end,:)).^1)
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
W_diff = data.abl.vib_mean(1:midRow,:) - data.ctrl.vib_mean(1:midRow,:);
% W_diff(abs(W_diff)<5) = 0;
imagesc(time,freq,W_diff);
set(gca,'ydir','normal','tickdir','out','clim',[-50 50])
box off
title('W_{abl} - W_{ctrl}')
ylabel({'HEAD';'Freq (Hz)'})
ch = colorbar;
set(ch,'location','North')

axes(axH(2))
W_diff = data.abl.vib_mean(midRow+1:end,:) - data.ctrl.vib_mean(midRow+1:end,:);
% W_diff(abs(W_diff)<5) = 0;
imagesc(time,freq,W_diff);
set(gca,'ydir','normal','tickdir','out','clim',[-50 50])
box off
xlabel('Time (ms)')
ylabel({'TAIL';'Freq (Hz)'})


% --- Fig 3: Variablity plot for ctrl and abl ---
figure('Name','Reponse variability across fish')
plot(sort(data.ctrl.vib_corrVec),'.-')
hold on
plot(sort(data.abl.vib_corrVec),'ro-')
set(gca,'tickdir','out'), box off
xlabel('Fish #')
ylabel('Response correlation to average response')
lh = legend('Control','Ablated');
set(lh,'location','best')
title('Response variability across fish for control and ablated groups')

