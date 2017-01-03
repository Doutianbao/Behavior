
clear all

%% Setting inputs
xlsPath = 'S:\Avinash\Ablations and behavior\Ablation data summary.xlsx';
saveDir = 'S:\Avinash\Ablations and behavior\GrpData';
xLim_tap = [-50 650];
xLim_dark = [-100 300];
onsetAlign = 1;
plotOrNot = 0;
sigmaXY = 100; %(default = NaN);
cLim = [0 300];
freqRange_dark = [5 45];

%% Ctrl, tap
disp('Getting ctrl, tap data...')
[paths,paramVals, filtInds] = GetFilteredPathsFromXLS(xlsPath);
for pp  = 1:length(paths)
    [path_curr, file_curr] = fileparts(paths{pp});
    if ~strcmpi(file_curr,'proc')
        paths{pp} = fullfile(paths{pp},'proc');
    end
end

% out = AnalyzeFreeSwims_nCycles_batch(paths(11:end), 'xLim',[0 750],'paramList',...
%     {'bodyAmp'});

[~, procData] = GetFishWaves_group(paths(), 'saveToProc',1,'xLim',xLim_tap,...
    'onsetAlign',1,'sigmaXY',sigmaXY,'plotOrNot',plotOrNot,'cLim',cLim,...
    'traceType','both','diffOrNot',0);

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
%     dispNext = input('Display next set of figures?(y/n): ','s');
%     if strcmpi(dispNext,'n')
%         break
%     end
end
disp('Completed Ctrl, vib')

%% Abl tap
disp('Getting abl, tap data...')
[paths,~] = GetFilteredPathsFromXLS(xlsPath,filtInds);
for pp  = 1:length(paths)
    [path_curr, file_curr] = fileparts(paths{pp});
    if ~strcmpi(file_curr,'proc')
        paths{pp} = fullfile(paths{pp},'proc');
    end
end

% out = AnalyzeFreeSwims_nCycles_batch(paths(1:end), 'xLim',[0 750],'paramList',...
%     {'bodyAmp'});

[~, procData] = GetFishWaves_group(paths(1:end-1), 'saveToProc',1,'xLim',xLim_tap,...
    'onsetAlign',1,'sigmaXY',sigmaXY,'plotOrNot',plotOrNot,'cLim',cLim,...
    'traceType','both','diffOrNot',0);


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
%     dispNext = input('Display next set of figures?(y/n): ','s');
%     if strcmpi(dispNext,'n')
%         break
%     end
end
disp('Completed abl, vib')

%% Ctrl, dark flash
disp('Getting ctrl, dark flash data...')
[paths,~] = GetFilteredPathsFromXLS(xlsPath,filtInds);
for pp  = 1:length(paths)
    [path_curr, file_curr] = fileparts(paths{pp});
    if ~strcmpi(file_curr,'proc')
        paths{pp} = fullfile(paths{pp},'proc');
    end
end

% out = AnalyzeFreeSwims_nCycles_batch(paths(9:end), 'xLim',[0 1500],'paramList',...
%     {'bodyAmp'});

[~, procData] = GetFishWaves_group(paths, 'saveToProc',1,'xLim',xLim_dark,...
    'onsetAlign',onsetAlign,'sigmaXY',sigmaXY,'plotOrNot',plotOrNot,'cLim',[0 25],...
    'freqRange',freqRange_dark,'traceType','curv','diffOrNot',1);

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
%     dispNext = input('Display next set of figures?(y/n): ','s');
%     if strcmpi(dispNext,'n')
%         break
%     end
end
disp('Completed ctrl, dark flash')

%% Abl, dark flash
disp('Getting abl, dark flash data...')
[paths,~] = GetFilteredPathsFromXLS(xlsPath,filtInds);
for pp  = 1:length(paths)
    [path_curr, file_curr] = fileparts(paths{pp});
    if ~strcmpi(file_curr,'proc')
        paths{pp} = fullfile(paths{pp},'proc');
    end
end

% out = AnalyzeFreeSwims_nCycles_batch(paths, 'xLim',[0 1500],'paramList',...
%     {'bodyAmp'});

[~, procData] = GetFishWaves_group(paths(), 'saveToProc',1,'xLim',xLim_dark,...
    'onsetAlign',onsetAlign,'sigmaXY',sigmaXY,'plotOrNot',plotOrNot,'cLim',[0 25],...
    'freqRange',freqRange_dark,'traceType','curv','diffOrNot',1);

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
%     dispNext = input('Display next set of figures?(y/n): ','s');   
%     if strcmpi(dispNext,'n')
%         break
%     end
end
disp('Completed abl, dark flash')

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
procData = data.abl.vib.procData;
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
delInd = [];
disp('Avg WTs from all control fish...')
procData = data.ctrl.vib.procData;
nFish = length(procData);
fishInds = 1:nFish;
fishInds(delInd) = [];
nFish = length(fishInds);
for fishNum = fishInds(:)'
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
count = 0;
for fishNum = fishInds(:)'
    count = count + 1;
    corrVec(fishNum) = corr2(data.ctrl.vib.mean,W_all{1}(:,:,count));
end
data.ctrl.vib.corrVec = corrVec;

%--- Abl ---
delInd = [];
disp('Avg WTs from all ablated fish...')
procData = data.abl.vib.procData;
nFish = length(procData);
fishInds = 1:nFish;
fishInds(delInd)=[];
nFish = length(fishInds);
for fishNum = fishInds(:)'
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
count = 0;
for fishNum = fishInds(:)'
    count = count + 1;
    corrVec(fishNum) = corr2(data.abl.vib.mean,W_all{2}(:,:,count));
end
data.abl.vib.corrVec = corrVec;
toc

%% The average WT across fish from average WTs for fish across trials - for dark flash stim and curvatures
tic
W_all = cell(2,1);

%--- Ctrl ---
disp('Avg WTs from all control fish...')
procData = data.ctrl.dark.procData;
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
procData = data.abl.dark.procData;
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


%% Vib Stim
% --- Unraveling trls ---
stimType = 'vib';
disp(stimType)
var = {};
var{1} = data.ctrl.(stimType).procData;
var{2} = data.abl.(stimType).procData;
W_all = cell(size(var));
wDims = [];
for trtmnt = 1:2
    disp(['Stim # ' num2str(trtmnt)])
    for fishNum = 1:length(var{trtmnt})
        disp(['Fish # ' num2str(fishNum)])
        W = var{trtmnt}{fishNum}.W;
        if fishNum == 1
            wDims = size(W.curv.coeff{1});
            nTrls = numel(W.trlList);
            W = [W.curv.coeff{:}];
            W = reshape(W,[wDims nTrls]);
            W_all{trtmnt} =  W;
        else
            nTrls = numel(W.trlList);
            W = [W.curv.coeff{:}];
            W = reshape(W,[wDims nTrls]);
            W_all{trtmnt} =  cat(3,W_all{trtmnt},W);
        end
    end
end

% --- Correlations ---
tic
disp('Getting correlations...')
disp('Ctrl')
% D = struct;
nTrls = size(W_all{1},3);
blah = nan(nTrls,2);
trlList = 1:nTrls;
for trl = trlList(:)'
    w1 = abs(W_all{1}(:,:,trl));
    otherTrls = setdiff(trlList,trl);
    w_mu1 = mean(abs(W_all{1}(:,:,otherTrls)),3);
    w_mu2 = mean(abs(W_all{2}),3);
    blah(trl,1) = corr2(w1,w_mu1);
    blah(trl,2) = corr2(w1,w_mu2);
end
data.ctrl.(stimType).corrMat = blah;
toc

disp('Getting correlations...')
disp('Abl')
nTrls = size(W_all{2},3);
blah = nan(nTrls,2);
trlList = 1:nTrls;
for trl = trlList(:)'
    w1 = abs(W_all{2}(:,:,trl));
    otherTrls = setdiff(trlList,trl);
    w_mu1 = mean(abs(W_all{2}(:,:,otherTrls)),3);
    w_mu2 = mean(abs(W_all{1}),3);
    blah(trl,1) = corr2(w1,w_mu1);
    blah(trl,2) = corr2(w1,w_mu2);
end
data.abl.(stimType).corrMat = blah;
toc


%% Dark flash stim
% --- Unraveling trls ---
stimType = 'dark';
var = {};
var{1} = data.ctrl.(stimType).procData;
var{2} = data.abl.(stimType).procData;

W_all = cell(size(var));
wDims = [];
for trtmnt = 1:2
    disp(['Stim # ' num2str(trtmnt)])
    for fishNum = 1:length(var{trtmnt})
        disp(['Fish # ' num2str(fishNum)])
        W = var{trtmnt}{fishNum}.W;
        if fishNum == 1
            wDims = size(W.curv.coeff{1});
            nTrls = numel(W.trlList);
            W = [W.curv.coeff{:}];
            W = reshape(W,[wDims nTrls]);
            W_all{trtmnt} =  W;
        else
            nTrls = numel(W.trlList);
            W = [W.curv.coeff{:}];
            W = reshape(W,[wDims nTrls]);
            W_all{trtmnt} =  cat(3,W_all{trtmnt},W);
        end
    end
end

% --- Correlations ---
tic
disp('Getting correlations...')
disp('Ctrl')
D = struct;
nTrls = size(W_all{1},3);
blah = nan(nTrls,2);
trlList = 1:nTrls;
for trl = trlList(:)'
    w1 = abs(W_all{1}(:,:,trl));
    otherTrls = setdiff(trlList,trl);
    w_mu1 = mean(abs(W_all{1}(:,:,otherTrls)),3);
    w_mu2 = mean(abs(W_all{2}),3);
    blah(trl,1) = corr2(w1,w_mu1);
    blah(trl,2) = corr2(w1,w_mu2);
end
data.ctrl.(stimType).corrMat = blah;
toc

disp('Getting correlations...')
disp('Abl')
nTrls = size(W_all{2},3);
blah = nan(nTrls,2);
trlList = 1:nTrls;
for trl = trlList(:)'
    w1 = abs(W_all{2}(:,:,trl));
    otherTrls = setdiff(trlList,trl);
    w_mu1 = mean(abs(W_all{2}(:,:,otherTrls)),3);
    w_mu2 = mean(abs(W_all{1}),3);
    blah(trl,1) = corr2(w1,w_mu1);
    blah(trl,2) = corr2(w1,w_mu2);
end
data.abl.(stimType).corrMat = blah;
toc


%% Saving data
saveOrNot = input('Save group data?(y/n): ', 's');
if strcmpi(saveOrNot,'y')
    disp('Saving data...')
    tic
    timeStamp = datestr(now,30);
    fName = ['procData_grp_' data.ablationType '_' timeStamp '.mat'];
    save(fullfile(saveDir,fName),'data');
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
fh = figure('Name','Difference WT maps');
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

% ---- Vib Stim -----
% -- Fig 1: Curvature WT for ctrl and abl : arrangement 1 ---
stimType = 'vib';
W1 = data.ctrl.(stimType).mean;
W2 = data.abl.(stimType).mean;
[r_max,c_max] = find(abs(W1) == max(abs(W1(:))));
ax = {};
ax{1} = [0.4 0.8 0 0.2]; % WT, Ctrl
ax{2} = [0.4 0.8 0.4 0.2]; % WT, Abl
ax{3} = [0.2 0.8 0.8 0.2]; % GPS, Ctrl & Abl
ax{4} = [0.4 0.2 0 0]; % Inst mean pow, Ctrl
ax{5} = [0.4 0.2 0.4 0]; % Inst mean pow, Abl

fh = figure('Name','Avg WT for ctrl and abl');
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3},ax{4},ax{5});
W = data.ctrl.(stimType).procData{1}.W;
freq = W.freq;
time = W.time;
cl = W.cLim;

%--- Avg WT, Ctrl ---
axes(axH(1))
imagesc(time,freq,abs(W1))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
title('CONTROL')
ylabel('Freq (Hz)')

%--- Avg WT, Abl ---
axes(axH(2))
imagesc(time,freq,abs(W2))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','ytick',[],'xtick',[])
box off
title('ABLATED')

%--- GPS for Ctrl & Abl --
axes(axH(3))
plot(mean(abs(W1),2),freq,'c')
hold on
plot(mean(abs(W2),2),freq,'m')
set(gca,'ytick',[],'tickdir','out','color','k')
box off
xlim([-inf inf])
ylim([min(freq), max(freq)])
xlabel('Mean Pow')
lh = legend('Ctrl','Abl');
set(lh,'color',[0.7 0.7 0.7],'location','best')

%--- Inst mean pow, Ctrl --
var1 = mean(abs(W1),1);
var2 = mean(abs(W2),1);
yL =[0 max(max(var1),max(var2))];
axes(axH(4))
plot(time,mean(abs(W1),1),'c')
hold on
set(gca,'color','k','tickdir','out')
box off
xlim([time(1) time(end)])
ylim(yL)
xlabel('Time (ms)')

axes(axH(5))
plot(time,mean(abs(W2),1),'m')
hold on
set(gca,'color','k','tickdir','out','ytick',[])
box off
xlim([time(1) time(end)])
ylim(yL)
xlabel('Time (ms)')

% ---- Vib Stim -----
% -- Fig 2: Curvature WT for ctrl and abl : arrangement 2 ---
ax = {};
ax{1} = [0.8 0.4 0 0.6]; % WT, Ctrl
ax{2} = [0.8 0.4 0 0.2]; % WT, Abl
ax{3} = [0.2 0.4 0.8 0.6]; % GPS, Ctrl
ax{4} = [0.2 0.4 0.8 0.2]; % GPS, Abl
ax{5} = [0.8 0.2 0 0]; % Inst mean pow

fh = figure('Name','Avg WT for ctrl and abl');
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3},ax{4},ax{5});

%--- Avg WT, Ctrl ---
axes(axH(1))
imagesc(time,freq,abs(W1))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'CONTROL','Freq (Hz)'})
title('Average Wavelet Transforms')


%--- Avg WT, Abl ---
axes(axH(2))
imagesc(time,freq,abs(W2))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'ABLATED','Freq (Hz)'})


%--- GPS for Ctrl & Abl --
var1 = mean(abs(W1),2);
var2 = mean(abs(data.abl.(stimType).mean),2);
xL = [min(min(var1),min(var2)) max(max(var1),max(var2))];

axes(axH(3))
plot(var1,freq,'c')
box off
xlim(xL)
ylim([min(freq),max(freq)])
set(gca,'color','k','ytick',[],'xtick',[])
title('GPS')

axes(axH(4))
plot(var2,freq,'m')
box off
xlim(xL)
ylim([min(freq),max(freq)])
xlabel('Mean Pow')
set(gca,'color','k','ytick',[])

% -- Inst mean freq --
axes(axH(5))
mf1 = instantaneouswavefreq(W1,freq);
mp1 = instantaneouswavepow(W1);

mf2 = instantaneouswavefreq(data.abl.(stimType).mean,freq);
mp2 = instantaneouswavepow(data.abl.(stimType).mean);
sf = max(max(mp1),max(mp2));
mp1 = mp1/sf;
mp2 = mp2/sf;
mp = max([mp1(:),mp2(:)],[],2);
t = time;
lowInds = find(mp<0.025);
mf1(lowInds) = nan;
mf2(lowInds) = nan;
t(lowInds) = nan;

beforeInds = time<0;
mf1(beforeInds) = nan;
mf2(beforeInds) = nan;
t(beforeInds) = nan;

plot(t,mf1,'c.')
hold on
plot(t,mf2,'m.')
set(gca,'color','k')
xlim([time(1) time(end)])
mf_all  = [mf1(:); mf2(:)];
ylim([min(mf_all) max(mf_all)])
ylabel({'Mean', 'Freq (Hz)'})
lh = legend('Ctrl','Abl');
set(lh,'location','best','color',[0.7 0.7 0.7])

linkaxes(axH([1 2 5]),'x')

%--- Fig 3: Difference WT maps ---
fh = figure('Name','Difference WT maps');
imagesc(time,freq, Standardize(W2-W1)-0.5)
box off
xlabel('Time (ms)')
ylabel('Freq (Hz)')
set(gca,'ydir','normal','tickdir','out','clim',[-0.5 0.5]);
ch =colorbar;
set(ch,'ytick',[-0.3 0.3],'yticklabel',{'Ctrl > Abl','Abl > Ctrl'})
title('$\overline{W}_{abl} - \overline{W}_{ctrl}$','interpreter','latex','fontsize',16)


% --- Fig 3: Variablity plot for ctrl and abl ---
figure('Name','Reponse variability across fish')
plot(sort(data.ctrl.(stimType).corrVec),'.-')
hold on
plot(sort(data.abl.(stimType).corrVec),'ro-')
set(gca,'tickdir','out'), box off
ylim([0.5 1])
xlabel('Fish #')
ylabel('Response correlation to average response')
lh = legend('Control','Ablated');
set(lh,'location','best')
title('Response variability across fish for control and ablated groups')

% ---- Dark flash Stim -----
% -- Fig 1: Curvature WT for ctrl and abl : arrangement 1 ---
stimType = 'dark';
W1 = data.ctrl.(stimType).mean;
W2 = data.abl.(stimType).mean;
[r_max,c_max] = find(abs(W1) == max(abs(W1(:))));
ax = {};
ax{1} = [0.4 0.8 0 0.2]; % WT, Ctrl
ax{2} = [0.4 0.8 0.4 0.2]; % WT, Abl
ax{3} = [0.2 0.8 0.8 0.2]; % GPS, Ctrl & Abl
ax{4} = [0.4 0.2 0 0]; % Inst mean pow, Ctrl
ax{5} = [0.4 0.2 0.4 0]; % Inst mean pow, Abl

fh = figure('Name','Avg WT for ctrl and abl');
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3},ax{4},ax{5});
W = data.ctrl.(stimType).procData{1}.W;
freq = W.freq;
time = W.time;
cl = W.cLim;

%--- Avg WT, Ctrl ---
axes(axH(1))
imagesc(time,freq,abs(W1))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
title('CONTROL')
ylabel('Freq (Hz)')

%--- Avg WT, Abl ---
axes(axH(2))
imagesc(time,freq,abs(W2))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','ytick',[],'xtick',[])
box off
title('ABLATED')

%--- GPS for Ctrl & Abl --
axes(axH(3))
plot(mean(abs(W1),2),freq,'c')
hold on
plot(mean(abs(W2),2),freq,'m')
set(gca,'ytick',[],'tickdir','out','color','k')
box off
xlim([-inf inf])
ylim([min(freq), max(freq)])
xlabel('Mean Pow')
lh = legend('Ctrl','Abl');
set(lh,'color',[0.7 0.7 0.7],'location','best')

%--- Inst mean pow, Ctrl --
var1 = mean(abs(W1),1);
var2 = mean(abs(W2),1);
yL =[0 max(max(var1),max(var2))];
axes(axH(4))
plot(time,mean(abs(W1),1),'c')
hold on
set(gca,'color','k','tickdir','out')
box off
xlim([time(1) time(end)])
ylim(yL)
xlabel('Time (ms)')

axes(axH(5))
plot(time,mean(abs(W2),1),'m')
hold on
set(gca,'color','k','tickdir','out','ytick',[])
box off
xlim([time(1) time(end)])
ylim(yL)
xlabel('Time (ms)')

% ---- Dark flash stim-----
% -- Fig 2: Curvature WT for ctrl and abl : arrangement 2 ---
ax = {};
ax{1} = [0.8 0.4 0 0.6]; % WT, Ctrl
ax{2} = [0.8 0.4 0 0.2]; % WT, Abl
ax{3} = [0.2 0.4 0.8 0.6]; % GPS, Ctrl
ax{4} = [0.2 0.4 0.8 0.2]; % GPS, Abl
ax{5} = [0.8 0.2 0 0]; % Inst mean pow

fh = figure('Name','Avg WT for ctrl and abl');
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3},ax{4},ax{5});

%--- Avg WT, Ctrl ---
axes(axH(1))
imagesc(time,freq,abs(W1))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'CONTROL','Freq (Hz)'})
title('Average Wavelet Transforms')


%--- Avg WT, Abl ---
axes(axH(2))
imagesc(time,freq,abs(W2))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'ABLATED','Freq (Hz)'})


%--- GPS for Ctrl & Abl --
var1 = mean(abs(W1),2);
var2 = mean(abs(W2),2);
xL = [min(min(var1),min(var2)) max(max(var1),max(var2))];


axes(axH(3))
plot(var1,freq,'c')
box off
xlim(xL)
ylim([min(freq),max(freq)])
set(gca,'color','k','ytick',[],'xtick',[])
title('GPS')

axes(axH(4))
plot(var2,freq,'m')
box off
xlim(xL)
ylim([min(freq),max(freq)])
xlabel('Mean Pow')
set(gca,'color','k','ytick',[])


% --- Inst mean freq ---
axes(axH(5))
mf1 = instantaneouswavefreq(W1,freq);
mp1 = instantaneouswavepow(W1);

mf2 = instantaneouswavefreq(W2,freq);
mp2 = instantaneouswavepow(W2);
sf = max(max(mp1),max(mp2));
mp1 = mp1/sf;
mp2 = mp2/sf;
mp = max([mp1(:),mp2(:)],[],2);
t = time;
lowInds = find(mp<0.075);
mf1(lowInds) = nan;
mf2(lowInds) = nan;
t(lowInds) = nan;

beforeInds = time<0;
mf1(beforeInds) = nan;
mf2(beforeInds) = nan;
t(beforeInds) = nan;

plot(t,mf1,'c.')
hold on
plot(t,mf2,'m.')
set(gca,'color','k')
xlim([time(1) time(end)])
mf_all  = [mf1(:); mf2(:)];
ylim([min(mf_all) max(mf_all)])
ylabel({'Mean', 'Freq (Hz)'})
lh = legend('Ctrl','Abl');
set(lh,'location','best','color',[0.7 0.7 0.7])

linkaxes(axH([1 2 5]),'x')

%--- Fig 3: Difference WT maps ---
fh = figure('Name','Difference WT maps');
imagesc(time,freq, Standardize(W2-W1)-0.5)
box off
xlabel('Time (ms)')
ylabel('Freq (Hz)')
set(gca,'ydir','normal','tickdir','out','clim',[-0.5 0.5]);
ch =colorbar;
set(ch,'ytick',[-0.3 0.3],'yticklabel',{'Ctrl > Abl','Abl > Ctrl'})
title('$|\overline{W}_{abl} - \overline{W}_{ctrl}|$','interpreter','latex','fontsize',16)


% --- Fig 5: Variablity plot for ctrl and abl ---
figure('Name','Reponse variability across fish')
plot(sort(data.ctrl.(stimType).corrVec),'.-')
hold on
plot(sort(data.abl.(stimType).corrVec),'ro-')
set(gca,'tickdir','out'), box off
ylim([0.5 1])
xlabel('Fish #')
ylabel('Response correlation to average response')
lh = legend('Control','Ablated');
set(lh,'location','best')
title('Response variability across fish for control and ablated groups')

%--- Fig 6: Difference WT maps ater maxnorm ---
fh = figure('Name','Difference WT maps - maxnorm');
imagesc(time,freq, Standardize(abs(W2))-Standardize(abs(W1)))
box off
xlabel('Time (ms)')
ylabel('Freq (Hz)')
set(gca,'ydir','normal','tickdir','out','clim',[-0.5 0.5]);
ch =colorbar;
set(ch,'ytick',[-0.3 0.3],'yticklabel',{'Ctrl > Abl','Abl > Ctrl'})
title('$|\overline{W}_{abl}| - |\overline{W}_{ctrl}|$','interpreter','latex','fontsize',16)


%% Plot Figs - for body curvatures (compact format)

% ---- VIB STIM -----
% -- Fig 1: Curvature WT for ctrl and abl : arrangement 1 ---
stimType = 'vib';
cl = [0 500];
xl = [-50 300];
W1 = data.ctrl.(stimType).mean;
W2 = data.abl.(stimType).mean;
[r_max,c_max] = find(abs(W1) == max(abs(W1(:))));
ax = {};
ax{1} = [0.33 1 0 0]; % WT, Ctrl
ax{2} = [0.33 1 0.33 0]; % WT, Abl
ax{3} = [0.33 1 0.66 0]; % WT, Abl - Ctrl

fh = figure('Name','Avg WT for ctrl and abl');
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3});
W = data.ctrl.(stimType).procData{1}.W;
freq = W.freq;
time = W.time;
% cl = W.cLim;

%--- Avg WT, Ctrl ---
axes(axH(1))
imagesc(time,freq,abs(W1))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out')
box off
title('CONTROL')
ylabel('Freq (Hz)')
% xlim([time(1) time(end)])
xlim(xl)

%--- Avg WT, Abl ---
axes(axH(2))
imagesc(time,freq,abs(W2))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','ytick',[])
box off
title('ABLATED')
xlabel('Time (ms)')
% xlim([time(1) time(end)])
xlim(xl)
cb = colorbar;
set(cb,'location','east','ycolor','w','ytick',[0 250 500],'box','off')

%--- Diff maps, Abl - Ctrl ---
axes(axH(3))
imagesc(time,freq,W2-W1)
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',[-cl(2)*0.4 cl(2)*0.4],'ydir','normal','tickdir','out','ytick',[])
box off
title('ABL - CTRL')
% xlim([time(1) time(end)])
xlim(xl)
cb = colorbar;
set(cb,'location','east','ytick',[-200 0 200],'box','off')

% -- Fig 2: Curvature WT for ctrl and abl : arrangement 2 ---
ax = {};
ax{1} = [1 0.33 0 0.66]; % WT, Ctrl
ax{2} = [1 0.33 0 0.33]; % WT, Abl
ax{3} = [1 0.33 0 0]; % WT, Abl - Ctrl

fh = figure('Name','Avg WT for ctrl and abl');
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3});
W = data.ctrl.(stimType).procData{1}.W;
freq = W.freq;
time = W.time;
cl = W.cLim;
cl = [0 500];

%--- Avg WT, Ctrl ---
axes(axH(1))
imagesc(time,freq,abs(W1))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'CONTROL';'Freq (Hz)'})
% xlim([time(1) time(end)])
xlim(xl)
title(['Average wavelet spectra of ' stimType '-elicited responses for control and ablated fish'])

%--- Avg WT, Abl ---
axes(axH(2))
imagesc(time,freq,abs(W2))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'ABLATED','Freq (Hz)'})
xlabel('Time (ms)')
% xlim([time(1) time(end)])
xlim(xl)
cb = colorbar;
set(cb,'location','east','ycolor','w','ytick',[0 250 500],'box','off')

%--- Diff map, Abl - Ctrl --
axes(axH(3))
imagesc(time,freq,W2-W1)
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',[-cl(2)*0.4 cl(2)*0.4],'ydir','normal','tickdir','out')
box off
ylabel({'ABL - CTRL'; 'Freq (Hz)'})
xlabel('Time (ms)')
% xlim([time(1) time(end)])
xlim(xl)
cb = colorbar;
set(cb,'location','east','ytick',[-200 0 200],'box','off')


% ---- DARK FLASH STIM -----
% -- Fig 1: Curvature WT for ctrl and abl : arrangement 1 ---
stimType = 'dark';
mult = 15;
% cl = [0 50];
W1 = data.ctrl.(stimType).mean * mult;
W2 = data.abl.(stimType).mean * mult;
[r_max,c_max] = find(abs(W1) == max(abs(W1(:))));
ax = {};
ax{1} = [0.33 1 0 0]; % WT, Ctrl
ax{2} = [0.33 1 0.33 0]; % WT, Abl
ax{3} = [0.33 1 0.66 0]; % WT, Abl - Ctrl

fh = figure('Name','Avg WT for ctrl and abl');
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3});
W = data.ctrl.(stimType).procData{1}.W;
freq = W.freq;
time = W.time;
% cl = W.cLim;


%--- Avg WT, Ctrl ---
axes(axH(1))
imagesc(time,freq,abs(W1))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out')
box off
title('CONTROL')
ylabel('Freq (Hz)')
% xlim([time(1) time(end)])
xlim(xl)

%--- Avg WT, Abl ---
axes(axH(2))
imagesc(time,freq,abs(W2))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','ytick',[])
box off
title('ABLATED')
xlabel('Time (ms)')
% xlim([time(1) time(end)])
xlim(xl)
cb = colorbar;
set(cb,'location','east','ycolor','w','ytick',[0 250 500],'box','off')

%--- Diff map, Abl - Ctrl --
axes(axH(3))
imagesc(time,freq,W2-W1)
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',[-cl(2)*0.4 cl(2)*0.4],'ydir','normal','tickdir','out','ytick',[])
box off
title('ABL - CTRL')
% xlim([time(1) time(end)])
xlim(xl)
cb = colorbar;
set(cb,'location','east','ytick',[-200 0 200],'box','off')

% -- Fig 2: Curvature WT for ctrl and abl : arrangement 2 ---
ax = {};
ax{1} = [1 0.33 0 0.66]; % WT, Ctrl
ax{2} = [1 0.33 0 0.33]; % WT, Abl
ax{3} = [1 0.33 0 0]; % WT, Abl - Ctrl

fh = figure('Name','Avg WT for ctrl and abl');
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3});
W = data.ctrl.(stimType).procData{1}.W;

%--- Avg WT, Ctrl ---
axes(axH(1))
imagesc(time,freq,abs(W1))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'CONTROL';'Freq (Hz)'})
% xlim([time(1) time(end)])
xlim(xl)
title(['Average wavelet spectra of ' stimType '-elicited responses for control and ablated fish'])

%--- Avg WT, Abl ---
axes(axH(2))
imagesc(time,freq,abs(W2))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'ABLATED','Freq (Hz)'})
xlabel('Time (ms)')
% xlim([time(1) time(end)])
xlim(xl)
cb = colorbar;
set(cb,'location','east','ycolor','w','ytick',[0 250 500],'box','off')

%--- Diff map, Abl - Ctrl --
axes(axH(3))
imagesc(time,freq,W2-W1)
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',[-cl(2)*0.4 cl(2)*0.4],'ydir','normal','tickdir','out')
box off
ylabel({'ABL - CTRL'; 'Freq (Hz)'})
xlabel('Time (ms)')
% xlim([time(1) time(end)])
xlim(xl)
cb = colorbar;
set(cb,'location','east','ytick',[-200 0 200],'box','off')

