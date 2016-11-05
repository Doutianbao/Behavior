

outDir = 'S:\Avinash\Ablations and behavior\GrpData';
xLim_tap = [-50 650];
xLim_dark = [-100 300];
onsetAlign = 1;
plotOrNot = 0;
data = struct;

%% Ctrl, tap
disp('Getting ctrl, tap data...')
[paths,paramVals] = GetFilteredPathsFromXLS();
for pp  = 1:length(paths)
    [path_curr, file_curr] = fileparts(paths{pp});
    if ~strcmpi(file_curr,'proc')
        paths{pp} = fullfile(paths{pp},'proc');
    end
end
[~, procData] = GetFishWaves_group(paths, 'saveToProc',1,'xLim',xLim_tap,'onsetAlign',onsetAlign,'sigmaXY',nan,'plotOrNot',0);
data.ctrl.vib = procData;

for fish = 1:length(procData);
    [outDir, ~] = fileparts(procData{fish}.Properties.Source);
    subDirName = 'WT';
    outDir = fullfile(outDir,subDirName);
    if exist(outDir)~=7
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
[paths,paramVals] = GetFilteredPathsFromXLS();
[~, procData] = GetFishWaves_group(paths, 'saveToProc',1,'xLim',xLim,'onsetAlign',onsetAlign,'sigmaXY',nan,'plotOrNot',0);
data.abl.vib = procData;


%% Ctrl, dark flash
xLim = xLim_dark;
disp('Getting ctrl, dark flash data...')
% [paths,paramVals] = GetFilteredPathsFromXLS();
% paths(GetNanIndsFromCellArray(paths))=[];
[~, procData] = GetFishWaves_group(paths(1:2), 'saveToProc',1,'xLim',xLim,'onsetAlign',1,'sigmaXY',nan,'plotOrNot',0);
data.ctrl.dark = procData;


%% Abl, dark flash
disp('Getting abl, dark flash data...')
[paths,paramVals] = GetFilteredPathsFromXLS();
paths(GetNanIndsFromCellArray(paths))=[];
[~, procData] = GetFishWaves_group(paths, 'saveToProc',1,'xLim',xLim,'onsetAlign',1,'sigmaXY',nan,'plotOrNot',0);
data.abl.dark = procData;

data.ablationType = paramVals.AblationType;

%% Saving data
saveOrNot = input('Save group data?(y/n): ', 's');
if strcmpi(saveOrNot,'y')
    disp('Saving data...')
    tic
    timeStamp = datestr(now,30);
    fName = ['procData_grp_BilInter_' timeStamp '.mat'];
    save(fullfile(outDir,fName),'data');
    toc
end



%% Plotting a few SLC and a few LLC responses in ctrl and M-hom abl fish 
procData_grp = struct;
[paths_abl, paramVals] = GetFilteredPathsFromXLS();

ablInds = 3;

[paths_ctrl,paramVals] = GetFilteredPathsFromXLS();



