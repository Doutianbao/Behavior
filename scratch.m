


outDir = 'S:\Avinash\Ablations and behavior\GrpData';

disp('Processing Bilalteral Intermediate group...')
tic
stimTypes = {'vib','dark'};
data = struct;
for ablOrNot = 1:2
    if ablOrNot ==1
        for stimType = 1:length(stimTypes)
            if strcmpi(stimTypes{stimType},'vib')
                disp('Processing controls...')
                disp('Processing vib stimulus...')
                xLim  = [-50 750];
                paths = GetFilteredPathsFromXLS();
                [~, procData] = GetFishWaves_group(paths, 'saveToProc',1,'xLim',xLim,'onsetAlign',1,'sigmaXY',nan,'plotOrNot',0);
                data.ctrl.vib = procData;
            else
                disp('Processing controls...')
                disp('Processing dark flash stimulus...')
                xLim =[-99 300];
                paths = GetFilteredPathsFromXLS();
                [~, procData] = GetFishWaves_group(paths, 'saveToProc',1,'xLim',xLim,'onsetAlign',1,'sigmaXY',nan,'plotOrNot',0);
                data.ctrl.dark = procData;
            end
        end
    else
        for stimType = 1:length(stimTypes)
            if strcmpi(stimTypes{stimType},'vib')
                disp('Processing ablated...')
                disp('Processing vib stimulus...')
                xLim  = [-50 750];
                paths = GetFilteredPathsFromXLS();
                pause(5)
                [~, procData] = GetFishWaves_group(paths, 'saveToProc',1,'xLim',xLim,'onsetAlign',1,'sigmaXY',nan,'plotOrNot',0);
                data.abl.vib = procData;
            else
                disp('Processing ablated...')
                disp('Processing dark flash stimulus...')
                xLim =[-99 300];
                paths = GetFilteredPathsFromXLS();
                pause(5)
                [~, procData] = GetFishWaves_group(paths, 'saveToProc',1,'xLim',xLim,'onsetAlign',1,'sigmaXY',nan,'plotOrNot',0);
                data.abl.dark = procData;
            end
        end
    end
end
toc
break;

%% Saving data
timeStamp = datestr(now,30);
fName = ['procData_grp_BilInter_' timeStamp '.mat'];
save(fullfile(outDir,fName),'data');
toc