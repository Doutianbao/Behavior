
%  --- Takes grpData and unrolls into format suitable for excel sheet

%% Inputs

saveDir = 'S:\Avinash\Ablations and behavior\GrpData';
grps = fieldnames(grpData);
stimTypes = {'vib', 'dark'};
trtmnts = {'ctrl','abl'};
flds = {'bendAmp','bendPer','onset'};

%% Creating data struct
tic
data_xls = UnrolledDataFrmGrpData(grpData);
toc

% --- Append to grpData for faster retrieval and save
grpData.data_xls = data_xls;
disp('Saving group ablation data...')
ts = datestr(now,30);
fName = ['Ablation data for all groups_' ts];
save(fullfile(saveDir,fName),'grpData','-v7.3')
disp('Done!')

%% Creating cell array and writing to xl sheet

disp('Writing to xl sheet...')
fldNames = fieldnames(data_xls);
data_array = cell(numel(data_xls.bendAmp)+1,length(fieldnames(data_xls)));
for fld = 1:length(fldNames)
    fldName = fldNames{fld};
    fldName(1) = upper(fldName(1));
    data_array{1,fld} = fldName;
    vars = data_xls.(fldNames{fld});
    if iscell(vars(1))
        data_array(2:end,fld) = data_xls.(fldNames{fld});
    else
        data_array(2:end,fld) = num2cell(data_xls.(fldNames{fld}));
    end    
end

fName = ['Peak Info for Group Data_' ts];
xlswrite(fullfile(saveDir,fName),data_array)
disp(['Written at ' saveDir])


