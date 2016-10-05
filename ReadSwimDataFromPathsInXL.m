
%% Path to xl sheet
clear, close all
path = 'S:\Avinash\Ablations and behavior';
fName = 'Ablation data summary.xlsx';
fPath = fullfile(path,fName);

%% Read data
[data,pData, dataMat, dimLbls] = ReadSwimDataFromPaths(fPath);
grpData = struct;
grpData.data = data;
grpData.paths = pData;
grpData.dataMat = dataMat;
grpData.dimLbls = dimLbls;


%% Save data as .mat file
saveOrNot = input('Save data as .mat file? (y/n) ','s');
if strcmpi(saveOrNot,'y')
    fName2 =['ElicitedSwimData_' data.ablationType '_' datestr(now,30) '.mat'];
    save(fullfile(path,fName2), 'grpData');
else
    disp('Data not saved!')
end


%% Save data as .csv file
tic
saveOrNot = input('Save data as .csv file? (y/n) ','s');
if strcmpi(saveOrNot,'y')
    disp('Creating xls type array')
    data_cell = CreateCellArray(data,dataMat);
    fName2 =['ElicitedSwimData_' data.ablationType '_' datestr(now,30) '.csv'];
    xlswrite(fName2,data_cell)
    disp('Saved!')
else
    disp('Data not saved!')
end
toc




