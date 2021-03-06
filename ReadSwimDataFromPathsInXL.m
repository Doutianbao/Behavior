
%% Path to xl sheet
clear, close all
fPath = 'S:\Avinash\Ablations and behavior';
fName = 'Ablation data summary.xlsx';
fullPath = fullfile(fPath,fName);

%% Read data
[data,pData, dataMat, dimLbls] = ReadSwimDataFromPaths(fullPath);
grpData = struct;
grpData.data = data;
grpData.paths = pData;
grpData.dataMat = dataMat;
grpData.dimLbls = dimLbls;



%% Save data as .mat file
saveOrNot = input('Save data as .mat file? (y/n) ','s');
if strcmpi(saveOrNot,'y')
    fName2 =['ElicitedSwimData_' data.ablationType '_' datestr(now,30) '.mat'];
    save(fullfile(fPath,fName2), 'grpData');
else
    disp('Data not saved!')
end



%% Save data as .csv file
tic
saveOrNot = input('Save data as .csv file? (y/n) ','s');
if strcmpi(saveOrNot,'y')
    disp('Creating xls type array')
    data_cell = CreateCellArray_fast(data,dataMat);
    fName2 =['ElicitedSwimData_' data.ablationType '_' datestr(now,30) '.xlsx'];
    xlswrite(fullfile(fPath,fName2),data_cell)
    disp('Saved!')
else
    disp('Data not saved!')
end
toc




