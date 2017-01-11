
%  --- Takes grpData and unrolls into format suitable for excel sheet

%% Inputs

saveDir = 'S:\Avinash\Ablations and behavior\GrpData';
grps = fieldnames(grpData);
stimTypes = {'vib', 'dark'};
trtmnts = {'ctrl','abl'};
flds = {'bendAmp','bendPer','onset'};

%% Creating data struct

nFlds = length(flds);
data_xls = struct;
count = 0;
count_pt = 0;
count_entry = 0;
tic
for grp = 1:length(grps)
    for trtmnt = 1:length(trtmnts)
        for stim = 1:length(stimTypes)
            blah = grpData.(grps{grp}).(trtmnts{trtmnt}).(stimTypes{stim}).procData;
            for fish = 1:length(blah)
                disp([upper(grps{grp}) ', ' trtmnts{trtmnt} ', ' stimTypes{stim} ', fish # ' num2str(fish)])
                toc
                foo_fish = blah{fish}.elicitedSwimInfo;
                if ~isempty(foo_fish)
                    nTrls = length(foo_fish.(flds{1}));
                    for fld = 1:nFlds
                        count = count_entry;
                        fldName = flds{fld};
                        foo_field = foo_fish.(flds{fld});
                        nBends = 0;
                        for trl = 1:nTrls
                            if length(foo_field) >= trl
                                nBends = max(nBends,length(foo_field{trl}));
                            end
                        end
                        for trl = 1:nTrls
                            if length(foo_field) >=trl
                                foo_trl = foo_field{trl};
                                for bend = 1:nBends
                                    count = count + 1;
                                    data_xls.GrpName{count} = grps{grp};
                                    data_xls.Treatment{count} = trtmnts{trtmnt};
                                    data_xls.StimType{count} = stimTypes{stim};
                                    data_xls.FishNum(count) = fish;
                                    data_xls.TrlNum(count) = trl;
                                    data_xls.BendNum(count) = bend;
                                    if bend > length(foo_trl)
                                        data_xls.(flds{fld})(count) = NaN;
                                    else
                                        data_xls.(flds{fld})(count) = foo_trl(bend);
                                    end
                                    
                                end
                            end
                        end
                    end
                    count_entry = count;
                else
                    error('No elicited swim info')
                end
            end
        end
    end
end
toc

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

ts = datestr(now,30);
fName = ['Peak Info for Group Data_' ts];
xlswrite(fullfile(saveDir,fName),data_array)
disp(['Written at ' saveDir])


