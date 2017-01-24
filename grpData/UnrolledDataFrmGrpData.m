function data_xls = UnrolledDataFrmGrpData(grpData)
%UnrolledDataFrmGrpData - Given grpData that contains all the ablation
%   information for all the ablation groups, returns a data structure whose
%   fields contain all the data unrolled into a format that makes it easy
%   to write into a xl sheet that can later be read into R.
% out = UnrolledGrpData(grpData);
% 
% Inputs:
% grpData - Contains all the ablation data for all the ablation groups
% Outputs:
% out - Unrolled data

grps = fieldnames(grpData);
stimTypes = {'vib', 'dark'};
trtmnts = {'ctrl','abl'};
flds = {'bendAmp','bendPer','onset'};

%% Creating data struct
nFlds = length(flds);
data_xls = struct;
count = 0;
count_entry = 0;
for grp = 1:length(grps)
    for trtmnt = 1:length(trtmnts)
        for stim = 1:length(stimTypes)
            blah = grpData.(grps{grp}).(trtmnts{trtmnt}).(stimTypes{stim}).procData;
            for fish = 1:length(blah)
                disp([upper(grps{grp}) ', ' trtmnts{trtmnt} ', ' stimTypes{stim} ', fish # ' num2str(fish)])                
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



end

