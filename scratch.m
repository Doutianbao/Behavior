

paramNames = fieldnames(data.abl.vib{1})';
for paramNum = 1:length(paramNames)
    paramNames{paramNum}(1) = upper(paramNames{paramNum}(1));
end
data_cell ={'AblationType', 'AblationBool','StimType','FishNum','SessionNum','TrlNum','BendNum'};
data_cell = [data_cell, paramNames];
ablationType = data.ablationType;
ablationType(strfind(ablationType, ' '))=[];
ablatedOrNotGrps = {'abl','ctrl'};
totalFishCount = 0;
sessionCount = 0;
abGrpNum = 0;
tic
for abGrp = ablatedOrNotGrps
    abGrpNum  = abGrpNum + 1;
    disp(abGrp{1});
    if strcmpi(abGrp{1},'abl')
        ablatedBool = true;
    elseif strcmpi(abGrp{1},'ctrl')
        ablatedBool = false;
    end
    var = data.(abGrp{1});
    stimTypes = fieldnames(var);
    for ss = 1:length(stimTypes)
        stimType = stimTypes{ss};
        disp(stimType)
        var = data.(abGrp{1}).(stimType);
        sessionCount = sessionCount + 1;
        for fishNum  = 1:length(var)
            totalFishCount = totalFishCount + 1;
            var = data.(abGrp{1}).(stimType){fishNum};
            params = fieldnames(var);
            for paramNum = 1:length(params)
                param = params{paramNum};
                disp(param)
                var = data.(abGrp{1}).(stimType){fishNum}.(param);
                for trlNum = 1:length(var)
                    if iscell(var)
                        var_bend = var{trlNum};
                    else
                        var_bend = var;
                    end
                    if size(data_cell,1)==1
                        valArray = cell(1,length(params));
                    end
                    valArray{paramNum} = [valArray{paramNum};  var_bend(:)];
                    if ~isempty(var_bend)
                        for bendNum = 1:length(var_bend)
                            paramVals = squeeze(dataMat(abGrpNum,ss,fishNum,:,trlNum,bendNum));
                            paramVals = paramVals(:)';
                            valArray = cell(size(paramVals));
                            for ii = 1:length(valArray)
                                valArray{ii} = paramVals(ii);
                            end
                            data_line = [{ablationType,ablatedBool,stimType,totalFishCount,sessionCount,trlNum,bendNum}, valArray];
                            data_cell = [data_cell; data_line];
                        end
                    else
                        paramVals = nan(1,length(params));
                        valArray = cell(size(paramVals));
                        for ii = 1:length(valArray)
                            valArray{ii} = paramVals(ii);
                        end
                        bendNum = nan;
                        data_line = [{ablationType,ablatedBool,stimType,totalFishCount,sessionCount,trlNum,bendNum}, valArray];
                        data_cell = [data_cell; data_line];
                    end                    
                end
            end
        end
    end    
end
toc

