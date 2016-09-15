% 
% %% Plot tail curvatures atop moving fish
% trl = 2;
% pauseDur = 0.1;
% inds = (trl-1)*750 + 1: (trl-1)*750 + 750;
% var1 = IM_proc_crop(:,:,inds);
% var2 = repmat([71 71],numel(inds),1);
% mlInds2 = midlineInds(inds);
% tC = tailCurv(:,:,inds);
% 
% PlayFishTracking2(IM_proc_crop,'fishPos',[71 71], 'midlineInds',midlineInds,'tailAngles',tA,'plotCurv',1,'tailCurv',tailCurv,'frameInds',inds,'pauseDur',pauseDur);

% 
% 
% 
% %% Writing sequences
% outDir = 'c:\users\pujalaa\documents\outDir'
% procDir = fullfile(imgDir,'trl1');
% mkdir(procDir);
% disp('Writing processed images...')
% tic
% trl{1} = Standardize(trl{1});
% for jj = 1:size(trl{1},3)
%     fName = ['Img_' sprintf('%0.4d', jj) '.bmp'];    
%     imwrite(trl{1}(:,:,jj),fullfile(procDir,fName),'bmp')
%     disp(num2str(jj))
% end
% toc
% 
% 
% procDir = fullfile(imgDir,'trl2');
% mkdir(procDir);
% disp('Writing processed images...')
% tic
% trl{2} = Standardize(trl{2});
% for jj = 1:size(trl{2},3)
%     fName = ['Img_' sprintf('%0.4d', jj) '.bmp'];    
%     imwrite(trl{2}(:,:,jj),fullfile(procDir,fName),'bmp')
%     disp(num2str(jj))
% end
% toc

%% Multiband XWT of curvature timeseries

% trl = 1;
% y = curv_trl(:,trl);
% dj = 1/24;
% 
% Wxy = ComputeMultibandXWT([y(:), y(:)], time_trl(:,1),'freqRange',[15 80],'dj',1/48,'timeRange',[0 0.6]);


%% Analyze data
% grps = {'ctrl','abl'};
% parNames = ['amp','per','angVel'];
% par = struct;
% testMat = 0;
% for grpNum = 1:length(grps) % Dim 1   
%     grp = grps{grpNum};
%     blah = data.(grp);    
%     disp(['Group: ' grp])
%     stimTypes = fieldnames(blah);
%     for stimNum = 1:length(stimTypes) % Dim 2      
%         stimType = stimTypes{stimNum};
%         disp(['Stim type: ' stimType])
%         blah2 = blah.(stimType);
%         for fishNum = 1:length(blah2) % Dim 3            
%             blah3 = blah2{fishNum};
%             parNames = fieldnames(blah3);
%             for parNum = 1:length(parNames) % Dim 4                
%                 parName = parNames{parNum};
% %                 eval([parName '.(stimType)= {};']);
%                 blah4 = blah3.(parName);             
%                 for trlNum = 1:length(blah4) % Dim 5                 
%                     if iscell(blah4(trlNum))
%                         blah5 = blah4{trlNum};                       
%                     else
%                         blah5 = blah4(trlNum);
%                     end                    
%                    if ~isempty(blah5)                        
%                        for pkNum = 1:length(blah5) % Dim 6                             
%                            if iscell(blah4)
%                                eval([parName '.(stimType){fishNum} =  TransposeCell(blah4);']);
%                            else
%                                eval([parName '.(stimType){fishNum} =  blah4;'])
%                            end                          
%                        end
%                    end
%                 end                
%             end
%         end        
%     end
% end
% 


%% Analyze data - Dimensions of multidimensional matrix
grps = {'ctrl','abl'};
%# First need to determine dimensions of multidimensional matrix
dim = zeros(1,6);
for grpNum = 1:length(grps) % Dim 1  
    dim(1) = max(dim(1),grpNum);
    grp = grps{grpNum};
    blah = data.(grp);    
    disp(['Group: ' grp])
    stimTypes = fieldnames(blah);
    for stimNum = 1:length(stimTypes) % Dim 2  
        dim(2)= max(dim(2),stimNum);
        stimType = stimTypes{stimNum};
        disp(['Stim type: ' stimType])
        blah2 = blah.(stimType);
        for fishNum = 1:length(blah2) % Dim 3  
            dim(3) = max(dim(3),fishNum);
            blah3 = blah2{fishNum};
            parNames = fieldnames(blah3);
            for parNum = 1:length(parNames) % Dim 4
                dim(4) = max(dim(4),parNum);
                parName = parNames{parNum};
                blah4 = blah3.(parName);             
                for trlNum = 1:length(blah4) % Dim 5
                    dim(5) = max(dim(5),trlNum);
                    if iscell(blah4(trlNum))
                        blah5 = blah4{trlNum};                       
                    else
                        blah5 = blah4(trlNum);
                    end                    
                   if ~isempty(blah5)                        
                       for pkNum = 1:length(blah5) % Dim 6  
                           dim(6) = max(dim(6),pkNum);                                                 
                       end
                   end
                end                
            end
        end        
    end
end
disp(['Dimentions of data matrix = [ ' num2str(dim) ']']);

%% Analyze data - Filling in multidimensional matrix

%# First need to determine dimensions of multidimensional matrix
dataMat = nan(dim);
%##################################
%## ndims(dataMat) = 6 ;
%## Dim1 = # of grps - 'ctrl','abl'
%## Dim2 = # of stim types - 'dark','vib'
%## Dim3 = # of fish
%## Dim4 = # of params - 'bendAmp','bendPer','onset','bendAngVel'
%## Dim5 = # of trls
%## Dim6 = # of pks
%#################################

disp('Getting data matrix...')
for grpNum = 1:length(grps) % Dim 1
    grp = grps{grpNum};
    blah = data.(grp);
    disp(['Group: ' grp])
    stimTypes = fieldnames(blah);
    for stimNum = 1:length(stimTypes) % Dim 2
        stimType = stimTypes{stimNum};
        disp(['Stim type: ' stimType])
        blah2 = blah.(stimType);
        for fishNum = 1:length(blah2) % Dim 3
            blah3 = blah2{fishNum};
            parNames = fieldnames(blah3);
            for parNum = 1:length(parNames) % Dim 4
                parName = parNames{parNum};
                blah4 = blah3.(parName);
                for trlNum = 1:length(blah4) % Dim 5
                    if iscell(blah4(trlNum))
                        blah5 = blah4{trlNum};
                    else
                        blah5 = blah4(trlNum);
                    end
                    if ~isempty(blah5)
                        for pkNum = 1:length(blah5) % Dim 6
                            dataMat(grpNum,stimNum,fishNum,parNum,trlNum,pkNum) = blah5(pkNum);
                        end
                    end
                end
            end
        end
    end
end
disp('Done!')

%%


