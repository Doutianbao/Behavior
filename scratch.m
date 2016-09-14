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

trl = 1;
y = curv_trl(:,trl);
dj = 1/24;

Wxy = ComputeMultibandXWT([y(:), y(:)], time_trl(:,1),'freqRange',[15 80],'dj',1/48,'timeRange',[0 0.6]);



