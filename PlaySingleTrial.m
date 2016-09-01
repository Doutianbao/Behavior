%% Read info from procData
[fileName,pathName] = uigetfile('*.mat');
disp('Reading data...')
tic
procData = matfile(fullfile(pathName,fileName));
tailCurv = procData.tailCurv;
midlineInds = procData.midlineInds;
IM_proc_crop = procData.IM_proc_crop;
toc

%%
TrializeTailBendsAndPlot

%% Plot tail curvatures atop moving fish
trl = 15;
pauseDur = [];
inds = (trl-1)*750 + 1: (trl-1)*750 + 750;
% var1 = IM_proc_crop(:,:,inds);
% var2 = repmat([71 71],numel(inds),1);
% mlInds2 = midlineInds(inds);
% tC = tailCurv(:,:,inds);
fp = (size(IM_proc_crop)-1)/2 + 1;
fp = fp(1:2);
PlayFishTracking(IM_proc_crop,'fishPos',fp, 'midlineInds',midlineInds,'tailAngles',tA,...
    'plotCurv',1,'tailCurv',tailCurv,'frameInds',inds,'pauseDur',pauseDur,'lineClr','r');
