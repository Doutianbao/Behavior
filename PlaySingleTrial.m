%% Read info from procData
% clear
procData = OpenMatFile();
tailCurv = procData.tailCurv;
midlineInds = procData.midlineInds;
disp('Reading cropped image stack...')
IM_proc_crop = procData.IM_proc_crop;
if sum(strcmpi(fieldnames(procData),'dsVecs'))
   dsVecs = procData.dsVecs;
end
toc

%%
TrializeTailBendsAndPlot

%% Plot tail curvatures atop moving fish
trl = 14;
plotCurv = 1;
writeVideo = 0;
pauseDur = 0;
nFramesInTrl = procData.nFramesInTrl;
inds = (trl-1)*nFramesInTrl + 1: (trl-1)*nFramesInTrl + nFramesInTrl;
% var1 = IM_proc_crop(:,:,inds);
% var2 = repmat([71 71],numel(inds),1);
% mlInds2 = midlineInds(inds);
% tC = tailCurv(:,:,inds);
fp = (size(IM_proc_crop)-1)/2 + 1;
fp = fp(1:2);
PlayFishTracking(IM_proc_crop,'fishPos',fp, 'midlineInds',midlineInds,'tailAngles',tA,...
    'plotCurv',plotCurv,'tailCurv',tailCurv,'frameInds',inds,'pauseDur',pauseDur,'lineClr','r','writeVideo',writeVideo);
