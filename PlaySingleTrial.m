
%% Plot tail curvatures atop moving fish
trl = 9;
pauseDur = 0.1;
inds = (trl-1)*750 + 1: (trl-1)*750 + 750;
var1 = IM_proc_crop(:,:,inds);
var2 = repmat([71 71],numel(inds),1);
mlInds2 = midlineInds(inds);
tC = tailCurv(:,:,inds);

PlayFishTracking(IM_proc_crop,'fishPos',[71 71], 'midlineInds',midlineInds,'tailAngles',tA,'plotCurv',1,'tailCurv',tailCurv,'frameInds',inds,'pauseDur',pauseDur);
