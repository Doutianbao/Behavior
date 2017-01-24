
%% Vib stim tracking example
% Data source = S:\Avinash\Ablations and behavior\Ventral RS\20160929\20160929_behavior\Fish1_ctrl1\fastDir_09-29-16-191808\vib\proc
% First load procData into matlab workspace as
%   procData = load(fullfile(procDir, procFileName));

saveDir = 'S:\Avinash\Ablations and behavior\GrpData\Session 20170121\FishTracking Videos'

frameInds = 1:350;
pauseDur = 0.01;


PlayFishTrials(procData,13,'bendAngleDispMode','total','nDispPts',125,'frameInds',frameInds,...
    'writeVideo',1,'saveDir',saveDir,'pauseDur',pauseDur);


%% Dark flash stim tracking example
% Data source = S:\Avinash\Ablations and behavior\Ventral RS\20160929\20160929_behavior\Fish1_ctrl1\fastDir_09-29-16-191808\dark\proc\procData_20161115T203028.mat


PlayFishTrials(procData,14,'bendAngleDispMode','total','nDispPts',125,'frameInds',frameInds,...
    'writeVideo',1,'saveDir',saveDir,'pauseDur',pauseDur);
