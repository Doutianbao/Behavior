function procData = AppendPksToProc(procData)
%AppendPksToProc - Get pks (sorted by first or second, left or right, tap
%   or dark) and append to procData

fishPos = procData.fishPos;
midlineInds = procData.midlineInds;
orientation = procData.orientation;
ref = procData.orientation;

%% Inputs for fast
fps = 500; % (30 for slow, 500 for fast)
nFramesInTrl = 750; %(1800 for slow, 750 for fast)
preStimPer = 0.1; % In seconds (1.5 for slow)
mult = 1; % Determines whether preStimPer gets added or subtracted while plotting rasters (-1 for slow)
nTrls = size(fishPos,1)/nFramesInTrl;
tapTrls = 1:2:nTrls;
flashTrls = 2:2:nTrls;
yOff = [220 200];
seg = 3; % [1 - head, 2 = tail, 3 - combined]
motionThr = 4;
peakDetThr = 0.2;

%## Sort all frame indices based on whether they are tap or dark trials
tapInds = [];
for trl = tapTrls
    tapInds = [tapInds, (trl-1)*nFramesInTrl + 1 : trl*nFramesInTrl];
end
darkInds = setdiff(1:size(fishPos,1),tapInds);

%% Getting motion info
imgDims = [size(ref), size(fishPos,1)];
motionInfo = GetMotionInfo(fishPos,orientation,imgDims(1));
time = (0:size(fishPos,1)-1)/fps;
dS_all = motionInfo.dS_all;
trlFrames = (0:nTrls-1)*nFramesInTrl + 1;
dS_all(trlFrames)=0; % Zeroing artifactual distances resulting from trial transitions

%% Correct fish pos, midlineInds, orientation
dS_thr = 50;
S = @(x)(sum(diff(x,[],1).^2,2)).^0.5;
badInds = find(dS_all>dS_thr);
badInds(badInds==1)=[];
badInds(badInds >= (length(dS_all)-1))=[];
iter = 0;
while numel(badInds)>0
    iter = iter + 1;
    disp(['Iter # ' num2str(iter)])
    if iter > 10
        break
    end
    fishPos(badInds,:) = 0.5*(fishPos(badInds-1,:) + fishPos(badInds+1,:));
    dS_all = [0; S(fishPos)];
    dS_all(trlFrames)=0;
    badInds = find(dS_all>dS_thr);
    badInds(badInds==1)=[];
    badInds(badInds >= (length(dS_all)-1))=[];    
end

for indNum = 1:length(badInds)
    ind = badInds(indNum);
    for seg = 1:length(midlineInds{ind})
        midlineInds{ind}{seg} = 0.5*(midlineInds{ind-1}{seg} + midlineInds{ind+1}{seg});
    end
end
orientation_old = orientation;
orientation = GetFishOrientationFromMidlineInds(midlineInds,imgDims(1:2),'s');
orientation = orientation';
motionInfo = GetMotionInfo(fishPos,orientation,imgDims(1));
curv = motionInfo.curv;

%% Getting peaks
minSwimInt = 150e-3;
maxBurstInt = 100e-3;
minBurstInt = 5e-3;
pkThr_1st = 1.5;
pkThr_2nd = 1;
minSwimPts = round(minSwimInt*fps);
maxBurstPts = round(maxBurstInt*fps);
minBurstPts = round(minBurstInt*fps);

ker = gausswin(8); ker = ker/sum(ker);
curv_sm = conv2(curv(:,3),ker(:),'same');
pks.all = GetPks(zscore(curv_sm),'polarity',0,'peakThr',pkThr_1st,'thrType','rel');


dPks = diff(pks.all);
inds = find(dPks >= minSwimPts) + 1;
pks.first =  pks.all(find(dPks>=minSwimPts)+1);
pks.first = [pks.all(1); pks.first(:)];
pks.all = GetPks(zscore(curv_sm),'polarity',0,'peakThr',pkThr_2nd,'thrType','rel');

pks.second = nan(size(pks.first));
for pk = 1:length(pks.first)
    if pk ==16
        a = 1;
    end
    dPks = pks.all-pks.first(pk);
    dPks(dPks<=0)=inf;
    dPks(dPks >=maxBurstPts)=inf;
    dPks(dPks <= minBurstPts)=inf;
    if ~ all(isinf(dPks))
        [~, ind] = min(dPks);
        pks.second(pk) = pks.all(ind);
    end
end
pks.second(isnan(pks.second))=[];

pks.first_tap = intersect(pks.first,tapInds);
pks.first_tap_left = pks.first_tap(curv_sm(pks.first_tap)>0);
pks.first_tap_right = pks.first_tap(curv_sm(pks.first_tap)<0);
pks.first_dark = intersect(pks.first,darkInds);
pks.first_dark_left = pks.first_dark(curv_sm(pks.first_dark)>0);
pks.first_dark_right = pks.first_dark(curv_sm(pks.first_dark)<0);
pks.second_tap = intersect(pks.second,tapInds);
pks.second_tap_left = pks.second_tap(curv_sm(pks.second_tap)>0);
pks.second_tap_right = pks.second_tap(curv_sm(pks.second_tap)<0);
pks.second_dark = intersect(pks.second,darkInds);
pks.second_dark_left = pks.second_dark(curv_sm(pks.second_dark)>0);
pks.second_dark_right = pks.second_dark(curv_sm(pks.second_dark)<0);
pks.curv = curv_sm;
pks.mu_lbl = {'1st_tap_left', '1st_tap_right','1st_dark_left','1st_dark_right',...
    '2nd_tap_left', '2nd_tap_right','2nd_dark_left','2nd_dark_right'};
pks.mu =[mean(curv_sm(pks.first_tap_left)), mean(abs(curv_sm(pks.first_tap_right))),...
    mean(curv_sm(pks.first_dark_left)), mean(abs(curv_sm(pks.first_dark_right))),...
    mean(curv_sm(pks.second_tap_left)), mean(abs(curv_sm(pks.second_tap_right))),...
    mean(curv_sm(pks.second_dark_left)), mean(abs(curv_sm(pks.second_dark_right)))];

procData.Properties.Writable = true;
procData.pks = pks;

%% Plots

%## Pks on full trace
figure('Name','All pks')
plot(time,curv_sm)
hold on
set(gca,'color','k')
plot(time(pks.all),curv_sm(pks.all),'yo')
plot(time(pks.first),curv_sm(pks.first),'g+')
plot(time(pks.second),curv_sm(pks.second),'m+')


%## First pks for tap and dark
figure('Name','First peaks')
plot(1,abs(pks.curv(pks.first_tap_left)),'b.')
hold on
plot(2,abs(pks.curv(pks.first_tap_right)),'r.')
plot(3,abs(pks.curv(pks.first_dark_left)),'g.')
plot(4,abs(pks.curv(pks.first_dark_right)),'m.')
plot(1:4,pks.mu(1:4),'yo-')
set(gca,'color','k')
xlim([0 5])
box off
set(gca,'xtick',1:4,'xticklabel',pks.mu_lbl(1:4))
title('1st peaks for tap and dark trls sorted by left and right turn')

%## Second pks for tap and dark
figure('Name','Second peaks')
plot(1,abs(pks.curv(pks.second_tap_left)),'b.')
hold on
plot(2,abs(pks.curv(pks.second_tap_right)),'r.')
plot(3,abs(pks.curv(pks.second_dark_left)),'g.')
plot(4,abs(pks.curv(pks.second_dark_right)),'m.')
plot(1:4,pks.mu(5:8),'yo-')
set(gca,'color','k')
xlim([0 5])
box off
set(gca,'xtick',1:4,'xticklabel',pks.mu_lbl(5:8))
title('2nd peaks for tap and dark trls sorted by left and right turn')


%## 1st and 2nd pks on one plot
figure('Name','First & Second peaks')
% plot(1,abs(pks.curv(pks.first_tap_left)),'b.')
% hold on
% plot(2,abs(pks.curv(pks.first_tap_right)),'r.')
% plot(3,abs(pks.curv(pks.first_dark_left)),'g.')
% plot(4,abs(pks.curv(pks.first_dark_right)),'m.')
plot(1:4,pks.mu(1:4),'yo-')
set(gca,'color','k')
xlim([0 5])
box off
% set(gca,'xtick',1:4,'xticklabel',pks.mu_lbl(1:4))
% title('Second peaks for tap and dark trls sorted by left and right turn')
% 
% plot(1.5,abs(pks.curv(pks.second_tap_left)),'b+')
hold on
% plot(2,abs(pks.curv(pks.second_tap_right)),'r+')
% plot(3.5,abs(pks.curv(pks.second_dark_left)),'g+')
% plot(4.5,abs(pks.curv(pks.second_dark_right)),'m+')
plot(1:4,pks.mu(5:8),'g+-')
set(gca,'color','k')
xlim([0 5])
lh = legend('1st pks','2nd pks');
set(lh,'color','w')
box off
% set(gca,'xtick',1:4,'xticklabel',pks.mu_lbl(5:8))
% title('Second peaks for tap and dark trls sorted by left and right turn')



end

