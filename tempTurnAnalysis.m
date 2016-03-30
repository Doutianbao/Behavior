%
% %% Swim analysis
% x = tracexy(:,1);
% y = tracexy(:,2);
% dx = diff(x);
% dy = diff(y);
% dS = sqrt(dx.^2 + dy.^2);
% noMovInds  = find(dS < 3);
% movInds = setdiff(1:length(dx),noMovInds);
% startInds = [];
% initInd = movInds(1);
% for jj = 2:length(movInds)
%     if movInds(jj)-movInds(jj-1) > 1
%         initInd = movInds(jj-1);
%         nextInd = movInds(jj);
%         startInds = [startInds; initInd;nextInd];
%     end
% end
% startInds = unique(startInds);
% x_mov = x(startInds);
% y_mov = y(startInds);
% dS_mov = sqrt(diff(x_mov).^2 + diff(y_mov).^2);
%
% orientation = orientation;
% orient_mov = orientation(startInds);
% [oX,oY] = pol2cart(orient_mov*pi/180,1);
% clear i
% orVec = oX + oY*i;
% dOr = angle(orVec(1:end-1,:).*conj(orVec(2:end,:)))*180/pi;
% straightInds = find(dOr > -5 & dOr < 5);
% leftInds = find(dOr >5);
% rightInds = find(dOr <-5);
% dOr_str = dOr(straightInds);
% dOr_left = dOr(leftInds);
% dOr_right = -dOr(rightInds);
% [orCount_left,orVals] = hist(dOr_left,[10:10:360]);
% [orCount_right,~] = hist(dOr_right,[10:10:360]);
% orMat = [orCount_left(:),orCount_right(:)];
% figure('Name','Turn histogram')
% bar(orVals,orMat,1)
% xlim([0 130])
% box off
% legend('Left','Right')
% set(gca,'tickdir','out','xtick',orVals)
% xlabel('Turn angle')
% ylabel('Count')
% title('Turn angle histogram')
%
% %% Swim distance by turn
% dS_left = dS_mov(leftInds);
% dS_right = dS_mov(rightInds);
% dS_straight = dS_mov(straightInds);
%
% [dslCount,dslVals] = hist(dS_left,20);
% [dsrCount,dsrVals] = hist(dS_right,20);
% [dssCount,dssVals] = hist(dS_straight,20);
% ker = gausswin(5);
% ker = ker/sum(ker);
% dslCount = conv2(dslCount(:),ker(:),'same');
% dsrCount = conv2(dsrCount(:),ker(:),'same');
% dssCount = conv2(dssCount(:),ker(:),'same');
% figure('Name','Swim dist histogram')
% plot(dslVals,dslCount,'linewidth',2)
% hold on
% plot(dsrVals,dsrCount,'r','linewidth',2)
% plot(dssVals,dssCount,'k','linewidth',2)
% xlim([-inf inf])
% ylim([-inf inf])
% box off
% set(gca,'tickdir','out')
% xlabel('Swim distance')
% ylabel('Count')
% legend('Left','Right','Straight')
% title('Swim distance by turn distribution')
%
% %% Velocity stuff
% timeVec = (0:length(tracexy)-1)*(1/30); % Since frame rate = 30fps
% timeVec_mov = timeVec(startInds);
% dT_mov = diff(timeVec_mov);
% vel_mov = dS_mov(:)./dT_mov(:);
% vel_left = vel_mov(leftInds);
% vel_right = vel_mov(rightInds);
% vel_straight = vel_mov(straightInds);
% [vlCount,vlVals] = hist(vel_left,50);
% [vrCount,vrVals] = hist(vel_right,50);
% [vsCount,vsVals] = hist(vel_straight,50);
% vlCount = conv2(vlCount(:),ker(:),'same');
% vrCount = conv2(vrCount(:),ker(:),'same');
% vsCount = conv2(vsCount(:),ker(:),'same');
% figure('Name','Velocity by turn distribution')
% plot(vlVals,vlCount,'linewidth',2)
% hold on
% plot(vrVals,vrCount,'r','linewidth',2)
% % plot(vsVals,vsCount,'k','linewidth',2)
% box off
% xlim([50 700])
% ylim([-inf inf])
% xlabel('Swim Vel (plx/sec)')
% ylabel('Count')
% title('Swim vel by turn distribution')
% set(gca,'tickdir','out')
% legend('Left','Right')
% % legend('Left','Right','Straight')
%
%
%
% break

%% Plotting all position- and orientation-adjusted spont swim trajectories with velocity modulation of hue or alpha

%## Get all motion info
% motionInfo = GetMotionInfo(tracexy,orientation,size(IM_proc,1),5);
%  motionInfo = GetMotionInfo(tracexy,orientation,imgDims(1),5)
%## Plot all adjusted trajectories
figure('Name','Position & orientation adjusted spont swim trajectories')
set(gca,'tickdir','out','color','k')
x_max= 0;
y_max = 0;
y_min = 0;
for ep = 2:length(motionInfo.traj_adj)
    %     plot(motionInfo.traj_adj{ep}(:,1),motionInfo.traj_adj{ep}(:,2),'-')
    x = motionInfo.traj_adj{ep}(:,1);
    y = motionInfo.traj_adj{ep}(:,2);
    x_max = max(x_max,max(abs(x)));
    y_max = max(y_max,max(y));
    y_min = min(y_min,min(y));
    patchline(x,y,'edgealpha',0.2,'linewidth',1.1,'edgecolor','g')    
    hold on
end
% plot(0,0,'ro','markersize',10,'linewidth',2)
plot([-x_max x_max],[0 0],'r--')
plot([0 0],[y_min,y_max],'r--')
box off
axis image
xlim([-x_max x_max])
ylim([y_min,y_max])
title('Position & orientation adjusted spont swim trajectories')

 PlotAngularHist(motionInfo.traj_angle*pi/180,50)
 
 figure
 polar(motionInfo.traj_angle_lim*pi/180,motionInfo.traj_speed,'.')
 title('Traj angle vs traj speed (limited)')
 
 figure
 polar(motionInfo.traj_angle_lim*pi/180,motionInfo.traj_vel,'.')
 title('Traj angle vs traj vel (limited)')
  
 figure
 polar(motionInfo.traj_angle_lim*pi/180,motionInfo.traj_angVel,'.')
  title('Traj angle vs traj angular vel (limited)')
  
 %% Binned trajectory speed vs traj angle
figure('Name','Binned trajectory speed vs traj angle')
trajAngle = motionInfo.traj_angle_lim(:);
trajAngle(trajAngle<0) = trajAngle(trajAngle<0) + 360;
trajSpeed = motionInfo.traj_speed(:);
unlikelyInds = find((trajAngle> 225) & (trajAngle < 315));
trajAngle(unlikelyInds) = [];
trajSpeed(unlikelyInds) = [];
nanInds = find(isnan(trajSpeed));
trajSpeed(nanInds) = [];
trajAngle(nanInds) = [];

binAngles = 0:360/30:360;

binnedSpeeds = cell(size(binAngles));
binnedAngles = binnedSpeeds;
numberInBin = zeros(size(binAngles));
binCtr = nan*numberInBin;
binMean = nan*numberInBin;
binStd = nan*numberInBin;
for jj = 1:length(binAngles)-1
    binInds = find((trajAngle >= binAngles(jj)) & (trajAngle < binAngles(jj+1)));
    if ~isempty(binInds)
    binnedSpeeds{jj} = trajSpeed(binInds);
    binnedAngles{jj} = trajAngle(binInds);
    numberInBin(jj) = numel(binInds);
    binCtr(jj) = mean([binAngles(jj) binAngles(jj+1)]);
    binMean(jj) = mean(trajSpeed(binInds));
    binStd(jj) = mean(trajSpeed(binInds));
    else       
        binnedSpeeds{jj} = nan;
        binnedAngles{jj}  = nan;
    end
end
binCtr = [binCtr(:); binCtr(1)];
binMean = [binMean(:); binMean(1)]
unlikelyInds = find((binAngles > 225) & (binAngles<315));
nanInds = find(isnan(binCtr));
remInds = setdiff(nanInds,unlikelyInds);
oneInds = find(numberInBin==1);
remInds = union(remInds,oneInds);
binCtr(remInds) = [];
binMean(remInds) = [];
 polar(binCtr*pi/180, binMean,'.-')
 
 binCtr_neg = binCtr;
 binCtr_neg(binCtr>180) = binCtr_neg(binCtr>180)-360;
 rInds = find(binCtr_neg>-90 & binCtr_neg<90);
 lInds = find(binCtr>90 & binCtr<270);
 

% figure('Name','Binned traj speed vs traj angle')
% polar(binCtr(rInds)*pi/180,binMean(rInds),'.-')
% hold on
% polar(mod(binCtr(lInds)+180,360)*pi/180,binMean(lInds),'ro-')
 

%% Binned trajectory vel vs traj angle
figure('Name','Binned trajectory vel vs traj angle')
trajAngle = motionInfo.traj_angle_lim(:);
trajAngle(trajAngle<0) = trajAngle(trajAngle<0) + 360;
trajSpeed = motionInfo.traj_vel(:);
unlikelyInds = find((trajAngle> 225) & (trajAngle < 315));
trajAngle(unlikelyInds) = [];
trajSpeed(unlikelyInds) = [];
nanInds = find(isnan(trajSpeed));
trajSpeed(nanInds) = [];
trajAngle(nanInds) = [];

binAngles = 0:360/30:360;

binnedSpeeds = cell(size(binAngles));
binnedAngles = binnedSpeeds;
numberInBin = zeros(size(binAngles));
binCtr = nan*numberInBin;
binMean = nan*numberInBin;
binStd = nan*numberInBin;
for jj = 1:length(binAngles)-1
    binInds = find((trajAngle >= binAngles(jj)) & (trajAngle < binAngles(jj+1)));
    if ~isempty(binInds)
    binnedSpeeds{jj} = trajSpeed(binInds);
    binnedAngles{jj} = trajAngle(binInds);
    numberInBin(jj) = numel(binInds);
    binCtr(jj) = mean([binAngles(jj) binAngles(jj+1)]);
    binMean(jj) = mean(trajSpeed(binInds));
    binStd(jj) = mean(trajSpeed(binInds));
    else       
        binnedSpeeds{jj} = nan;
        binnedAngles{jj}  = nan;
    end
end
polar(binCtr*pi/180, binMean,'.-')