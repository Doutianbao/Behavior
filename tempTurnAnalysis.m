
%% Swim analysis
x = tracexy(:,1);
y = tracexy(:,2);
dx = diff(x);
dy = diff(y);
dS = sqrt(dx.^2 + dy.^2);
noMovInds  = find(dS < 3);
movInds = setdiff(1:length(dx),noMovInds);
startInds = [];
initInd = movInds(1);
for jj = 2:length(movInds)
    if movInds(jj)-movInds(jj-1) > 1
        initInd = movInds(jj-1);
        nextInd = movInds(jj);
        startInds = [startInds; initInd;nextInd]; 
    end
end
startInds = unique(startInds);
x_mov = x(startInds);
y_mov = y(startInds);
dS_mov = sqrt(diff(x_mov).^2 + diff(y_mov).^2);

orientation = orientation;
orient_mov = orientation(startInds); 
[oX,oY] = pol2cart(orient_mov*pi/180,1);
clear i
orVec = oX + oY*i;
dOr = angle(orVec(1:end-1,:).*conj(orVec(2:end,:)))*180/pi;
straightInds = find(dOr > -5 & dOr < 5);
leftInds = find(dOr >5);
rightInds = find(dOr <-5);
dOr_str = dOr(straightInds);
dOr_left = dOr(leftInds);
dOr_right = -dOr(rightInds);
[orCount_left,orVals] = hist(dOr_left,[10:10:360]);
[orCount_right,~] = hist(dOr_right,[10:10:360]);
orMat = [orCount_left(:),orCount_right(:)];
figure('Name','Turn histogram')
bar(orVals,orMat,1)
xlim([0 130])
box off
legend('Left','Right')
set(gca,'tickdir','out','xtick',orVals)
xlabel('Turn angle')
ylabel('Count')
title('Turn angle histogram')

%% Swim distance by turn
dS_left = dS_mov(leftInds);
dS_right = dS_mov(rightInds);
dS_straight = dS_mov(straightInds);

[dslCount,dslVals] = hist(dS_left,20);
[dsrCount,dsrVals] = hist(dS_right,20);
[dssCount,dssVals] = hist(dS_straight,20);
ker = gausswin(5);
ker = ker/sum(ker);
dslCount = conv2(dslCount(:),ker(:),'same');
dsrCount = conv2(dsrCount(:),ker(:),'same');
dssCount = conv2(dssCount(:),ker(:),'same');
figure('Name','Swim dist histogram')
plot(dslVals,dslCount,'linewidth',2)
hold on
plot(dsrVals,dsrCount,'r','linewidth',2)
plot(dssVals,dssCount,'k','linewidth',2)
xlim([-inf inf])
ylim([-inf inf])
box off
set(gca,'tickdir','out')
xlabel('Swim distance')
ylabel('Count')
legend('Left','Right','Straight')
title('Swim distance by turn distribution')

%% Velocity stuff
timeVec = (0:length(tracexy)-1)*(1/30); % Since frame rate = 30fps
timeVec_mov = timeVec(startInds);
dT_mov = diff(timeVec_mov);
vel_mov = dS_mov(:)./dT_mov(:);
vel_left = vel_mov(leftInds);
vel_right = vel_mov(rightInds);
vel_straight = vel_mov(straightInds);
[vlCount,vlVals] = hist(vel_left,50);
[vrCount,vrVals] = hist(vel_right,50);
[vsCount,vsVals] = hist(vel_straight,50);
vlCount = conv2(vlCount(:),ker(:),'same');
vrCount = conv2(vrCount(:),ker(:),'same');
vsCount = conv2(vsCount(:),ker(:),'same');
figure('Name','Velocity by turn distribution')
plot(vlVals,vlCount,'linewidth',2)
hold on
plot(vrVals,vrCount,'r','linewidth',2)
% plot(vsVals,vsCount,'k','linewidth',2)
box off
xlim([50 700])
ylim([-inf inf])
xlabel('Swim Vel (plx/sec)')
ylabel('Count')
title('Swim vel by turn distribution')
set(gca,'tickdir','out')
legend('Left','Right')
% legend('Left','Right','Straight')



break


