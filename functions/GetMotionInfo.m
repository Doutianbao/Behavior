function motionInfo = GetMotionInfo(fishPos,orientation,imgLen)
%GetMotionInfo - Returns some info about fish motion
% motionInfo = GetMotionInfo(fishPos, orientation, imgLen)
% Inputs:
% fishPos - T x 2 matrix where T is the total number of time points (img
%   frames) and the 1st and the 2nd cols contain the x and y coordinates of the fish
% orientation - A vector containing the orientation of the fish at each
%   time point
% imgLen - Length of each image, i.e., size(IM,1)

trajPlot = 0;
motionThr = 4;

[motionFrames, swimStartFrames] = GetMotionFrames(fishPos,motionThr);

x = fishPos(swimStartFrames,1);
y = fishPos(swimStartFrames,2);
dS = sqrt(diff(x).^2 + diff(y).^2);

or = orientation(swimStartFrames);
or = mod(or+180,360);
or = 360-or;

traj = {};
traj_trans = {};
traj_trans_smooth = {};
traj_adj = {};
epInds = {};
ker  = gausswin(3); ker = ker/sum(ker);
traj_angle = zeros(size(swimStartFrames));
if trajPlot
    figure
%     axis image
end
for jj = 2:length(swimStartFrames)-1
    epInds{jj} = swimStartFrames(jj):swimStartFrames(jj+1)-1;
    xx = fishPos(epInds{jj},1);    
    yy  = imgLen-fishPos(epInds{jj},2);
    traj{jj} = [xx(:), yy(:)];
    traj_trans{jj} = [xx(:)-xx(1), yy(:)-yy(1)];
    [theta,rho] = cart2pol(traj_trans{jj}(:,1),traj_trans{jj}(:,2));
    theta = theta-(or(jj)*pi/180)-pi/2;
    [xx_or, yy_or] = pol2cart(theta,rho);    
    traj_adj{jj} = [xx_or(:) yy_or(:)];
    c = traj_adj{jj}(end,1) + traj_adj{jj}(end,2)*1i;
    if abs(c)== 0
        traj_angle(jj) = nan;
    else
        traj_angle(jj) = angle(c)*180/pi;
    end
    
    if trajPlot
        hold on
        plot(traj_adj{jj}(:,1),traj_adj{jj}(:,2),'.-','color',rand(1,3)), drawnow
%         pause(0.5)
        xlim([-inf inf])
        ylim([-inf inf])
        shg
    end
end
traj_angle(traj_angle < 0) = traj_angle(traj_angle<0) + 360;

%### Extracting turn info based on swim path and fish orientation
turnInfo = getTurnInfo(fishPos, motionFrames);
orientInfo = getOrInfo(orientation,motionFrames);

%### Creating motionInfo struct variable
motionInfo.motionFrames = motionFrames;
motionInfo.swimStartFrames = swimStartFrames;
motionInfo.dS = dS;
motionInfo.or = or;
motionInfo.epInds = epInds;
motionInfo.traj = traj;
motionInfo.traj_trans = traj_trans;
% motionInfo.traj_trans_smooth  = traj_trans_smooth;
motionInfo.traj_adj = traj_adj;
motionInfo.traj_angle = traj_angle;
motionInfo.turnInfo = turnInfo;
motionInfo.orientInfo = orientInfo;
end

function varargout = GetMotionFrames(fishPos, motionThr)
%GetMotionFrames Return indices of frames wherein fish is in motion
% motionFrames = GetMotionFrames(fishPos,motionThr);
% motionFrames, swimStartFrames = GetMotionFrames(...);
% Inputs:
% fishPos - T x 2 matrix where T is number of time points. 1st and 2nd cols
%   contain x and y coordinates of the fish
% motionThr - minimum distance that must be traversed between subsequent
%   frames for it to be considered motion as opposed to drift
% Outputs:
% motionFrames - All the frames in which the fish is in motion
% swimStartFrames - Subset of motion Frames corresponding to onset of swim
%   episode

x = fishPos(:,1);
y = fishPos(:,2);

dx1 = x(2:end)-x(1:end-1);
dy1 = y(2:end)-y(1:end-1);
dS1 = sqrt(dx1.^2 + dy1.^2);
dS1 = dS1(1:end-1);


dx2 = x(3:end)-x(1:end-2);
dy2 = y(3:end)-y(1:end-2);
dS2 = sqrt(dx2.^2 + dy2.^2);

motFrames = find((dS1 >= motionThr) & (dS2 >=2*motionThr));

startInds = [];
initInd = motFrames(1);
for jj = 2:length(motFrames)
    if motFrames(jj)-motFrames(jj-1) > 1
        initInd = motFrames(jj-1);
        nextInd = motFrames(jj);
        startInds = [startInds; initInd;nextInd]; 
    end
end
startInds = unique(startInds);

varargout{1} = motFrames;
varargout{2} = startInds;

end


function turnInfo = getTurnInfo(fishPos,motionFrames)
turnAngleThr = 5;

fishPos_motion = fishPos(motionFrames,:);
motionVecs = diff(fishPos_motion);
motionVecs = motionVecs(:,1) + motionVecs(:,2)*1i;
dTh = angle(motionVecs(1:end-1).*conj(motionVecs(2:end)));
dTh_motion = [0; dTh(:)]*180/pi;

rightInds = motionFrames(find(dTh_motion < -turnAngleThr)+1);
leftInds = motionFrames(find(dTh_motion > turnAngleThr)+1);
straightInds = motionFrames(find(abs(dTh_motion) <=5)+1);

turnInfo = {};
turnInfo{1} = leftInds;
turnInfo{2} = rightInds;
turnInfo{3} = straightInds;
turnInfo{4} = dTh_motion;
turnInfo{5} = {'Left Inds', 'Right Inds', 'Straight Inds','Turn angles'};

end

function orientInfo = getOrInfo(orientation,motionFrames)
orient_mov = orientation(motionFrames); 
[oX,oY] = pol2cart(orient_mov*pi/180,1);
orVec = oX + oY*1i;
dOr = angle(orVec(1:end-1,:).*conj(orVec(2:end,:)))*180/pi;
straightInds = find(dOr > -5 & dOr < 5);
leftInds = find(dOr >5);
rightInds = find(dOr <-5);

orientInfo = {};
orientInfo{1} = leftInds;
orientInfo{2} = rightInds;
orientInfo{3} = straightInds;
orientInfo{4} = dOr;
orientInfo{5} = {'leftInds','rightInds','straightInds','dOr'};
end
