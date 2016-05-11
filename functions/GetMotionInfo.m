function motionInfo = GetMotionInfo(fishPos,orientation,imgLen,varargin)
%GetMotionInfo - Returns some info about fish motion
% motionInfo = GetMotionInfo(fishPos, orientation, imgLen)
% motionInfo = GetMotionInfo(...,motionThr);
% motionInto = GetMotionInfo(...,trajPlot)
% Inputs:
% fishPos - T x 2 matrix where T is the total number of time points (img
%   frames) and the 1st and the 2nd cols contain the x and y coordinates of the fish
% orientation - A vector containing the orientation of the fish at each
%   time point
% imgLen - Length of each image, i.e., size(IM,1)
% motionThr - The fish has to move by this many pixels to be considered in
%   motion
% trajPlot - Trajectory plot boolean
% Outputs:
% motionInfo - Structure variable containing swimming info extracted from
%   trajectories.
%   .motionFrames - Frame indices where fish was in motion
%   .swimStartFrames - Frames where fish commenced a swim episode
%   .ds - ?? (later)
%   .epInds - Frames corresponding to swim episodes
%   .traj - Unadjusted traj
%   .traj_trans - Initial x,y position corrected traj
%   .traj_angle - Initial orientation corrected traj
%   .traj_adj - Initial x,y position and orientaion corrected traj
%   .traj_angle - Angle of vector connecting starting and end point of traj
%   .traj_speed - Mean speed over 1st three frames (or fewer) of traj
%   .traj_vel - Mean vel computed over 1st ...
%   .traj_angVel - Mean angular vel computer over 1st ...
%   .traj_angle_lim - Limited traj angle computed over 1st ...
%   .turnInfo - ?? (later)
%   .orientInfo - ?? (later)
% 
% Avinash Pujala, HHMI, 2016

trajPlot = 0; % 1 results in plotting of trajectories
motionThr = 5;
lpf = 50; % Low pass filter for timeseries (only if sampling rate greater than 100 fps)
frameRate = 500;

if nargin ==4
    motionThr = varargin{1};
elseif nargin ==5
    motionThr = varargin{1};
    trajPlot = varargin{2};
end
[motionFrames, swimStartFrames] = GetMotionFrames(fishPos,motionThr);

%## Distance info
D = @(v)sqrt(sum(diff(v,[],1).^2,2));
x = fishPos(swimStartFrames,1);
y = fishPos(swimStartFrames,2);
% dS = sqrt(diff(x).^2 + diff(y).^2);
dS = D(fishPos(swimStartFrames,:));
dS_all = D(fishPos);

%## Orientation info
if (size(orientation,1)==2) && (size(orientation,1) < size(orientation,2))
    orientation = orientation';
elseif (size(orientation,1) < size(orientation,2))
    orientation = orientation';
end
or = orientation(swimStartFrames,1);
or = mod(or+180,360);
or = 360-or;
ker = gausswin(6); ker = ker/sum(ker);
dOr = conv2((-DiffOrientation(orientation(:,1)))*180/pi, ker(:),'same');


%## Head curvature info
if size(orientation,2)>1
    curv = GetCurvInfo(orientation);
    curv = [curv, sum(curv,2)];
    curv = fix(curv/8)*8;
    if frameRate < 100
        for jj = 1:size(curv,2)
            curv(:,jj) = conv2(curv(:,jj),ker(:),'same');
        end
    else
         for jj = 1:size(curv,2)
            curv(:,jj) = chebfilt(curv(:,jj),1/500,lpf,'low');
        end
    end
else
    curv = NaN;
end


%## Motion vec info
motionVecInfo = GetMotionVecInfo(fishPos);

traj = {};
traj_trans = {};
traj_trans_smooth = {};
traj_adj = {};
epInds = {};
traj_angle = zeros(size(swimStartFrames));
[traj_angle_lim,speed,vel,angVel] = deal(traj_angle);

if trajPlot
    figure
    %     axis image
end
for jj = 2:length(swimStartFrames)-1
    epInds{jj} = swimStartFrames(jj):swimStartFrames(jj+1)-1;
    xx = fishPos(epInds{jj},1);
    yy  = imgLen-fishPos(epInds{jj},2);
    
    %## Compute starting position and orientation adjusted trajectory
    traj{jj} = [xx(:), yy(:)];
    traj_trans{jj} = [xx(:)-xx(1), yy(:)-yy(1)];
    [theta,rho] = cart2pol(traj_trans{jj}(:,1),traj_trans{jj}(:,2));
    theta = theta-(or(jj)*pi/180)-pi/2;
    [xx_or, yy_or] = pol2cart(theta,rho);
    traj_adj{jj} = [xx_or(:) yy_or(:)];
    
    %## Traj angle (angle of vector connecting start and end pt of traj)
    c = traj_adj{jj}(end,1) + traj_adj{jj}(end,2)*1i;
    if abs(c)== 0
        traj_angle(jj) = nan;
    else
        traj_angle(jj) = angle(c)*180/pi;
    end
    %     traj_angle(jj)  = GetTrajAngle(traj_adj);
    
    %## Traj speed (Mean of speed in first 3 frames or fewer; since frame rate is
    %## fixed, using dist as proxy for speed), vel, angular vel
    nPts = min(size(traj_adj{jj},1),4);
    blah = traj_adj{jj}(1:nPts,:);
    speed(jj) = sum(sqrt(sum((diff(blah)).^2,2)),1)/(nPts-1);
    vel(jj) = abs(blah(end,1) + blah(end,2)*1i)/(nPts-1);
    angVel(jj) = (angle(blah(end,1) + blah(end,2)*1i)*180/pi)/(nPts-1);
    if length(blah)>2
        traj_angle_lim(jj) = angle(blah(2,1) + blah(2,2)*1i)*180/pi;
    end
    
    
    %## Plot trajectories if specified
    if trajPlot
        hold on
        plot(traj_adj{jj}(:,1),traj_adj{jj}(:,2),'.-','color',rand(1,3)), drawnow
        axis image
        %         pause(0.1)
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
motionInfo.dS_all = [0; dS_all(:)];
motionInfo.or = or;
motionInfo.dOr = dOr;
motionInfo.curv = curv;
motionInfo.motionVecs = motionVecInfo;
motionInfo.epInds = epInds;
motionInfo.traj = traj;
motionInfo.traj_trans = traj_trans;
% motionInfo.traj_trans_smooth  = traj_trans_smooth;
motionInfo.traj_adj = traj_adj;
motionInfo.traj_angle = traj_angle(2:end);
motionInfo.traj_angle_lim = traj_angle_lim(2:end);
motionInfo.traj_speed = speed(2:end);
motionInfo.traj_vel = vel(2:end);
motionInfo.traj_angVel = angVel(2:end);
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

% motFrames = find((dS1 >= motionThr) & (dS2 >=2*motionThr));
motFrames = find(dS1 >= motionThr);

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
if (size(orientation,1)==2) && (size(orientation,1) < size(orientation,2))
    orientation = orientation';
end
orient_mov = orientation(motionFrames,1);
orient_mov = orient_mov(:);
[oX,oY] = pol2cart(orient_mov*pi/180,1);
orVec = oX + oY*1i;
dOr = angle(orVec(1:end-1,1).*conj(orVec(2:end,1)))*180/pi;
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

function curv = GetCurvInfo(orientation)
if (size(orientation,1)==2) && (size(orientation,1)<size(orientation,2))
    orientation = orientation';
end
curv = nan(size(orientation,1),size(orientation,2)-1);
toVec = @(x)pol2cart(x*pi/180,ones(size(x)));
for seg = 1:size(curv,2)
    [x1,y1] = toVec(orientation(:,seg));
    [x2,y2] = toVec(orientation(:,seg+1));
    curv(:,seg)= angle((x1(:) + y1(:)*1i).* conj(x2(:) + y2(:)*1i));   
end
curv = -curv*180/pi;
end

function traj_angle = GetTrajAngle(traj_adj)
c = traj_adj(end,1) + traj_adj(end,2)*1i;
if abs(c)== 0
    traj_angle(jj) = nan;
else
    traj_angle(jj) = angle(c)*180/pi;
end
end

function motionVecInfo = GetMotionVecInfo(fishPos)
% Gets angle and magnitude of motion vectors

[x,y] =  deal(fishPos(:,1),fishPos(:,2));
R = @(x)(x + 1*(rand(size(x))-0.5));
[x,y] = deal(R(x), R(y));
[dx,dy] = deal(diff(x),diff(y));
motionVecs = dx + dy*1i;
A = @(z)(angle(z(2:end).*conj(z(1:end-1))))*180/pi;
motionAngles = A(motionVecs);
motionLengths = abs(motionVecs(2:end));
motionVecInfo = [motionAngles, motionLengths];
end





