function motionInfo = GetMotionInfo(fishPos,orientation,swimStartFrames,imgLen)
%GetMotionInfo - Returns some info about fish motion

trajPlot = 1;

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
if trajPlot
    figure
    axis image
end
for jj = 2:length(swimStartFrames)-1
    epInds{jj} = swimStartFrames(jj):swimStartFrames(jj+1)-1;
    xx = fishPos(epInds{jj},1);    
    yy  = imgLen-fishPos(epInds{jj},2);
    traj{jj} = [xx(:), yy(:)];
    traj_trans{jj} = [xx(:)-xx(1), yy(:)-yy(1)];
%     xx = conv2(xx(:)-xx(1),ker(:),'same');
%     yy = conv2(yy(:)-yy(1),ker(:),'same');  
%     traj_trans_smooth{jj} = [xx(:), yy(:)];
    [theta,rho] = cart2pol(traj_trans{jj}(:,1),traj_trans{jj}(:,2));
    theta = theta-(or(jj)*pi/180)-pi/2;
    [xx_or, yy_or] = pol2cart(theta,rho);    
    traj_adj{jj} = [xx_or(:) yy_or(:)];
    if trajPlot
%         plot(traj{jj}(:,1), traj{jj}(:,2),'.-')
        hold on
        plot(traj_adj{jj}(:,1),traj{jj}(:,2),'r.-')
        xlim([-200 200])
        ylim([-200 200])
    end
end
motionInfo.dS = dS;
motionInfo.or = or;
motionInfo.epInds = epInds;
motionInfo.traj = traj;
motionInfo.traj_trans = traj_trans;
% motionInfo.traj_trans_smooth  = traj_trans_smooth;
motionInfo.traj_adj = traj_adj;

% orientation = orientation;
% orient_mov = orientation(startInds); 
% [oX,oY] = pol2cart(orient_mov*pi/180,1);
% clear i
% orVec = oX + oY*i;
% dOr = angle(orVec(1:end-1,:).*conj(orVec(2:end,:)))*180/pi;
% straightInds = find(dOr > -5 & dOr < 5);
% leftInds = find(dOr >5);
% rightInds = find(dOr <-5);


end



