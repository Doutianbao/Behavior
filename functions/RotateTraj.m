function traj_rot = RotateTraj(traj, theta)
%RotateTraj - Given a trajectory (N X 2), returns the trajectory after
%   rotation by specified angle
% traj_rot = RotateTraj(traj,theta);
% 
% Avinash Pujala, Koyama lab/HHMI, 2016

if size(traj,1)==2 && size(traj,2)~=2
    traj = traj';
end
theta = theta*(pi/180); % Convert to radians
T_rot =[[cos(theta), -sin(theta)]; [sin(theta) cos(theta)]]; % Rotation transform
traj_rot = (T_rot*traj')';

end

