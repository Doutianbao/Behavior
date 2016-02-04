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

motFrames = find((dS1 >= motionThr) & (dS2 >= 1.25*motionThr));

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

