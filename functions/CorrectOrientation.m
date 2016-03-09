
function orientation = CorrectOrientation(orientation,varargin)
% CorrectOrientation - Returns corrected fish orientation based on the
%   assumption that fish orientation cannot change by more than the
%   inputted steepAngle from one from to the immediately next.
% orientation =  CorrectOrientation(orientation, steepAngle)
%% dO vectorially
% orientation = (360-fish.orientation);
if nargin == 1
    steepAngle = 75;
elseif nargin == 2
    steepAngle = varargin{1};
end

orientation(1) = orientation(2);
[orX,orY] = pol2cart(orientation*pi/180,1);
clear i
orVec = orX + orY*i;
dOr = angle(orVec(1:end-1,:).*conj(orVec(2:end,:)))*180/pi;
jumpInds = find(abs(dOr)>=steepAngle);
for jj = 1:length(dOr)
    th = angle(orVec(jj,:).*conj(orVec(jj+1,:)))*180/pi;
    if abs(th) >= steepAngle
%         orientation(jj+1) = mod(180+orientation(jj+1),360);
       orientation(jj+1) = orientation(jj);
       orVec(jj+1,:) = orVec(jj,:);
    end
end

% Reobtain dOr
[orX,orY] = pol2cart(orientation*pi/180,1);
orVec = orX + orY*i;
dOr = angle(orVec(1:end-1,:).*conj(orVec(2:end,:)))*180/pi;
jumpInds = find(abs(dOr)>=130);
