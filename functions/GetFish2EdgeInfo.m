function varargout = GetFish2EdgeInfo(fishPos,orientation,edgeInds)
%GetFish2EdgeInfo Returns the distance between the fish and it nearest edge
%   as well as the angle between the fish's orientation and the tangent at
%   this edge point.
% [S, O] = GetFish2EdgeInfo(fishPos, orientation, edgeInds);
% 
% Inputs:
% fishPos - T x 2 matrix, where 1st and 2nd col are the x, and y
%   coordinates of the fish's head centroid for each of the T time points
% orientation - A vector of size T, where each entry is the fish's
%   orientation for the corresponding time point
% edgeInds - Indices of a circle fit to the edge of the arena
% 
% Outputs:
% S - Vector of fish's distance from nearest edge index
% O - Vector of the difference in angle between the fish's orientation and the tangenst at the nearest edge point
% E - Vector of edge inds that were closest to the fish
% 
% Avinash Pujala, HHMI, 2016

S = nan(size(fishPos,1),1);
O = S;
O2 = S;
E = nan(size(fishPos));
E2 = E;
eInds = edgeInds;
eInds = [eInds(:,1)-mean(eInds(:,1)), eInds(:,2)-mean(eInds(:,2))];
[th, rho] = cart2pol(eInds(:,1), eInds(:,2));
th = th*180/pi;
% th = PhaseToTangent(th*180/pi);
% th = mod(th+90, 90);
M = @(v)(sqrt(v(1)^2 + v(2)^2));
A = @(v1,v2)acos((dot(v1,v2)/(M(v1)*M(v2))))*180/pi;
disp('Getting dist and angle to nearest edge point...')
dispChunk = round(size(fishPos,1)/5);
tic
for jj = 1:size(fishPos,1)
    eInds = edgeInds;
    fp = fishPos(jj,:);
    [ind1, S(jj)] = Dist2Nearest(fp,eInds);
    v = eInds(ind1,:);
    tVec = v-fp;
    [theta,~] = cart2pol(tVec(1),tVec(2));
    [tVec(1), tVec(2)] = pol2cart((theta+pi/2),1);
    
    E(jj,:) = v;
  
    or = mod(orientation(jj)+180,360);

    [x,y] = pol2cart(or*pi/180,1);
%     [v1(1),v1(2)] = pol2cart(th(ind1)*pi/180,1);
    v2 = [x,y];
%     O(jj) = (180-A(v1,v2));
    O(jj) = 90-PhaseToTangent(A(v2,tVec));
    O2(jj) = angle(tVec(1) + tVec(2)*1i)*180/pi;
    
    if mod(jj,dispChunk)==0
        disp(['Frame # ' num2str(jj)])
    end
end
toc
varargout{1} = S;
varargout{2} = O;
varargout{3} = E;
varargout{4} = O2;
end

function [nearestInd,dist] = Dist2Nearest(v,V)
% Finds the 2D coordinates in V that are nearest to v, as well as the
%   distance
if size(v,1)==2 && size(v,2)~=2
    v= v';
end
if size(V,1)==2 && size(V,2)~=2
    V= V';
end
[~,nearestInd] = min(sqrt(sum((repmat(v,size(V,1),1)-V).^2,2)));

dist = sqrt(sum((V(nearestInd,:)-v).^2));

end

function c = Angle2Complex(angle)
% Given an angle in degrees, returns a unit length complex # pointing in the
% direction of the angle
[x,y] = pol2cart(angle*pi/180,1);
c = x + y*1i;
end

function theta = PhaseToTangent(theta)
theta = asin(sin(theta.*pi/180)).*180/pi; 
theta = abs(theta); 
theta = 90-theta;
end

