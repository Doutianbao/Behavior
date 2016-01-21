function dS = deltaDistance(xyCoord)
%deltaDistance - When given x and y coordinates, computes distance traveled
%   from previous point
%   dS = deltaDistace(xyCoordinates);
%   Inputs:
%   xyCoordinates - n-by-2 matrix where n is number of coordinate points, and the 2 columns
%        are x and y coordinates respectively

x = xyCoord(:,1);
y = xyCoord(:,2);
dS = zeros(size(x,1)-1,1);

dx = diff(x);
dy  = diff(y);
dS = (dx.^2 + dy.^2).^0.5;
end

