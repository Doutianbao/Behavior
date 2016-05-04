function dOr = DiffOrientation(orientation)
% DiffOrientation - Given a vector of angles (in degrees) which corresponds to fish
%   orientation returns the derivative of those angles
% dOr = DiffOrientation(orienttion)

orientation = orientation(:)*pi/180;
[x,y] = pol2cart(orientation,ones(size(orientation)));
c = x + 1i*y;
dOr = angle(c(2:end).*conj(c(1:end-1)));

end

