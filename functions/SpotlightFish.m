
function spotlightImg = SpotlightFish(im,varargin)
% SpotlightFish - Returns an image where pixels outside a radius from fish
%   head centroid are set to 0
% spotLightImg = SpotlightFish(im)
% radialImg = SpotlightFish(im,fishPos)
% radialImg = RadialFish(...,radius)
% Inputs:
% im - Image containing fish
% fishPos - x,y coordinates of head centroid. If not specified or [], then
%   automatically determines fish pos
% radius - Spotlight radius (default = 70)


radius = 70;
if nargin ==1
    fishPos = GetFishPos(im,30);
elseif nargin ==2
    fishPos = varargin{1};
elseif nargin == 3
    fishPos = varargin{1};
    radius = varargin{2};    
end
if isempty(fishPos)
    fishPos = GetFishPos(im,30);
end
[row,col] = find(ones(size(im)));
S = sqrt(sum([col-fishPos(1) row-fishPos(2)].^2,2));
outPxls = find(S>=radius);

spotlightImg = im;
spotlightImg(outPxls)=0;


end