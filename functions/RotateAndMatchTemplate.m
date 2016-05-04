
function varargout = RotateAndMatchTemplate(im,varargin)
% RotateAndMatchTemplate - Rotates a line of certain length around fish head centroid
%   and reconstructs fish pixels that way
% matchVec = RotateAndMatchTemplate(im,varargin)
% [matchVec, corrVec] = RotateAndMatchTemplate(im,varargin)
% varargin =  startPt,dTh,template
% Inputs:
% im - Image containing fish (Please note, the img must be such that eyes are bright not dim)
% fishPos - x,y coordinates of head centroid. If not specified or [], then
%   automatically determines fish pos
% dTh - Angles by which to rotate line (default = 5)
% template  - Custom template or that created by CreateTailTemplate.
% Outputs:
% matchVec - A vector of values that indicate degree of match after rotation
%   of template at various angles
% corrVec - A vector of values that indicate correlateion ....
% thetas - Angles (in degrees), by which template is rotated.

dTh = 5;
base = 6;
height = 25;
grad = -0.15;
if nargin ==1
    fishPos = GetFishPos(im,30);
    [T,~] = CreateTailTemplate(base,height,grad);
elseif nargin == 2
    fishPos = varargin{1};
    [T,~] = CreateTailTemplate(base,height,grad);
elseif nargin ==3
    fishPos = varargin{1};
    dTh = varargin{2};
    [T,~] = CreateTailTemplate(base,height,grad);
elseif nargin ==4
    fishPos = varargin{1};
    dTh = varargin{2};
    T = varargin{3};
end

if isempty(fishPos)
    fishPos = GetFishPos(im,30);
end
if isempty(dTh)
    dTh  = 5;
end

fishPos = fliplr(fishPos); % The form in which GetFishPos outputs fishPos is not in row, col coordinates.
[r,c] = find(T);
imCtr = [round(size(im,1)/2 + 0.4999),round(size(im,1)/2 + 0.4999)];
x = round(r - (size(T,1)/2) + imCtr(1));
y = c + imCtr(2);
overX = find(x<1 | x>size(im,1));
overY = find(y<2 | y>size(im,2));
x(overX)=[];
y(overY)=[];
r(overX)=[];
c(overY) =[];
T_trunc = T(min(r):max(r), min(c):max(c));
thetas = 0:dTh:360;
C = nan(size(thetas));
M = C;
im = TranslateImage(im,fishPos);
for rot = 1:length(thetas)
    img_rot = imrotate(im,thetas(rot),'crop');
    img_trunc = img_rot(min(x):max(x), min(y):max(y));
    C(rot) = corr2(img_trunc,T_trunc);
    blah = (img_trunc.*T_trunc);
    M(rot) = sum(blah(:))/sum(T_trunc(:));
end

varargout{1} = M;
varargout{2} = C;
varargout{3} = thetas;
end