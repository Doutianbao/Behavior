function varargout = GetFishPxls(im, fishPos, varargin)
%GetFishPxls When given a background subtracted image, and fish coordinates
%    extracts fish only pxls
%
%  fishPxlImg = GetFishPxls(im,fishPos)
% [fishPxlImg,fishPxls] = GetFishPxls(im,fishPos)
%       ...             = GetFishPxls(im,fishPos,intensityThr,searchRadius)
% Inputs:
% im  - Image in which to find fish (it is assumed that image has been 
%       processed such that fish pixels are brighter)
% fishPos - Position of fish in cartesian space (not image space)
% intensityThr - Intensity over which to threshold to find fish pixels. If
%   instensityThr = [], or it is not specified then finds threshold using
%   Otsu's method (multithresh MATLAB function)


iThr_def = multithresh(im,3);
iThr = iThr_def(2);
dThr = 105;
if nargin ==3
    iThr = varargin{1};    
elseif nargin == 4
    iThr = varargin{1};
    dThr = varargin{2};
end
if isempty(iThr)
    iThr = iThr_def(2);
end

ker = gausswin(5,2);
ker = ker(:)*ker(:)';
im(im<iThr)=0;
[r,c] = find(im);
x = c;
y = r;
X = nan*zeros(size(im));
Y  = X;
inds = sub2ind(size(im),r,c);
X(inds) = x;
Y(inds) = y;
S = sqrt((X-fishPos(1)).^2 + (Y-fishPos(2)).^2);
S(S>dThr)=0;
inds = find(S);

% S_conv = conv2(S,ker,'same');
% B = S_conv.*S;
se = strel('disk',4);
S(isnan(S))=0;
S_open = imopen(S,se);
B = S_open.*S;
B(B>0) = max(B(:));
B(isnan(B))=0;

if nargout < 2
    varargout{1} = B;
else
    varargout{2} = find(B);    
end

end

