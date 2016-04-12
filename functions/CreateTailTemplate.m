function varargout = CreateTailTemplate(base,height,varargin)
%CreateTailTemplate - Creates a triangular template for matching with
%   fish's tail using specified base and height
% [template, templateInds] = CreateTailTemplate(base,height);
% [template, templateInds] = CreateTailTemplate(..., grad);
% 
% Inputs:
% base - Base of the trianglular template in pxls
% height - Height of the triangular template in pxls
% Outputs:
% template - Binary image that tightly contains the template shape
% templateInds - N x 2 matrix, where 1st and 2nd col are row and col inds
%   respectively
% grad - Intensity gradient to apply to tail template. If grad < 0, then
%   template intensity from base of triangle to tip drops with slope equal
%   to grad
grad = -0.15;
if nargin==3;
    grad = varargin{1};
end
x_sub = round(base/2)*(0:1);
y = 0:height;
y_sub = linspace(0,height,length(x_sub));
x = interp1(y_sub,x_sub,y,'spline');
side1 = [sort(-x(:)),y(:)];
side2 = [x(:), flipud(y(:))];
side3 = [[side1(:,1); side2(:,1)],zeros(size([side1(:,1); side2(:,1)]))];
tempOutline = [side1; side2; side3];
img = zeros(2*(length(min(-x):max(x)))+1, 2*(max(y))+1);
getCtr = @(x)[round(size(x,1)/2 + 0.4999) round(size(x,2)/2 + 0.4999)];
imgCtr = getCtr(img);
tempOutline(:,1) = tempOutline(:,1) + imgCtr(1);
tempOutline(:,2) = tempOutline(:,2) + imgCtr(2);
[x_round, y_round] =  deal(round(tempOutline(:,1)), round(tempOutline(:,2)));
outlineInds = unique(sub2ind(size(img),x_round,y_round));
img(outlineInds)=1;
img = img(imgCtr(1)-ceil(base/2):imgCtr(1)+ceil(base/2), imgCtr(2):end);
[r,c] = find(ones(size(img)));
[x,y] = find(img); 
img(inpolygon(r,c,x,y)) = 1;
t = 1:size(img,2);
y = grad*t;
[X,~] = meshgrid(-grad*t,y);
X = X(1:size(img),:);
img = img.*X;
varargout{1} = img;
[r,c] = find(img);
varargout{2} = [r,c];
end