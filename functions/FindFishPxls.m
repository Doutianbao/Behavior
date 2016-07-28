function varargout = FindFishPxls(img,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


n = 5;
if nargin > 1
    n = varargin{1};
end

thr = multithresh(img,n);

img_quant = imquantize(img,thresh);

levels = unique(img_quant);
for level = 1:numel(levels)
    [x,y] = find(img_quant==levels(level));
    count =0;
    [xx,yy] = find()
end

end

