
function varargout = GetInnerPxls(im,periPxls)
% GetInnerPxls - Given an image and the indices of a perimeter, returns the
%   indices of pixels that are inside and outside the convex hull defined 
%   by the perimenter indices
% inPxls = GetInnerPxls(im, periPxls);
% [inPxls, outPxls] = GetInnerPxls(im,periPxls);
% Inputs:
% im - Input image
% periPxls - Indices of the perimeter of a shape within which to find the
%   inner pixels
% Outputs:
% inPxls - Inner pixels within the convex hull defined by the perimeter
% outPxls - Outer pixels    "           "           "
% Doesn't really work, and also there's MATLAB's inpolygon

if nargin < 2
    error('At least 2 inputs required!')
end
if all(size(periPxls)>1)
    if size(periPxls,1)==2
        periPxls = periPxls';
    end
    r = periPxls(:,1);
    c = periPxls(:,2);
    periPxls = sub2ind(size(im),r,c);
else
    [r,c] = ind2sub(size(im),periPxls);
end

[row,col] = find(ones(size(im)));
inPxls = [];
outPxls = [];
% refPt = [r(1) c(1)];
% refS = sum(sqrt(sum(([r(:) c(:)]- repmat(refPt,length(r),1)).^2,2)));
for jj = 1:length(row)
%     S = sum(sqrt(sum((repmat([row(jj) col(jj)],length(r),1) - [r(:) c(:)]).^2,2)));
%     if S < refS
%         inPxls = [inPxls;[row(jj) col(jj)]];
%     else
%         outPxls = [outPxls;[row(jj) col(jj)]];
%     end
  
end
inPxls = sub2ind(size(im), inPxls(:,1), inPxls(:,2));
outPxls = sub2ind(size(im),outPxls(:,1),outPxls(:,2));
outPxls = setdiff(outPxls,periPxls);

varargout{1} = inPxls;
varargout{2} = outPxls;

end