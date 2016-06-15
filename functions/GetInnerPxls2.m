
function varargout = GetInnerPxls2(im,periPxls)
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
if size(periPxls,2)==1
    periPxls = ind2sub(size(im),periPxls);
end

ctr  = round(mean(periPxls,1));
periPxls = fliplr(periPxls-repmat(ctr,size(periPxls,1),1)); % Converting from row-col to x-y inds
[th,rho] = cart2pol(periPxls(:,1),periPxls(:,2));
inPts =[];
N  = length(th);
dispChunk = round(N/10);
for thInd = 1:N
    if mod(thInd,dispChunk)==0
        disp([num2str(round(100*thInd/N)) '%'])
    end
    r = 1:rho(thInd);
    t = repmat(th(thInd),1,length(r));
    inPts = [inPts; t(:), r(:)];
end
[inCoords(:,1),inCoords(:,2)] = pol2cart(inPts(:,1),inPts(:,2));
inCoords = round(inCoords);
inPxls = fliplr(inCoords) + repmat(ctr,size(inCoords,1),1);
D = diff(inPxls,[],1);
f = find(D(:,1)==0 & D(:,2)==0);
inPxls(f,:)=[];
[outPxls(:,1),outPxls(:,2)] = ind2sub(size(im),setdiff(1:size(im,1)*size(im,2),sub2ind(size(im),inPxls(:,1),inPxls(:,2))));

varargout{1} = inPxls;
varargout{2} = outPxls;

end