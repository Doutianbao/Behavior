
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
% 
% Avinash Pujala, Koyama lab/HHMI, 2016

if nargin < 2
    error('At least 2 inputs required!')
end
if size(periPxls,2)==1
    periPxls = ind2sub(size(im),periPxls);
end
%# Find center of periPxls
ctr  = round(mean(periPxls,1));
%# Translate image so as to bring the center of periPxls to the origin
[r,c] = find(im);
r = r- ctr(1);
c = c - ctr(2);

%# Convert shifted tranlated image coordinates to Cartesian coordinates
x = c; y = r;
periPxls = fliplr(periPxls-repmat(ctr,size(periPxls,1),1)); % Converting from row-col to x-y inds
[thP,rhoP] = cart2pol(periPxls(:,1),periPxls(:,2));
[thI, rhoI] = cart2pol(x,y);
epsilon = abs(diff(thP));
epsilon(epsilon==0)=[];
epsilon = min(epsilon)*2;
inPts =[];
N  = length(thP);
dispChunk = round(N/10);
disp('Finding pixels within input region...')
for thInd = 1:N
    if mod(thInd,dispChunk)==0
        disp([num2str(round(100*thInd/N)) '%'])
    end
    th = thP(thInd);
    rho = rhoP(thInd);
    blah = find(abs(thI-th)< (2*epsilon));
    blah = blah(rhoI(blah) < rho);
    inPts = [inPts; [thI(blah) rhoI(blah)]];
    thI(blah) = [];
    rhoI(blah)=[];    
end
[inCoords(:,1),inCoords(:,2)] = pol2cart(inPts(:,1),inPts(:,2));
inCoords = round(inCoords);
inPxls = fliplr(inCoords) + repmat(ctr,size(inCoords,1),1);
% D = diff(inPxls,[],1);
% f = find(D(:,1)==0 & D(:,2)==0);
% inPxls(f,:)=[];
[outPxls(:,1),outPxls(:,2)] = ind2sub(size(im),setdiff(1:size(im,1)*size(im,2),sub2ind(size(im),inPxls(:,1),inPxls(:,2))));

varargout{1} = inPxls;
varargout{2} = outPxls;

end