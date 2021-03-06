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


iThr_def = 0.04*im(fishPos(2),fishPos(1));
iThr = iThr_def;
% iThr_def = multithresh(im,3);
% iThr = iThr_def(2);

dThr = 60;
if nargin ==3
    iThr = varargin{1};    
elseif nargin == 4
    iThr = varargin{1};
    dThr = varargin{2};
end
if isempty(iThr)
    iThr = iThr_def;
end

%#### Find all pixels >= iThr within the radius dThr from fish
%      position


[row,col] = find(ones(size(im)));
S = sqrt(sum([col-fishPos(1) row-fishPos(2)].^2,2));

outPxls = find(S>=dThr);
inPxls = find(S<dThr);
inBrightPxls = inPxls(im(inPxls)>=iThr);
im_fish = zeros(size(im));
im_fish(inBrightPxls)=1;
inDimPxls = setdiff(inPxls,inBrightPxls);
blah = im;
blah(inDimPxls)=-20;

%## Find blobs and get some properties
cc = bwconncomp(im_fish);
stats = regionprops(cc,'Area','Extent','Extrema');

%## Identify head blob
count = 0;
for jj = 1:length(cc.PixelIdxList)
    [row,col] = ind2sub(size(im),cc.PixelIdxList{jj});
    S = sqrt(sum([col-fishPos(1) row-fishPos(2)].^2,2));
    if sum(S==0)
        count = count+1;
        if count > 1
            error('Head blob detected twice!')
        end
        headBlobIdx = jj;
    end
end

%## Compare blobs to head blob and mark as not fish based on some
%   properties
blobInds = setdiff(1:length(cc.PixelIdxList),headBlobIdx);
nonFishBlobInds = [];
nonFishPxlInds = [];
for jj = blobInds(:)'
    if  (stats(jj).Area < (stats(headBlobIdx).Area/20)) || (stats(jj).Extent < (stats(headBlobIdx).Extent/10))
        nonFishBlobInds = [nonFishBlobInds; jj];
        pxlInds = cc.PixelIdxList{jj};
        nonFishPxlInds = [nonFishPxlInds; pxlInds(:)];
    end
end

%## Remove non-fish blobs
if ~isempty(nonFishPxlInds)
    im_fish(nonFishPxlInds) = 0;
end


%## Connecting the fish blobs
if ~isempty(blobInds)
ex1 = stats(headBlobIdx).Extrema;
distToHeadBlob = [];
extremaMat =[];
for jj = blobInds(:)'
    ex2 = stats(jj).Extrema;
    xMat = repmat(ex1(:,1),1,length(ex2(:,1)))-repmat(ex2(:,1)',length(ex1(:,1)),1);
    yMat =  repmat(ex1(:,2),1,length(ex2(:,2)))-repmat(ex2(:,2)',length(ex1(:,2)),1);
    sMat = sqrt(xMat.^2 + yMat.^2);
    [headMin,blobMin] = find(sMat==min(sMat(:)));
    if numel(headMin)>1
        headMin = headMin(1);
        blobMin = blobMin(1);
    end
    distToHeadBlob = [distToHeadBlob; sMat(headMin,blobMin)];
    extremaMat = [extremaMat;[headMin, blobMin]];
end

closestToHeadBlob = blobInds(find(distToHeadBlob==min(distToHeadBlob)));

ptsToConnect = [stats(headBlobIdx).Extrema(headMin,:); stats(closestToHeadBlob).Extrema(blobMin,:)];
rcMat = floor(fliplr(ptsToConnect));
% im_fish(rcMat(1,1):rcMat(2,1), rcMat(1,2):rcMat(2,2))=1;

lineCoords = LineBetween2Pts(rcMat(1,:),rcMat(2,:));
for jj = 1:size(lineCoords,1)
    im_fish(lineCoords(jj,1), lineCoords(jj,2))=1;
end
end

%## Close gap between fish blobs
% se = strel('disk',20);
% im_fish = imclose(im_fish,se);


%%
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
%     varargout{1} = B;
    varargout{1} = im_fish;
else
    varargout{2} = find(B);    
end

end

function pixelCoords = LineBetween2Pts(pt1,pt2)
if numel(pt1)~=2 & numel(pt2)~=2
    error('Input must be a point in 2D space');    
end
rows = pt1(1):pt2(1); 
cols = pt1(2):pt2(2);
nPts = max(numel(rows),numel(cols));
rows = round(linspace(pt1(1),pt2(1),nPts));
cols = round(linspace(pt1(2),pt2(2),nPts));
% if numel(rows) > numel(cols)
%     rows = round(linspace(pt1(1),pt2(1),numel(rows)));
% elseif numel(cols) > numel(rows)
%     cols = round(linspace(pt1(2),pt2(2),numel(cols)));
% end
pixelCoords = [rows(:) cols(:)];
end
