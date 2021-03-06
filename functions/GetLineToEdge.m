
function lineInds = GetLineToEdge(im,startPt,varargin)
% lineInds = GetLineToEdge(im,startInd,prevStartInd,dTh,lineLen)
dTh = 5;
lineLen = 15;
if nargin < 2
    error('At least 2 inputs required!')
elseif nargin ==2
    prevStartPt = [];
elseif nargin == 3
    prevStartPt = varargin{1};
elseif nargin ==4
    prevStartPt = varargin{1};
    dTh = varargin{2};
elseif nargin == 5
    prevStartPt = varargin{1};
    dTh = varargin{2};
    lineLen = varargin{3};
end

if isempty(dTh)
    dTh  = 5;
end

[rImg,indMat] = RadialFish(im,startPt,dTh,lineLen);

if isempty(prevStartPt)
    farInds = 1:size(indMat,1);
else
    lastInds = indMat(:,end);
    [rL,cL] = ind2sub(size(im),lastInds);
    x = cL; y = rL;
    S1 = sqrt(sum((startPt-prevStartPt).^2,2));
    S2 = sqrt(sum(([x y] - repmat(prevStartPt,length(rL),1)).^2,2));
    S2_90 = sqrt(S1^2 + lineLen^2);
    farInds = find(S2>1*S2_90);
end

%## Find lines that that move farther way from previous starting point and
%   that have pixels that are of smaller value that previous pixels along
%   the line

dR = diff(rImg,[],2);
[~,minInds] = min(dR,[],2);
% r = 1:length(minInds);
% negInds = dR(sub2ind(size(dR),r(:),minInds(:)));
% minNegInds = find(negInds<0);
% minNegFarInds =  intersect(farInds,minNegInds);

negSlope = rImg(:,end)-rImg(:,end-1);
thr = mean(negSlope(negSlope<0)) - 0.5*std(negSlope(negSlope<0));
edgeDropInds = find(negSlope < thr);

%## Find lines that pass through fish
backgroundInt= sort(rImg(:),'ascend');
backgroundInt(backgroundInt==0)=[];
backgroundInt = mean(backgroundInt(1:round(numel(rImg)/2)));
signalMat = rImg > 1*backgroundInt;
nPxls = sum(signalMat,2);
muPxls = mean(rImg,2);
muPxls = (muPxls-min(muPxls))/(max(muPxls)-min(muPxls));
throughFishInds = find((nPxls >= round(0.9*lineLen)) & (muPxls >0.2));


blah = repmat(1:size(signalMat,2),size(signalMat,1),1);
blah(~signalMat) = 0;
maxInds = max(blah,[],2);
edgeOffInds = find(maxInds >= (lineLen-1));

% edgeOnInds = find(maxInds < 5);

% mu = mean(rImg(minNegFarInds,:),2);
% sig = std(rImg(minNegFarInds,:),[],2);
% cv = sig./mu;
% [~,throughFishInds] = sort(cv,'ascend');
% throughFishInds = throughFishInds(1:min(numel(cv),5));

finalInds = intersect(intersect(intersect(farInds, edgeDropInds),edgeOffInds),throughFishInds);
if isempty(finalInds)
    finalInds =  intersect(intersect(farInds, edgeOffInds),throughFishInds);
    if isempty(finalInds)
        finalInds = intersect(farInds,throughFishInds);
        if isempty(finalInds)
            finalInds = throughFishInds;
        end
    end
end
[~, finalInd] = min(rImg(finalInds,end));
finalInd = finalInds(finalInd);

% minNegFarThroughInds = intersect(minNegFarInds,throughFishInds);
% minNegFarThroughInds = minNegFarInds(throughFishInds);

if isempty(finalInds)
    a = 1;
end

% inds =  minInds(minNegFarInds);
% inds = find(inds>= (lineLen-4));
% % [~,inds] = max(inds);
% inds = minNegFarInds(inds(1));

lineInds = indMat(finalInd,:);
lineInds = lineInds(:);

end