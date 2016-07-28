
function varargout = GetMidlinesByThinning(imgStack,varargin)
%GetMidlinesByThinnng Gets midline of fish in each image of imgStack by thinning
%   process
% mlInds = GetMidlinesByThinning(img);
% mlInds = GetMidlinesByThinning(img,'zThr',zThr,'mu',mu,'sigma',sigma,...
%   'minPxls',minPxls,'maxPxls',maxPxls,'fishPos',fishPos,'nIter',nIter);
% [mlInds,dsVecs]
% Inputs:
% img - Image with fish in it
% 'nThr' - Number of threshold levels to go through to identify fish (see GetMultiThr)
% 'minThr' - Minimum threshold value to use in determining fish pixels for
%   binarization (see GetMultiThr).
% 'fishPos' - Position of fish, so as to exlude midlines that are too far.
% 'mu' - Mean to subtract when computing z-scores. If not specified, mu is
%   computed for input image, but when available it is better to use the mu
%   computed for the entire image stack.
% 'sigma' - Standard deviation to use for computign z-scores. If not
%   specified, it is computed.
% 'minPxls' - Min number of presumed pxls for midline. Default - 40
% 'maxPxls' - Max # of presumed pxls for midline. Default - 100
% 'nIter' - Number of iteration to vary binarization threshold to try and
%   find real midline indices.
% 'kerSize' - Gaussian kernel size for smoothing of image(s).
% 'process' - Setting to 'parallel' results in parallel processing
% 'plotBool' - Seting to 1, results in diplaying of images with midlines
% 'pauseDur' - Pause interval between successive plots
% Outputs:
% mlInds - A cell array of the size of the image stack, where each cell
%   contains the fish midline indices for the corresponding image in the
%   stack
% dsVecs - A cell array  like mlInds, but where each cell contains the
%   distances from the head centroid of the fish to the midline points.
% 
% Avinash Pujala, Koyama lab/HHMI, 2016

nThr = 4;
minThr = 0;
mu = [];
sigma = [];
minPxls = 30;
maxPxls = 100;
fishPos = [];
nIter = 20;
process = 'serial';
plotBool = 0;
kerSize = 7;
pauseDur = 0;

for in = 1:numel(varargin)
    if ischar(varargin{in})
        switch lower(varargin{in})      ;
            case 'nthr'
                nThr = varargin{in+1};
            case 'minthr'
                minThr = varargin{in+1};
            case 'fishpos'
                fishPos = varargin{in+1};
            case 'mu'
                mu = varargin{in+1};
            case 'sigma'
                sigma = varargin{in+1};
            case 'minpxls'
                minPxls = varargin{in+1};
            case 'maxpxls'
                maxPxls = varargin{in+1};
            case 'niter'
                nIter = varargin{in+1};
            case 'process'
                process = varargin{in+1};
            case 'plotbool'
                plotBool = varargin{in+1};
            case 'kersize'
                kerSize = varargin{in+1};
            case 'pausedur'
                pauseDur = varargin{in+1};
        end
    end
end

if isempty(mu)
    Z = @(x) (x-mean(x(:)))/std(x(:));
    imgStack = Z(imgStack);
    mu = 0;
    sigma = 1;
else
    imgStack = (imgStack-mu)/sigma;
end
ker = gausswin(kerSize)*gausswin(kerSize)';
ker = ker/sum(ker);

mlInds = cell(size(imgStack,3),1);
dsVecs = mlInds;
imgInds = 1:size(imgStack,3);
dispChunk = round(size(imgStack,3)/5);
if plotBool
    figure('Name','Midline inds by thinning')
end
imgDims = size(imgStack);
if strcmpi(process,'serial')
    for tt = imgInds(:)'
        img = conv2(imgStack(:,:,tt),ker,'same');
        fp = fishPos(tt,:);
        blah = GetMidlineByThinning(img,'nThr',nThr,'minThr',minThr,'mu',mu,'sigma',sigma,...
            'minPxls',minPxls,'maxPxls',maxPxls,'fishPos',fp);
        [blah,dsVecs{tt}] = GetMidlineCaudalToFishPos(fp,blah,imgDims(1:2));
        mlInds{tt} = blah;
        if plotBool
            cla
            img(mlInds{tt})= mu + 20*sigma;
            imagesc(img),axis image
            hold on
            plot(fp(1),fp(2),'r+','markersize',10)
            title(['Img #' num2str(tt)])
            shg
            if isempty(pauseDur)
                pause()
            else
                pause(pauseDur)
            end
        end
        if mod(tt,dispChunk)==0
            disp(['Img # ' num2str(tt)])
        end
    end
else
    parfor tt = imgInds(:)'
        img = conv2(imgStack(:,:,tt),ker,'same');
        fp = fishPos(tt,:);        
       blah = GetMidlineByThinning(img,'nThr',nThr,'minThr',minThr,'mu',mu,'sigma',sigma,...
            'minPxls',minPxls,'maxPxls',maxPxls,'fishPos',fp);
        [blah,dsVecs{tt}] = GetMidlineCaudalToFishPos(fp,blah,imgDims(1:2));
        mlInds{tt} = blah;
        if mod(tt,dispChunk)==0
            disp(['Img # ' num2str(tt)])
        end
    end
end


if length(mlInds)==1
    mlInds = mlInds{1};
end
if length(dsVecs)==1
    dsVecs = dsVecs{1};
end

varargout{1} = mlInds;
varargout{2} = dsVecs;

end

function mlInds = GetMidlineByThinning(img,varargin)
%GetMidlineByThinnng Gets midline of fish in input image by thinning
%   process
% mlInds = GetMidlinesByThinning(img);
% mlInds = GetMidlinesByThinning(img,'zThr',zThr,'mu',mu,'sigma',sigma,...
%   'minPxls',minPxls,'maxPxls',maxPxls,'fishPos',fishPos,'nIter',nIter);
% Inputs:
% img - Image with fish in it
% 'zThr' - Threshold for image binarization after expressing pixel
%   intensities as z-score
% 'nThr' - Number of threshold levels
% 'minThr' - Minimum value for threshold
% 'fishPos' - Position of fish, so as to exlude midlines that are too far.
% 'mu' - Mean to subtract when computing z-scores. If not specified, mu is
%   computed for input image, but when available it is better to use the mu
%   computed for the entire image stack.
% 'sigma' - Standard deviation to use for computign z-scores. If not
%   specified, it is computed.
% 'minPxls' - Min number of presumed pxls for midline. Default - 40
% 'maxPxls' - Max # of presumed pxls for midline. Default - 100


for in = 1:numel(varargin)
    if ischar(varargin{in})
        switch lower(varargin{in})
            case 'nthr'
                nThr = varargin{in+1};
            case 'minthr'
                minThr = varargin{in+1};
            case 'fishpos'
                fishPos = varargin{in+1};
            case 'mu'
                mu = varargin{in+1};
            case 'sigma'
                sigma = varargin{in+1};
            case 'minpxls'
                minPxls = varargin{in+1};
            case 'maxpxls'
                maxPxls = varargin{in+1};
        end
    end
end

mlInds = GetMLBT(img,nThr,minThr,minPxls,maxPxls,fishPos);

    function mlInds = GetMLBT(img,nThr,minThr,minPxls,maxPxls,fishPos)
        S = @(v1,v2)(sqrt(sum((v1-v2).^2,2)));
        S2 = @(v1,v2)sqrt(sum((repmat(v1,size(v2,1),1)-v2).^2,2));
        img_bw = zeros(size(img));
        [~,img_quant] = GetMultiThr(img,nThr,'minThr',minThr);
        img_bw(img_quant>0)=1;
        img_thin = bwmorph(bwmorph(img_bw,'thin',Inf),'spur');
        cc = bwconncomp(img_thin);
        rp = regionprops(cc,'Centroid');
        imgDims= size(img);
        for region = 1:length(cc.PixelIdxList)
            pxls = cc.PixelIdxList{region};
            if numel(pxls)<minPxls
                img_thin(pxls)=0;
            elseif numel(pxls) > maxPxls
                img_thin(pxls)=0;
            end
            v2 = [];
            if ~isempty(fishPos)
                s = S(rp(region).Centroid,fishPos);
                [v2(:,2),v2(:,1)] = ind2sub(imgDims,pxls);
                s2 = min(S2(fishPos,v2));
                if (s > minPxls) || (s2 > 10)
                    img_thin(pxls)=0;
                end
            end
        end
        mlInds = find(img_thin);
    end
end

function [mlInds, distVec] = GetMidlineCaudalToFishPos(fp,mlInds,imgDims)
[mlCoords(:,2),mlCoords(:,1)] = ind2sub(imgDims,mlInds);

%## Find midline pts rostral to head centroid pos and remove
[mlCoords_ord,s_ord] = IncreasinglyDistantPtsFromRef(fp,mlCoords);
%  jumpInd = find(s_ord > (mean(s_ord) + 1.5*std(s_ord)));
[~,jumpInd] = max(s_ord);
if jumpInd > length(s_ord)-jumpInd
    mlCoords_ord(jumpInd:end,:)=[];
    s_ord(jumpInd:end) = [];
else
    mlCoords_ord(1:jumpInd-1,:)=[];
    s_ord(1:jumpInd-1) = [];
end
%###

mlInds = sub2ind(imgDims,mlCoords_ord(:,2),mlCoords_ord(:,1));
distVec = cumsum(s_ord);

end

function [pts_order,s_order] = IncreasinglyDistantPtsFromRef(refPt,pts)
S = @(v,V) sqrt(sum((repmat(v,size(V,1),1)-V).^2,2));
pts_new = pts;
s = S(refPt,pts_new);
s_new = s;
count= 0;
pts_order = zeros(size(pts));
s_order = zeros(size(s));
while length(s_new)>1
    count = count + 1;
    [d, ind] = min(s_new);
    s_order(count)=d;
    pts_order(count,:) = pts_new(ind,:);
    refPt = pts_new(ind,:);
    %     s_new(ind)=[];
    pts_new(ind,:) = [];
    s_new = S(refPt,pts_new);
end
zerInds = find(pts_order(:,1)==0 | pts_order(:,2)==0);
pts_order(zerInds,:) = [];
s_order(zerInds) = [];
end



