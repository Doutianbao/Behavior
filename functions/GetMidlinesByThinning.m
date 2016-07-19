
function mlInds = GetMidlinesByThinning(imgStack,varargin)
%GetMidlinesByThinnng Gets midline of fish in each image of imgStack by thinning
%   process
% mlInds = GetMidlinesByThinning(img);
% mlInds = GetMidlinesByThinning(img,'zThr',zThr,'mu',mu,'sigma',sigma,...
%   'minPxls',minPxls,'maxPxls',maxPxls,'fishPos',fishPos,'nIter',nIter);
% Inputs:
% img - Image with fish in it
% 'zThr' - Threshold for image binarization after expressing pixel
%   intensities as z-score
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
% 'process' - Setting to 'parallel' results in parallel processing
% 'plotBool' - Seting to 1, results in diplaying of images with midlines.

zThr = 2;
mu = [];
sigma = [];
minPxls = 40;
maxPxls = 100;
fishPos = [];
nIter = 20;
process = 'serial';
plotBool = 0;

for in = 1:numel(varargin)
    if ischar(varargin{in})
        switch lower(varargin{in})
            case 'zthr'
                zThr = varargin{in+1};
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
        end
    end
end

if isempty(mu)
    Z = @(x) (x-mean(x(:)))/std(x(:));
    img = Z(imgStack);
else
    img = (imgStack-mu)/sigma;
end

mlInds = cell(size(imgStack,3),1);
imgInds = 1:size(imgStack,3);
dispChunk = round(size(imgStack,3)/5);
if plotBool
    figure('Name','Midline inds by thinning')
end
if strcmpi(process,'parallel')
    for tt = imgInds(:)'
        img = imgStack(:,:,tt);
        fp = fishPos(tt,:);
        mlInds{tt} = GetMidlineByThinning(img,'zThr',zThr,'mu',mu,'sigma',sigma,...
            'minPxls',minPxls,'maxPxls',maxPxls,'fishPos',fp,'nIter',nIter);
        if plotBool
            cla
            img(mlInds{tt})= mu + 20*sigma;
            imagesc(img),axis image
            title(['Img #' num2str(tt)])
        end
        if mod(tt,dispChunk)==0
            disp(['Img # ' num2str(tt)])
        end
    end
else
    parfor tt = imgInds(:)'
        img = imgStack(:,:,tt);
        fp = fishPos(tt,:);
        mlInds{tt} = GetMidlineByThinning(img,'zThr',zThr,'mu',mu,'sigma',sigma,...
            'minPxls',minPxls,'maxPxls',maxPxls,'fishPos',fp,'nIter',nIter);
        if mod(tt,dispChunk)==0
            disp(['Img # ' num2str(tt)])
        end
    end
end


if length(mlInds)==1
    mlInds = mlInds{1};
end

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


for in = 1:numel(varargin)
    if ischar(varargin{in})
        switch lower(varargin{in})
            case 'zthr'
                zThr = varargin{in+1};
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
        end
    end
end


thr = zThr;
mlInds = GetMLBT(img,thr,minPxls,maxPxls,fishPos);
count = 0;
while (numel(mlInds) < minPxls) && (count < nIter)
    count = count + 1;
    disp(['Lowering threshold..., iter # ' num2str(count)])
    mlInds = GetMLBT(img,thr,minPxls,maxPxls,fishPos);
    thr = 0.95*thr;
end

count = 0;
while (numel(mlInds) > maxPxls) && (count < nIter)
    count = count + 1;
    disp(['Raising threshold..., iter # ' num2str(count) ])
    mlInds = GetMLBT(img,thr);
    thr = 1.05*thr;
end


    function mlInds = GetMLBT(img,zThr,minPxls,maxPxls,fishPos)
        S = @(v1,v2)(sqrt(sum((v1-v2).^2,2)));
        S2 = @(v1,v2)sqrt(sum((repmat(v1,size(v2,1),1)-v2).^2,2));
        img_bw = img;
        img_bw(img_bw<=zThr)=0; img_bw(img_bw>zThr)=1;
        img_thin = bwmorph(img_bw,'thin',Inf);
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


