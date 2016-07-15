function mlInds = GetMidlineByThinning(img,varargin)
%GetMidlineByThinnng Gets midline of fish in input image by thinning
%   process
% mlInds = GetMidlinesByThinning(img);
% mlInds =
% GetMidlinesByThinning(img,'zThr',zThr,'mu',mu,'sigma',sigma,'minPxls',minPxls,'maxPxls',maxPxls,'nIter',nIter);
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

zThr = 2;
mu = [];
sigma = [];
minPxls = 40;
maxPxls = 100;
fishPos = [];
nIter = 20;

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

if isempty(mu)
    Z = @(x) (x-mean(x(:)))/std(x(:));
    img = Z(img);
else
    img = (img-mu)/sigma;
end


thr = zThr;
mlInds = GetMLBT(img,thr);
count = 0;
while (numel(mlInds) < minPxls) && (count < nIter)
    count = count + 1;
    disp(['Lowering threshold..., iter # ' num2str(count)])
    mlInds = GetMLBT(img,thr);
    thr = 0.95*thr;
end

count = 0;
while (numel(mlInds) > maxPxls) && (count < nIter)
    count = count + 1;
    disp(['Raising threshold..., iter # ' num2str(count) ])
    mlInds = GetMLBT(img,thr);
    thr = 1.05*thr;
end

    function mlInds = GetMLBT(img,zThr)
        img_bw = img;
        img_bw(img_bw<=zThr)=0; img_bw(img_bw>zThr)=1;
        img_thin = bwmorph(img_bw,'thin',Inf);
        cc = bwconncomp(img_thin);
        rp = regionprops(cc,'Centroid');
        for jj = 1:length(cc.PixelIdxList)
            pxls = cc.PixelIdxList{jj};
            if numel(pxls)<minPxls
                img_thin(pxls)=0;
            elseif numel(pxls) > maxPxls
                img_thin(pxls)=0;
            end
            if ~isempty(fishPos)
                S = @(v1,v2)(sqrt(sum((v1-v2).^2,1)));
                s = S(rp(jj).Centroid,fishPos);
                if s > minPxls
                    img_thin(pxls)=0;
                end
            end
        end
        
        mlInds = find(img_thin);
    end

end

