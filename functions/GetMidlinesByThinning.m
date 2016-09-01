
function varargout = GetMidlinesByThinning(imgStack,varargin)
%GetMidlinesByThinnng Gets midline of fish in each image of imgStack by thinning
%   process
% mlInds = GetMidlinesByThinning(img);
% mlInds = GetMidlinesByThinning(img,'zThr',zThr,'mu',mu,'sigma',sigma,...
%   'minPxls',minPxls,'maxPxls',maxPxls,'fishPos',fishPos,'nIter',nIter);
% [mlInds,dsVecs,failedInds] = GetMidlinesByThinning(...);
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
% Outputs:
% mlInds - Cell array of the size of image stack where each cell contains
%   the midline pixel indices of the fish
% dsVecs - Cell array like mlInds, where each cell contains distances for
%   between midline pixels and fish head centroid pixel
% failedInds - Indices of images in which midline detection failed.
%
% Avinash Pujala, Koyama lab/HHMI, 2016

nThr = 4;
minThr = 0;
mu = [];
sigma = [];
minPxls = 20;
maxPxls = 100;
fishPos = [];
nIter = 20;
process = 'serial';
plotBool = 0;
kerSize = 7;
pauseDur = 0;
nWorkers = 10;

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
ker = ker/sum(ker(:));

mlInds = cell(size(imgStack,3),1);
dsVecs = mlInds;
imgInds = 1:size(imgStack,3);
dispChunk = round(size(imgStack,3)/5);
% dispChunk = 1;
if plotBool && ~strcmpi(process,'parallel')
    figure('Name','Midline inds by thinning')
end
imgDims = size(imgStack);
if strcmpi(process,'serial')
    count = 0;
    for tt = imgInds(:)'
        img = conv2(imgStack(:,:,tt),ker,'same');
        fp = fishPos(tt,:);         
        blah = GetMidlineByThinning(img,'nThr',nThr,'minThr',minThr,'mu',mu,'sigma',sigma,...
            'minPxls',minPxls,'maxPxls',maxPxls,'fishPos',fp);     
        [mlInds{tt},dsVecs{tt}] = GetMidlineCaudalToFishPos(fp,blah,imgDims(1:2));       
        if plotBool
            cla
            img(mlInds{tt})= mu + 20*sigma;
            imagesc(img),axis image
            hold on
            plot(fp(1),fp(2),'r+','markersize',10)
            title(['Img #' num2str(tt)])
            drawnow
            shg
            if isempty(pauseDur)
                pause()
            else
                pause(pauseDur)
            end
        end
        count= count + 1;
        if mod(tt,dispChunk)==0
            disp([num2str(100*(count/numel(imgInds))) '% completed...'])
        end
    end
else
    if matlabpool('size')==0
        matlabpool(nWorkers);
    end
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

lens = zeros(length(mlInds),1);
for jj = 1:length(mlInds)
    if iscell(mlInds)
        lens(jj) = length(mlInds{jj});
    else
        lens(jj) = length(mlInds);
    end
    
end
zerInds  = find(lens==0);
nonZerInds = setdiff(1:length(mlInds),zerInds);
if length(mlInds)>1
    for jj = zerInds(:)'
        if jj==1
            mlInds{jj} = mlInds{nonZerInds(1)};
            dsVecs{jj} = dsVecs{nonZerInds(1)};
        elseif jj == length(mlInds)
            mlInds{jj} = mlInds{jj-1};
            dsVecs{jj} = dsVecs{nonZerInds(end)};
        else
            [r0,c0] = ind2sub(imgDims(1:2),mlInds{jj-1});
            [r2,c2] = ind2sub(imgDims(1:2),mlInds{jj+1});
            ds0 = dsVecs{jj-1};
            ds2 = dsVecs{jj+1};
            n = min([length(r0),length(r2)]);
            rc = round(0.5*([r0(1:n), c0(1:n)] + [r2(1:n), c2(1:n)]));
            dsVecs{jj} = 0.5*(ds0(1:n)+ ds2(1:n));
            mlInds{jj} = sub2ind(imgDims(1:2),rc(:,1),rc(:,2));
        end       
    end
end

% disp('Correcting midline point order...')  
% mlInds = CorrectOrder(mlInds,imgDims); %-->  Feel like it's better to
%   correct after smoothening


varargout{1} = mlInds;
varargout{2} = dsVecs;
varargout{3} = zerInds;

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

[mlInds,img_denoised] = GetMLBT(img,nThr,minThr,minPxls,maxPxls,fishPos);
% c = 0;
% while isempty(mlInds) && c < 10
%     [mlInds,img_denoised] = GetMLBT(img_denoised,nThr,minThr,minPxls,maxPxls,fishPos);
% end

    function varargout = GetMLBT(img,nThr,minThr,minPxls,maxPxls,fishPos)
        S = @(v1,v2)(sqrt(sum((v1-v2).^2,2)));
        S2 = @(v1,v2)sqrt(sum((repmat(v1,size(v2,1),1)-v2).^2,2));
        img_bw = zeros(size(img));
        [thr,img_quant] = GetMultiThr(img,nThr,'minThr',minThr);
        thr(thr<=minThr)=[];
        thr(thr==0)=[];
        img_denoised = img;
        oneInds = [];
        lvls = unique(img_quant(:));
        lvls(lvls==0)=[];       
        count = numel(lvls);
        for lvl = lvls(:)'
            inds = find(img_quant==lvl);
            if numel(inds)< 1000
                oneInds = [oneInds; inds(:)];
            else
                img_denoised(img < thr(max(count,1)))=0;
            end
            count = count -1;
        end
        img_bw(oneInds)=1;
%         img_bw(img_quant>0)=1;
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
                if (s > 2*minPxls) || (s2 > 10)
                    img_thin(pxls)=0;
                end
            end
        end
        mlInds = find(img_thin);     
        varargout{1} = mlInds;
        varargout{2} = img_denoised;
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
if ~isempty(mlCoords_ord)
    mlCoords_ord = [round(fp); mlCoords_ord];
    s_ord = [0; s_ord];
end

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

function mlInds = CorrectOrder(mlInds,imgDims)
for n = 3:length(mlInds)
    clear sub1 sub2
    [sub1(:,2),sub1(:,1)] = ind2sub(imgDims(1:2), mlInds{n-1});
    [sub2(:,2),sub2(:,1)] = ind2sub(imgDims(1:2),mlInds{n});
    len = min([size(sub1,1), size(sub2,1)]);    
    a1 = GetDiffAngles(sub1);
    a2 = GetDiffAngles(sub2);   
    if len~=0
        d = sum(sqrt(sum((sub1(1:len,:)- sub2(1:len,:)).^2,2)),1);
        d_flip = sum(sqrt(sum((sub1(1:len,:)- sub2(len:-1:1,:)).^2,2)),1);cla
        if (d_flip < d) || ((a1*a2 <0) && (abs(a1-a2)> abs(a1--a2))) && (abs(a1-a2)>360)
            mlInds{n} = flipud(mlInds{n});
            disp(['Corrected frame # ' num2str(n)])       
        end
    end    
end

function A = GetDiffAngles(vec)
dV = diff(vec,[],1);
C = dV(:,1) + dV(:,2)*i;
A = angle(C(1:end-1).*conj(C(2:end)));
A = sum(A)*180/pi;
end
end

