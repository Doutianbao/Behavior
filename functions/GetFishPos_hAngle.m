function varargout = GetFishPos(IM,varargin)
%GetFishPos Give an image stack returns the x,y coordinates of the fish in
%   each of the images of the stack
% fishPos = GetFishPos(IM);
% fishpos = GetFishPos(IM,nPixels)
% fishPos = GetFishPos(IM, nPixels,...)
% [fishPos,headOrientationVector] = GetFishPos(IM,nPixels,...);
% Inputs:
% IM - Image stack
% nPixels - # of bright pixels for centroid determination
% Optional input arguments
% 'method' - Centroid detection method: {'median'}, ['mean']
% 'filter' - Bandpass values, e.g., [15 25]; If bp is a scalar then creates
%   a gaussian kernel of this size and uses this to convolve the image
%   instead of using 'gaussianbpf'
% 'process' - 'serial' or 'parallel'; Process in serial or parallel
% 'lineLen' - Length of line segment to use to determine head orientation
%   (default: 15 pxls)
% 'pxlLim' - Limit of how many pixels the fish can traverse in the imaging
%   setup. This is can be used to reduce false positives. (Not yet implemented)
%
% Outputs:
% fishPos - x, y coordinates of fish head position
% headOrientationVector - x,y coordinates of the head orientation vector (length
%   is determined by the variable 'lineLen')
%
% Avinash Pujala, HHMI, 2016

nPixels = 30;
method = 'mean';
process = 'serial';
filterFlag = 0;
orFlag = 0;
nHood = 3;
poolSize = 10;

nArgs = length(varargin);
for jj = 1:nArgs
    if ischar(varargin{jj})
        switch lower(varargin{jj})
            case 'method'
                method =  varargin{jj+1};
            case 'filter'
                bp =  varargin{jj+1};
                if ~isempty(bp)
                    filterFlag = 1;
                end
            case 'process'
                process =  varargin{jj+1};
            case lower('pxlLim')
                pxlLim = varargin{jj+1};
            case lower('lineLen')
                lineLen = varargin{jj+1};
                orFlag = 1;
                hOr = cell(size(IM,3),1);
                if isempty(lineLen)
                    lineLen = 15;
                end
        end
    end
end

if nargin > 1
    nPixels = varargin{1};
end

x = zeros(1,size(IM,3));
y = x;

if filterFlag
    if numel(bp)==1
        ker = gausswin(bp);
        ker = ker(:)*ker(:)';
        ker = ker/numel(ker);
        filterFlag = 2;
    elseif numel(bp)==2
        [~,flt] = gaussianbpf(IM(:,:,1),bp(1),bp(2));
    else
        error('Check inputs, filter input must be a scalar or a 2-element vector!')
    end
end

dispChunk = round(size(IM,3)/50)+1;
% dispChunk = 1;
if strcmpi(process,'serial')
    disp('Tracking fish...')
    tic
    for jj=1:size(IM,3)
        img= IM(:,:,jj);
        if filterFlag ==1
            img = gaussianbpf(img,flt);
        elseif filterFlag ==2
            img = conv2(img,ker,'same');
        end
        [r,c] = FishPosInImg(img,nPixels,method);
        x(jj) = c;
        y(jj) = r;
        if orFlag
            blah= GetMidlines(img,[c,r],lineLen,'plotBool',1);
             hOr{jj} = SmoothenMidline(blah{1}{1},img,nHood);            
        end
        if mod(jj,dispChunk)==0
            disp(['Img # ' num2str(jj)])
            ShowFishPos(img,[r,c],jj)
        end
    end
    toc
elseif strcmpi(process, 'parallel')
    imgFrames = 1:size(IM,3);
    tic
    disp('Tracking fish...')
    if matlabpool('size')==0
        matlabpool(poolSize)
    end
    parfor jj=imgFrames
        if filterFlag ==1
            img= gaussianbpf(IM(:,:,jj),flt);
        elseif filterFlag ==2
            img = conv2(IM(:,:,jj),ker,'same');
        else
            img = IM(:,:,jj);
        end
        [r,c] = FishPosInImg(img,nPixels,method);
        x(jj) = c;
        y(jj) = r;
        if orFlag
            blah= GetMidlines(img,[c,r],lineLen);
            hOr{jj} = SmoothenMidline(blah,img,nHood);
        end
        if mod(jj,dispChunk)==0
            disp(['Img # ' num2str(jj)])
        end
    end
end
fishPos = [x; y]';
varargout{1} = fishPos;
varargout{2} = hOr;
end

function [r,c] = FishPosInImg(img,nPixels,method)
[~,maxInds] = sort(img(:),'descend');
maxInds = maxInds(1:nPixels);
[r,c] = ind2sub(size(img),maxInds);
if strcmpi(method,'median')
    r = round(median(r));
    c = round(median(c));
elseif strcmpi(method,'mean')
    r = round(sum(r(:).*img(maxInds))/sum(img(maxInds)));
    c = round(sum(c(:).*img(maxInds))/sum(img(maxInds)));
end

end
function ShowFishPos(img,fishPos,imgNum)
cla
imagesc(img), axis image, colormap(gray)
hold on
plot(fishPos(2), fishPos(1),'ro')
title(['Frame # ' num2str(imgNum)])
drawnow
shg
end

function tailCurv = SmoothenMidline(mlInds,img,nHood)
tailCurv = zeros(length(mlInds),2);
[r,c] = ind2sub(size(img),mlInds);
[C,R] = meshgrid(1:size(img,2),1:size(img,1));
for jj = 1:size(r,1)
    rInds = r(jj)-nHood:r(jj)+nHood;
    rInds(rInds<=0)=[];
    cInds = c(jj)-nHood:c(jj)+nHood;
    cInds(cInds<=0)=[];
    rNeighbors = R(rInds,cInds);
    cNeighbors = C(rInds,cInds);
    wts = img(rInds,cInds);
%     ker = ones(length(wts));
%     wts = conv2(wts,ker,'same');
    pxlInd(1) = sum(rNeighbors(:).*wts(:))/sum(wts(:));
    pxlInd(2) = sum(cNeighbors(:).*wts(:))/sum(wts(:));
    tailCurv(jj,:)= fliplr(round(pxlInd*10)/10); % Flipping to give in x-y rather than row-col coordinates
end
end

function tailCurv_spline = SplineTailCurv(tailCurv,smoothFactor)
% y = [tailCurv(1,:); tailCurv; tailCurv(end,:)];
y = tailCurv;
t = 1:size(y,1);
ts = linspace(1,size(y,1),size(y,1)/smoothFactor);
ys(:,1) = spline(t,y(:,1),ts);
ys(:,2) = spline(t,y(:,2),ts);
y(:,1) = spline(ts,ys(:,1),t);
y(:,2) = spline(ts,ys(:,2),t);
tailCurv_spline = y;
% tailCurv_spline = [xx(2:end-1); yy(2:end-1)]';

end
