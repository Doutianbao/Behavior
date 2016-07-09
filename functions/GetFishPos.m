function varargout = GetFishPos(IM,nPxls,varargin)
%GetFishPos Give an image stack returns the x,y coordinates of the fish in
%   each of the images of the stack
% fishPos = GetFishPos(IM);
% fishpos = GetFishPos(IM,nPixels)
% fishPos = GetFishPos(IM, nPixels,...)
% [fishPos,headOrientationVector] = GetFishPos(IM,nPixels,...);
% [fishPos,headOr]
% Inputs:
% IM - Image stack
% nPixels - # of bright pixels for centroid determination. Setting nPixels
%   to [], results in default value of 30.
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

method = 'mean';
process = 'serial';
filterFlag = 0;
fltOrKer = [];
orFlag = 0;
nHood = 3;
plotBool = 1;
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
            case 'plotbool'
                plotBool = varargin{jj+1};
        end
    end
end

if isempty(nPxls)
    nPxls = 30;
end

if filterFlag
    if numel(bp)==1
        ker = gausswin(bp);
        ker = ker(:)*ker(:)';
        ker = ker/numel(ker);
        fltOrKer = ker;
        filterFlag = 2;
    elseif numel(bp)==2
        [~,flt] = gaussianbpf(IM(:,:,1),bp(1),bp(2));
        fltOrKer = flt;
    else
        error('Check inputs, filter input must be a scalar or a 2-element vector!')
    end
end

dispChunk = round(size(IM,3)/50)+1;
% dispChunk = 1;
fishPos = zeros(size(IM,3),2);
tic
disp('Tracking fish...')
if strcmpi(process,'serial')
    for jj=1:size(IM,3)
        [fishPos(jj,:),hOr{jj},img] = FishPosAndHeadVec(IM(:,:,jj),filterFlag,orFlag,...
            fltOrKer,nPxls,method,lineLen);
        if mod(jj,dispChunk)==0
            disp(['Img # ' num2str(jj)])
        end
        if plotBool
            ShowFishPos(img,[fishPos(jj,2),fishPos(jj,1)],hOr{jj},jj)
        end
    end
elseif strcmpi(process, 'parallel')
    imgFrames = 1:size(IM,3);
    if matlabpool('size')==0
        matlabpool(poolSize)
    end
    parfor jj=imgFrames
        [fishPos(jj,:),hOr{jj},img] = FishPosAndHeadVec(IM(:,:,jj),filterFlag,orFlag,...
            fltOrKer,nPxls,method,lineLen);
        if mod(jj,dispChunk)==0
            disp(['Img # ' num2str(jj)])
        end
    end
end
disp('Correcting head orientation for jumps...')
hOr_corr = JumpCorrectHeadOr(hOr);
toc
varargout{1} = fishPos;
varargout{2} = hOr_corr;
end

function varargout = FishPosAndHeadVec(img,filterFlag,orFlag,fltOrKer,nPxls,method,lineLen)
if filterFlag ==1
    img = gaussianbpf(img,fltOrKer);
elseif filterFlag ==2
    img = conv2(img,fltOrKer,'same');
end
[r,c] = FishPosInImg(img,nPxls,method);
x = c;
y = r;
if orFlag
    blah= GetMidlines(img,[c,r],lineLen,'plotBool',0);
    %             hOr{jj} = SmoothenMidline(blah{1}{1},img,nHood);
    inds = blah{1}{1};
    [yy,xx] = ind2sub(size(img),inds);
    hOr = [xx,yy];
    
else
    hOr = [];
end
varargout{1} = [x,y];
varargout{2} = hOr;
varargout{3} = img;
end

function [r,c] = FishPosInImg(img,nPxls,method)
[~,maxInds] = sort(img(:),'descend');
maxInds = maxInds(1:nPxls);
[r,c] = ind2sub(size(img),maxInds);
if strcmpi(method,'median')
    r = round(median(r));
    c = round(median(c));
elseif strcmpi(method,'mean')
    r = round(sum(r(:).*img(maxInds))/sum(img(maxInds)));
    c = round(sum(c(:).*img(maxInds))/sum(img(maxInds)));
end

end

function ShowFishPos(img,fishPos,hOr,imgNum)
cla
if size(hOr,2)==1
    [y,x] = ind2sub(size(img),hOr);
else
    x = hOr(:,1);
    y = hOr(:,2);
end

imagesc(img), axis image, colormap(gray)
hold on
plot(fishPos(2), fishPos(1),'ro')
plot(x,y,'g.')
title(['Frame # ' num2str(imgNum)])
drawnow
shg
end

function hOr_corr = JumpCorrectHeadOr(hOr)
A = @(v1,v2)angle(v1*conj(v2))*180/pi;
hOr_corr = hOr;
for n = 2:length(hOr_corr)-1   
    c_prev  = Line2Vec(hOr_corr{n-1});
    c_now = Line2Vec(hOr_corr{n});
    c_next = Line2Vec(hOr_corr{n+1});
    a1 = A(c_prev,c_now);
    a2 = A(c_now,c_next);
    a3 = A(c_prev,c_next);
    if (abs(a1)> abs(a3)) && (abs(a2)> abs(a3)) && (mean([abs(a1), abs(a2)])>100)
        hOr_corr{n}  = round((hOr_corr{n-1} + hOr_corr{n+1})/2);
        disp(['Corrected # ' num2str(n)])
    end
end
    function c = Line2Vec(xy)
        a = xy(1,:)-xy(end,:);
        c = a(1) + a(2)*1i;
    end    
end
