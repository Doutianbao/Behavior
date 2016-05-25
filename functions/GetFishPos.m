function fishPos = GetFishPos(IM,varargin)
%GetFishPos Give an image stack returns the x,y coordinates of the fish in
%   each of the images of the stack
% fishPos = GetFishPos(IM);
% fishpos = GetFishPos(IM,nPixels)
% fishPos = GetFishPos(IM, nPixels,...)
% Inputs:
% IM - Image stack
% nPixels - # of bright pixels for centroid determination
% Optional input arguments
% 'method' - Centroid detection method: {'median'}, ['mean']
% 'filter' - Bandpass values, e.g., [15 25];
% 'process' - 'serial' or 'parallel'; Process in serial or parallel
% Avinash Pujala, HHMI, 2016

nPixels = 30;
method = 'mean';
process = 'serial';
filterFlag = 0;
poolSize = 12;

nArgs = length(varargin);
for jj = 1:nArgs
    if strcmpi(varargin{jj},'method')
        method =  varargin{jj+1};
    end
    if strcmpi(varargin{jj},'filter')
        bp =  varargin{jj+1};
        filterFlag = 1;
    end
    if strcmpi(varargin{jj},'process')
        process =  varargin{jj+1};
    end
end
if nargin > 1
    nPixels = varargin{1};
end


x = zeros(1,size(IM,3));
y = x;
if filterFlag
    [~,flt] = gaussianbpf(IM(:,:,1),bp(1),bp(2));
end
dispChunk = round(size(IM,3)/50)+1;
if strcmpi(process,'serial')
    disp('Tracking fish...')
    tic
    for jj=1:size(IM,3)
        img= IM(:,:,jj);
        if filterFlag
            img = gaussianbpf(img,flt);
        end
        [r,c] = FishPosInImg(img,nPixels,method);
        x(jj) = c;
        y(jj) = r;
        if mod(jj,dispChunk)==0
            disp(['Img # ' num2str(jj)])
            ShowFishPos(img,[r,c],jj)
        end
    end
    fishPos = [x; y]';
    toc
elseif strcmpi(process, 'parallel')
    imgFrames = 1:size(IM,3);
    tic
    disp('Tracking fish...')
    if matlabpool('size')==0
        matlabpool(poolSize)
    end
    if filterFlag
        parfor jj=imgFrames
            img= gaussianbpf(IM(:,:,jj),flt);
            [r,c] = FishPosInImg(img,nPixels,method);
            x(jj) = c;
            y(jj) = r;
            if mod(jj,dispChunk)==0
                disp(['Img # ' num2str(jj)])
            end
        end
    else
        parfor jj=imgFrames
            img= IM(:,:,jj);
            [r,c] = FishPosInImg(img,nPixels,method);
            x(jj) = c;
            y(jj) = r;  
            disp(['Img # ' num2str(jj)])
        end
    end    
    fishPos = [x; y]';
end

end

function [r,c] = FishPosInImg(img,nPixels,method)
[~,maxInds] = sort(img(:),'descend');
maxInds = maxInds(1:nPixels);
[r,c] = ind2sub(size(img),maxInds);
if strcmpi(method,'median')
    r = round(median(r));
    c = round(median(c));
elseif strcmpi(method,'mean')
    %             r = round(mean(r));
    %             c = round(mean(c));
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
