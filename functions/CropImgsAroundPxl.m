
function IM_crop = CropImgsAroundPxl(IM,pxlPos,nPxls,varargin)
% CropImgsAroundPxl - Returns an image where pixels outside a radius from fish
%   head centroid are set to 0
% IM_crop = CropImgsAroundPxl(im)
% IM_crop = CropImgsAroundPxl(im,pxlPos)
% IM_crop = CropImgsAroundPxl(...,nPxls)
% IM_crop = CropImgsAroundPxl(...,'plotBool',plotBool,'pauseDur',pauseDur,'procType',procType,'poolSize',poolSize)
% Inputs:
% IM - Image stack to be cropped, where the 3rd dim is time, channel, etc
% pxlPos - x,y coordinates of pxl around which to crop images. Size of
%   pxlPos should be 2-by-N, where N = size(IM,3);
% nPxls - # of pxls on either side of specified pxl, thus dimensions of
%   IM_crop will be (2*nPxls+1)-by-(2*nPxls+1). Default is minimum of
%   half of smallest img dimension or 50.
% plotBool - Setting to true results in displaying of images. Default =
%   false
% pauseDur - Duration for which to pause between displaying images. Default
%   is zero.
% procType - Setting to 'parallel' results in cropping in parallel. Default
%   is 'serial'
% poolSize = Matlab poolsize for llel processing, default = 10
% 
% Avinash Pujala, Koyama lab/HHMI, 2016

plotBool = 0;
pauseDur = 0;
procType = 'serial';
poolSize = 10;

if nargin < 3
    error('Minimum 3 inputs required!')
end
for jj = 1:numel(varargin)
    if ischar(varargin{jj})
        switch lower(varargin{jj})
            case 'plotbool'
                plotBool= varargin{jj+1};
            case 'pausedur'
                pauseDur = varargin{jj+1};
            case 'proctype'
                procType = varargin{jj+1};
            case 'poolsize'
                poolSize = varargin{jj_1};
        end
    end
end
minDim = round(min([size(IM,1), size(IM,2)])/2);
if nPxls > minDim
    nPxls = minDim;
end
rowInds = 1:size(IM,1);
colInds = 1:size(IM,2);
IM_crop = zeros(nPxls*2+1, nPxls*2+1, size(IM,3));
dispChunk = ceil(size(IM,3)/5);
if plotBool
    figure('Name','Cropped images')
    cLim = [min(IM(:)),max(IM(:))];
end
imgInds = 1:size(IM,3);
if strcmpi(procType,'parallel')
    if matlabpool('size')==0
        matlabpool(poolSize)
    end
    parfor iNum = imgInds
        x = circshift(rowInds,[0,-pxlPos(iNum,2)]);
        x = [x(end-nPxls:end), x(1:nPxls)];
        y = circshift(colInds,[0,-pxlPos(iNum,1)]);
        y = [y(end-nPxls:end), y(1:nPxls)];
        IM_crop(:,:,iNum) = IM(x,y,iNum);
        if mod(iNum,dispChunk)==0
            disp(['Img #' num2str(iNum)])
        end
    end
else
    for iNum = imgInds
        x = circshift(rowInds,[0,-pxlPos(iNum,2)]);
        x = [x(end-nPxls:end), x(1:nPxls)];
        y = circshift(colInds,[0,-pxlPos(iNum,1)]);
        y = [y(end-nPxls:end), y(1:nPxls)];
        IM_crop(:,:,iNum) = IM(x,y,iNum);
        if plotBool
            imagesc(IM_crop(:,:,iNum)),axis image
            set(gca,'clim',cLim)
            hold on
            title(['Img #' num2str(iNum)])
            pause(pauseDur)
        end
        if mod(iNum,dispChunk)==0
            disp(['Img #' num2str(iNum)])
        end
    end
end
for iNum = 1:size(IM,3)
    x = circshift(rowInds,[0,-pxlPos(iNum,2)]);
    x = [x(end-nPxls:end), x(1:nPxls)];
    y = circshift(colInds,[0,-pxlPos(iNum,1)]);
    y = [y(end-nPxls:end), y(1:nPxls)];
    IM_crop(:,:,iNum) = IM(x,y,iNum);
    if plotBool
        imagesc(IM_crop(:,:,iNum)),axis image
        set(gca,'clim',cLim)
        hold on
        title(['Img #' num2str(iNum)])
        pause(pauseDur)
    end
    if mod(iNum,dispChunk)==0
        disp(['Img #' num2str(iNum)])
    end
end

end