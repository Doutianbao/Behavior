function fishPos = GetFishPos_parallel(IM,varargin)
%GetFishPos Give an image stack returns the x,y coordinates of the fish in
%   each of the images of the stack
% fishPos = GetFishPos(IM);
% fishpos = GetFishPos(IM,nPixels,method)
% Inputs:
% IM - Image stack
% nPixels - # of bright pixels for centroid determination (default: 40)
% method - Centroid detection method: {'median'}, ['mean']
% 
% Avinash Pujala, HHMI, 2016

nPixels = 40;
method = 'median';
poolSize = 10;
if nargin == 2
    nPixels = varargin{1};
elseif nargin == 3
    nPixels = varargin{1};
    method = varargin{2};
end


x = zeros(1,size(IM,3));
y = x;
imgFrames = 1:size(IM,3);
dispChunk= round(numel(imgFrames)/30);
imgHeight = size(IM,1);
imgWidth = size(IM,2);
tic
disp('Tracking fish...')
openedPool = 0;
if matlabpool('size')==0
    matlabpool(poolSize)
    openedPool = 1;
end
parfor jj=imgFrames
    ii= IM(:,:,jj);  
    [~,maxInds] = sort(ii(:),'descend');
    maxInds = maxInds(1:nPixels);
    [r,c] = ind2sub([imgHeight, imgWidth],maxInds);
    if strcmpi(method,'median')
        r = round(median(r));
        c = round(median(c));
    elseif strcmpi(method,'mean')
        r = round(mean(r));
        c = round(mean(c));
    end    
    x(jj) = r;
    y(jj) = c;
    if mod(jj,dispChunk)==0
        disp(['Img # ' num2str(jj)])
    end
end
fishPos = [y; x]';
if openedPool
    matlabpool close
end
end

