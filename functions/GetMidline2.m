function midlineInds = GetMidline2(IM,varargin)
% GetMidLine - Given an image series or an image dir containing an image
%   series, returns a series of indices corresponding to the midline of the
%   the fish in each image
% midlineInds = GetMidline(IM);
% midlineInds = GetMidline(..., fishPos);
% midlineInds = GetMidline(IM,fishPos, lineLens);
% Inputs:
% IM - Either an image stack (3D matrix where 3rd dim is time) or
%   image dir with image series. If IM is image dir, then assumes all
%   images are jpegs.
% fishPos - T x 2 matrix where T is the number of time points and the
%   columns are x and y positions of the fish. If not specified, or if
%   empty ([]), then reads the fish position automatically with default
%   parameters.
% lineLens - A vector where each value indicates the length of the line
%   segment in the series of line segments to fit to the fish's body
% Outputs:
% midlineInds - Indices of the midline in the image series; L x T matrix
%   where T is the number of images in the series and L is number of line
%   indices. L = sum(lineLens).

lineLens =  [18 16 14 10 8 8];
if nargin == 1
    if isdir(IM)
        IM  = ReadImgSequence(IM);
    end
    fishPos = GetFishPos(ProcessImages(IM),50);
elseif nargin == 2
    fishPos = varargin{1};
elseif nargin ==3
    fishPos = varargin{1};
    lineLens = varargin{2};
end

if isempty(fishPos)
    fishPos = GetFishPos(ProcessImages(IM),50);
end

midlineInds = {};
for imgNum = 1:size(IM,3)
    img = IM(:,:,imgNum);
    img = max(img(:))-img;
    lineInds = {};
    for jj = 1:length(lineLens)
        if jj ==1
            lineInds{jj} = GetML(blah,startPt,[],[],lineLens(jj));
        else
            lineInds{jj} = GetML(blah,startPt,prevStartPt,[],lineLens(jj));
        end
    end
    si = lineInds{jj}(end);
    [r,c] = ind2sub(size(blah),si);
    x = c;
    y = r;
    prevStartPt = startPt;
    startPt = [x,y];
    img = blah;
    midlineInds{imgNum} = lineInds;
    cla
    for kk = 1:length(lineInds)
        img(lineInds{kk}) = max(img(:));
        imagesc(img),axis image
    end
    hold on
    plot(fishPos(imgNum,1),fishPos(imgNum,2),'k*')
    title(num2str(imgNum))
    shg
    % pause(0.2)
end
end

function lineInds = GetML(im,startPt,varargin)
% lineInds = GetMidline(im,startInd,prevStartInd,dTh,lineLen)

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
    preVec = startPt-prevStartPt;
    preVec = preVec(1) + preVec(2)*1i;
    
    postVecs = [x y] - repmat(startPt,length(x),1);
    postVecs = postVecs(:,1) + postVecs(:,2)*1i;
    dAngles = angle(repmat(preVec,length(x),1).*conj(postVecs)) *(180/pi);
    farInds = find(abs(dAngles)<=90);
end


%## Find lines that pass through fish
backgroundInt= sort(rImg(:),'ascend');
backgroundInt(backgroundInt==0)=[];
backgroundInt = mean(backgroundInt(1:round(numel(rImg)/2)));
signalMat = rImg > 1*backgroundInt;
nPxls = sum(signalMat,2);
muPxls = mean(rImg,2);
muPxls = (muPxls-min(muPxls))/(max(muPxls)-min(muPxls));
slopes = rImg(:,end)-rImg(:,1);
nms = (nPxls(:).^1.5).*muPxls(:).*(1./abs(slopes(:)));
nms = nms(farInds);
nms(find(isinf(nms)))=1;
nms = (nms-min(nms(:)))/(max(nms)-min(nms));
nmsInds = find(nms>0.4);
blahInds = farInds(nmsInds);
nms = nms(nmsInds);

%## Find lines that are not contiguous blocks (i.e. islands) and eliminate
[blockSizes,blockInds] = GetContiguousBlocks(blahInds);
blockEndInds = blockInds(2:end)-1;
blockEndInds = [blockEndInds(:); length(blahInds)];
blockMaxes = zeros(size(blockInds));
for jj = 1:length(blockInds)
    blockMaxes(jj) = max(nms(blockInds(jj):blockEndInds(jj)));
end

%## For 1st line segment choose smaller block because head block will be
%## bigger than tail block

% if numel(farInds)== size(rImg,1);
%     [~, bigBlock] = min(blockSizes);
% else
% %     [~, bigBlock] = max(blockSizes);
%     [~,bigBlock] = max(blockSizes.*blockMaxes);
% end
[~,bigBlock]= max(blockSizes.*blockMaxes);

if numel(blockInds)> bigBlock
    keepInds = blockInds(bigBlock):blockInds(bigBlock+1);
else
    keepInds = blockInds(bigBlock):length(blahInds);
end

blahInds = blahInds(keepInds);
nms = nms(keepInds);
zerInds = find(blahInds==0);
blahInds(zerInds)=[];
nms(zerInds) = [];
throughFishInds = find(nPxls > 0.9*lineLen);

%## Find the line where the sum of pxls is largest
% finalInd = round(median(blahInds)); % This used to be a reasonable approx
% [~,maxInds] = max(sum(rImg(blahInds,:),2)); % Buggy too
[~,maxInds] = max(nms);

finalInd = blahInds(round(median(maxInds)));
lineInds = indMat(finalInd,:);
lineInds = lineInds(:);


    function [blockSizes, blockInds] = GetContiguousBlocks(values)
        % Given a set of values, returns the sizes of contiguous blocks and the
        %   block start indices
        blockInds = 1;
        blockSizes = [];
        count = 1;
        for jj = 1:length(values)-1
            if (values(jj+1)-values(jj))==1
                count = count+1;
            else
                blockSizes = [blockSizes; count];
                count = 0;
                blockInds =[blockInds; jj+1];
            end
        end
        blockSizes = [blockSizes; count];
    end

end

end
