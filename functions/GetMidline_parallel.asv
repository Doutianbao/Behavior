function midlineInds = GetMidline_parallel(IM,varargin)
% GetMidLine_parallel - Given an image series or an image dir containing an image
%   series, returns a series of indices corresponding to the midline of the
%   the fish in each image
% midlineInds = GetMidline(IM);
% midlineInds = GetMidline(..., fishPos);
% midlineInds = GetMidline(IM,fishPos,lineLens);
% midlineInds = GetMidline(IM,fishPos,lineLens, imgExt);
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
% imgExt - Image extension (default: 'jpg')
% Outputs:
% midlineInds - Indices of the midline in the image series; L x T matrix
%   where T is the number of images in the series and L is number of line
%   indices. L = sum(lineLens).

heights =  [18 16 14 10 8 8];
imgExt = 'jpg';
dTh = 4;
if nargin == 1
    if ischar(IM) && isdir(IM)
        IM  = ReadImgSequence(IM,imgExt);
    end
    fishPos = GetFishPos(ProcessImages(IM),50);
elseif nargin == 2
    fishPos = varargin{1};
elseif nargin ==3
    fishPos = varargin{1};
    heights = varargin{2};
elseif nargin ==4
    fishPos = varargin{1};
    heights = varargin{2};
    imgExt = varargin{3};
end

if isempty(fishPos)
    fishPos = GetFishPos(ProcessImages(IM),50);
end

if matlabpool('size')==0
    poolSize = 10;
    matlabpool(poolSize);
end
imgInds = 1:size(IM,3);
midlineInds = cell(1,length(heights));
lineEndInds = nan(size(fishPos));
tic  
for seg = 1:length(heights)
    disp(['Obtaining orientation of segment ' num2str(seg) '...'])    
    lineInds = cell(size(imgInds));      
    if seg ==1
        startPts_now = fishPos;
        startPts_prev = nan(size(fishPos)) ;
    else
        startPts_prev = startPts_now;
        startPts_now = lineEndInds;
    end
    height = heights(seg);
    dispChunk = ceil(size(IM,3)/5);
    parfor imgNum = imgInds
        if mod(imgNum,dispChunk)==0
            disp(imgNum)
        end
        img = IM(:,:,imgNum);
        img = max(img(:))-img;
        lineInds{imgNum} = GetML(img,startPts_now(imgNum,:), startPts_prev(imgNum,:), dTh, height);
        [row,col] = ind2sub(size(img),lineInds{imgNum}(end));
        lineEndInds(imgNum,:) = [col,row];          
    end     
    midlineInds{seg} = lineInds;
end
midlineInds = Reconfigure(midlineInds);
toc
% matlabpool close
end


function lineInds = GetML(im,startPt,varargin)
% lineInds = GetMidline(im,startInd,prevStartInd,dTh,lineLen)

dTh = 4;
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
    dTh  = 4;
end

[rImg,indMat] = RadialFish(im,startPt,dTh,lineLen);
ker = gausswin(round(dTh))*gausswin(round(lineLen/4))'; ker = ker/sum(ker(:));
rImg = conv2(rImg,ker,'same');

if (isempty(prevStartPt)) || (any(isnan(prevStartPt)))
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
Standardize = @(x)(x-min(x))/(max(x)-min(x));
muPxls = mean(rImg,2);
backgroundInt = mean(muPxls);
signalMat = rImg > 1*backgroundInt;
rImg(rImg<0.5*backgroundInt)=min(rImg(:));
nPxls = sum(signalMat,2);
[lps,mr] = GetLineProfileSpread(rImg);
nml = (nPxls(:).^1.5).*muPxls(:).*lps(:);
nml(isinf(nml))= max(nml);
nml = nml(farInds);
mr = mr(farInds);
nml = Standardize(nml);
mr = Standardize(mr);
mr = nml;
probInds = find(nml>0.5);
nml(probInds) = nml(probInds).*mr(probInds);
thr = 0.5;
probInds= find(nml>thr);
count = 0;
while (numel(probInds)<2) && (count <10)
    thr = thr*0.9;
    probInds = find(nml > thr);
    disp('Lowering threshold to find segment...')
    count = count + 1;
end
if isempty(probInds)
    blahInds = farInds; % At the moment, not really dealing with a segment not being found!
else
    blahInds = farInds(probInds);
end
nml = nml(probInds);

%## Find lines that are not contiguous blocks (i.e. islands) and eliminate
[blockSizes,blockInds] = GetContiguousBlocks(blahInds);

%##Do not uncomment these lines, this could give an erorr
%###############################################
% blockInds(blockSizes==1)=[]; 
% blockSizes(blockSizes==1)=[];
%###############################################

blockEndInds = blockInds + blockSizes - 1;
blockMaxes = zeros(size(blockInds));
for kk = 1:length(blockInds)
    blockMaxes(kk) = max(nml(blockInds(kk):blockEndInds(kk)));
end

%## For 1st line segment choose smaller block because head block will be
%## bigger than tail block

%# Choose a block that has at at least a few lines and high max int, but weight max int more.
[~,bigBlock]= max(blockSizes.*(2*blockMaxes)); 

keepInds = blockInds(bigBlock):blockEndInds(bigBlock);
blahInds = blahInds(keepInds);
nml = nml(keepInds);
zerInds = find(blahInds==0);
blahInds(zerInds)=[];
nml(zerInds) = [];
ctrOfMassInd = round(sum((1:length(nml)).*nml(:)')/sum(nml));

realInd = blahInds(ctrOfMassInd);
lineInds = indMat(realInd,:)';


function [blockSizes, blockInds] = GetContiguousBlocks(values)
        % Given a set of values, returns the sizes of contiguous blocks and the
        %   block start indices
        blockInds = 1;
        blockSizes = [];
        count = 1;
        for kk = 1:length(values)-1
            if (values(kk+1)-values(kk))==1
                count = count+1;
            else
                blockSizes = [blockSizes; count];
                count = 1;
                blockInds =[blockInds; kk+1];
            end
        end
        blockSizes = [blockSizes; count];
  end
  
end

function midlineInds = Reconfigure(midlineInds)
% Reconfigures midlineInds to make compatible with
%   GetFishOrientationFromMidlineInds()
mlInds =cell(length(midlineInds{1}),1);
for imgNum  = 1:length(midlineInds{1})
    segInds = cell(length(midlineInds),1);
    for seg = 1:length(midlineInds)
        segInds{seg} = midlineInds{seg}{imgNum};
    end
    mlInds{imgNum} = segInds;
end
midlineInds = mlInds;
end

  function [lps,mr] = GetLineProfileSpread(rImg)
  Y =  sort(rImg,2,'descend');  
  blah = mean(Y(:,1:round(size(rImg,2)/5)),2);
  lps = sum(rImg,2)./blah; 
  
  Y = rImg';
  X = [ones(size(Y,1),1), [1:size(Y,1)]'];  
  B = X\Y;
  m = abs(B(2,:));
  b = B(1,:) ;
  Y_est = X*B;  
  res = sqrt(sum((Y_est-Y).^2,1));
  m = max(m)-m;
  res =max(res)-res;
  mr = m(:).*res(:);
  end
