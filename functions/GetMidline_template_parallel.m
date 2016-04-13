function midlineInds = GetMidline_template(IM,varargin)
% GetMidLine - Given an image series or an image dir containing an image
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
templateGradient = -0.2;
dTh = 5;
poolSize = 12;

imgExt = 'jpg';
if nargin == 1
    if isstr(IM) && isdir(IM)
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

if matlabpool('size') == 0
    matlabpool(poolSize)
end
imgInds = 1:size(IM,3);
% startPts_now = nan(size(fishPos));
% startPts_prev = startPts_now;
midlineInds = cell(1,length(heights));
lineEndInds = nan(size(fishPos));
tic  
for seg = 1:length(heights)
    disp(['Obtaining orientation of segment ' num2str(seg) '...'])
    base = min(1,round(heights(seg)/5));
    T = CreateTailTemplate(base,heights(seg),templateGradient);
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
        lineInds{imgNum} = GetML(img,startPts_now(imgNum,:), startPts_prev(imgNum,:), dTh, height, T);
        [row,col] = ind2sub(size(img),lineInds{imgNum}(end));
        lineEndInds(imgNum,:) = [row,col];          
    end     
    midlineInds{seg} = lineInds;
end
toc
end

function lineInds = GetML(im,startPt,varargin)
% lineInds = GetMidline(im,startInd,prevStartInd,dTh,lineLen,matchTemplate)

dTh = 4;
lineLen = 15;
if nargin < 2
    error('At least 2 inputs required!')
elseif nargin ==2
    prevStartPt = [];
    T = CreateTailTemplate(6,30,-0.15);
elseif nargin == 3
    prevStartPt = varargin{1};
    T = CreateTailTemplate(6,30,-0.15);
elseif nargin ==4
    prevStartPt = varargin{1};
    dTh = varargin{2};
    T = CreateTailTemplate(6,30,-0.15);
elseif nargin == 5
    prevStartPt = varargin{1};
    dTh = varargin{2};
    lineLen = varargin{3};
    T = CreateTailTemplate(6,30,-0.15);
elseif nargin ==6;
    prevStartPt = varargin{1};
    dTh = varargin{2};
    lineLen = varargin{3};
    T = varargin{4};
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
[mV,~] = RotateAndMatchTemplate(im,startPt,dTh,T);
[lps,~] = GetLineProfileSpread(rImg);
nml = lps(:).*mV(:).^2;
% nml = mV.^2;
nml = nml(farInds);

nml = Standardize(nml);
thr = 0.4;
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
if isempty(blockSizes)
    blockSizes = 1;
    blockInds = 1;
end

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
[~,bigBlock]= max((1./blockSizes).*(2*blockMaxes));

keepInds = blockInds(bigBlock):blockEndInds(bigBlock);


blahInds = blahInds(keepInds);

nml = nml(keepInds);
zerInds = find(blahInds==0);
blahInds(zerInds)=[];
nml(zerInds) = [];
ctrOfMassInd = round(sum((1:length(nml)).*nml(:)')/sum(nml));
[~, maxInd] = max(nml);

realInd = blahInds(ctrOfMassInd);
% realInd = blahInds(maxInd);
lineInds = indMat(realInd,:)';

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
                count = 1;
                blockInds =[blockInds; jj+1];
            end
        end
        blockSizes = [blockSizes; count];
    end

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
