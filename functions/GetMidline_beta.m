function midlineInds = GetMidline_beta(IM,varargin)
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
dTh = 4;
imgExt = 'jpg';
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

midlineInds = {};
for imgNum = 1:size(IM,3)
    img = IM(:,:,imgNum);
    img = max(img(:))-img;
    lineInds = {};
    startPt = fishPos(imgNum,:);
     
    [lineInds_all,optMat,indMat,parentMap] = GetMLs(img,startPt,dTh,heights); 
 
%     lineInds = cell(size(optMat,1),1);
%     lineInds{1} = lineInds_all{1}(:,indMat{1});
%     for seg = 2:length(optMat)  
%         ind = indMat{seg};
%         lineInds{seg} = lineInds_all{seg}(:,ind);       
%     end

   lineInds = cell(length(parentMap{end}),1);

    midlineInds{imgNum} = lineInds;     
    cla
    for kk = 1:length(lineInds)
        img(lineInds{kk}) = max(img(:));
        imagesc(img),axis image
    end
    hold on
    plot(fishPos(imgNum,1),fishPos(imgNum,2),'ko')
    title(num2str(imgNum))
    shg
    pause(0.15)
end
end

function [lineInds_all,optMat,indMat,parentMap] = GetMLs(im,fishPos,dTh,lineLens)
startPt = fishPos;
prevStartPt  =[];
lineInds_all = cell(numel(lineLens),1);
numMap = lineInds_all;
parentMap = numMap;
lineInds_all{1} = GetML(im,startPt,prevStartPt, dTh, lineLens(1));
numMap{1} = size(lineInds_all{1},2);
for kk = 1:size(lineInds_all{1},2);
    parentMap{1}{kk} = num2str(kk);
end
for seg = 2:numel(lineLens)
    blah = lineInds_all{seg-1};
    lineInds_sub =[];
    numMap{seg}=[];
    count = 0;
    for ln = 1:size(blah,2)        
        lInds = blah(:,ln);
        startPt = IndToSub(im,lInds(end));
        prevStartPt = IndToSub(im,lInds(1));
        temp = GetML(im,startPt,prevStartPt,dTh,lineLens(seg));
        for kk = 1:size(temp,2)
            count = count + 1;
            parentMap{seg}{count} = [parentMap{seg-1}{ln} num2str(kk)];
        end
        numMap{seg} = [numMap{seg}, size(temp,2)];
        lineInds_sub = [lineInds_sub, temp];        
    end
    numMap{seg}  = cumsum(numMap{seg});
    lineInds_all{seg} = lineInds_sub;
end
[optMat,indMat] = BestLineInds(lineInds_all,numMap);

    function [optMat,indMat] = BestLineInds(lineInds_all,numMap)
        optMat = cell(size(lineInds_all));
        indMat = optMat;              
        for ln1 = 1:numel(lineInds_all(1))
            l1 = im(lineInds_all{1}(:,ln1));
%             lineGrad = 0.2*(1:size(l1,1));
            lineGrad = ones(size(l1));          
            muPxls = repmat(lineGrad(:),1,size(l1,2));
            muPxls = mean(muPxls.*l1,1);
            lps = GetLineProfileSpread(l1');
            nml = muPxls(:).*lps(:);
            optMat{1}(ln1) = nml;            
        end
        for seg = 2:numel(lineInds_all)
            lineAppends = [];
            count = 0;
            for ln1 = 1:numel(numMap(seg-1))
                if ln1 ==1
                    inds = 1:numMap{seg}(ln1);
                else
                    inds = numMap{seg}(ln1-1)+1:numMap{seg}(ln1);
                end
                count = count + 1;                
                for ln2 = inds(:)'
                    %                     l1 = im(lineInds_all{seg-1}(:,ln1));
                    l2  = im(lineInds_all{seg}(:,ln2));
                    line_conc = [l1;l2];
                    lineGrad = 0.2*(1:size(line_conc,1));
                    muPxls = repmat(lineGrad(:),1,size(line_conc,2));
                    muPxls = mean(muPxls.*line_conc,1);
                    lps = GetLineProfileSpread(line_conc');
                    nml = muPxls(:).*lps(:);
                    optMat{seg}(ln2) = nml;                   
                end
                [~,indMat{seg}(ln1)] = max(optMat{seg});
            end
        end
        [~,indMat{1}] = max(optMat{1});
    end
    function sub  = IndToSub(im,ind)
        [r,c] = ind2sub(size(im),ind);
        x = c; y = r;
        sub = [x,y];
    end
end


function lineInds = GetML(im,startPt,varargin)
% lineInds = GetMidline(im,startInd,prevStartPt,dTh,lineLen)

dTh = 4;
lineLen = 15;
grad = 0.2;
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
    farInds = find(abs(dAngles)<=130);
end


%## Find lines that pass through fish
Standardize = @(x)(x-min(x))/(max(x)-min(x));
lineGrad = grad.*(1:size(rImg,2));
muPxls = mean(rImg.*repmat(lineGrad,size(rImg,1),1),2);
backgroundInt = mean(muPxls);

[lps,~] = GetLineProfileSpread(rImg);
nml = muPxls(:).*lps(:);
nml = nml(farInds);
nml = Standardize(nml);

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
if isempty(blockSizes)
    blockSizes = 1;
    blockInds = 1;
end
blockEndInds = blockInds + blockSizes - 1;

if numel(blockSizes)>1
    a = 1;
end
comInds = nan(size(blockSizes));
for blk = 1:numel(blockSizes)
    blkInds = blockInds(blk):blockEndInds(blk);
    comInds(blk) = blahInds(round(sum(blkInds(:).*nml(blkInds))/sum(nml(blkInds))));
end

lineInds = indMat(comInds,:)';


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
    function C = CorrC(mtrx, vctr)
        C = nan(size(mtrx,1),1);
        for row = 1:size(mtrx,1)
            blah = corrcoef(mtrx(row,:),vctr(:)');
            C(row) = blah(2);
        end
    end
end

function [lps,mr] = GetLineProfileSpread(rImg)
% Y =  sort(rImg,2,'descend');
Y = sort(rImg,2,'ascend');
blah = mean(Y(:,1:round(size(rImg,2)/2)),2);
lps = blah/sum(rImg(:));

Y = rImg';
X = [ones(size(Y,1),1), [1:size(Y,1)]'];
B = X\Y;
m = abs(B(2,:));
% b = B(1,:) ;
Y_est = X*B;
res = sqrt(sum((Y_est-Y).^2,1));
m = max(m)-m;
res =max(res)-res;
mr = m(:).*res(:);
end
