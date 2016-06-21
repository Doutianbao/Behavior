function midlineInds = GetMidlines(IM,varargin)
% GetMidLines - Given an image series or an image dir containing an image
%   series, returns a series of indices corresponding to the midline of the
%   the fish in each image
% midlineInds = GetMidline(IM);
% midlineInds = GetMidline(..., fishPos);
% midlineInds = GetMidline(IM,fishPos,lineLens);
% midlineInds = GetMidline(IM,fishPos,lineLens, imgExt,'ref',refImg);
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
% 'extraArenaInds', extraArenaInds - Image indices to be treated as outside
%   the arena and set to low intensity value.
%   values
% Outputs:
% midlineInds - Indices of the midline in the image series; L x T matrix
%   where T is the number of images in the series and L is number of line
%   indices. L = sum(lineLens).
%
% Avinash Pujala, Koyama lab/HHMI, 2016

heights =  [18 16 14 10];
dTh = 1;
imgExt = 'jpg';
extraArenaInds = [];
procType = 'serial';

if nargin == 1
    if ischar(IM) && isdir(IM)
        IM  = ReadImgSequence(IM,imgExt);
    end
    fishPos = GetFishPos(ProcessImages(IM),50);
elseif nargin == 2
    fishPos = varargin{1};
elseif nargin == 3
    fishPos = varargin{1};
    heights = varargin{2};
elseif nargin > 3
    fishPos = varargin{1};
    heights = varargin{2};
    imgExt = varargin{3};
end

for jj =  3:numel(varargin)
    if isstr(varargin{jj})
        switch lower(varargin{jj})
            case 'extraarenainds'
                extraArenaInds = varargin{jj+1};
            case 'proctype'
                procType = varargin{jj+1};
        end
    end
end

if isempty(fishPos)
    fishPos = GetFishPos(ProcessImages(IM),50);
end
midlineInds = cell(size(IM,3),1);

imgInds = 1:size(IM,3);
%# Getting outside arena points to minimize drawing of midline segments on them
if ~isempty(extraArenaInds)
    disp('Getting extra-arena points...')
%     rp = randperm(length(imgInds));
%     inds = imgInds(rp(1:min(length(rp),100)));
%     IM_sub = IM(:,:,inds);
%     minInt = min(IM_sub(:));
%     filtSize = ceil(mean([size(IM,1) size(IM,2)])/50);
%     ref = bfilter2(Standardize(ref),filtSize,[filtSize/2 0.5]);
%     extraArenaInds = ref <= 0.15;    
    
    rp = randperm(length(imgInds));
    inds = imgInds(rp(1:min(length(rp),100)));
    IM_sub = IM(:,:,inds);
    minInt = min(IM_sub(:));
    
else
    minInt =1;
end

if strcmpi(procType,'parallel')
    if matlabpool('size')==0
        matlabpool(10)
    end
    disp('Img # ')
    parfor imgNum = imgInds;
        img = IM(:,:,imgNum);
        img(extraArenaInds) = minInt;
        
        [lineMat, ~] = GetMLs(img,fishPos(imgNum,:),dTh,heights);
        
        [midlineInds{imgNum},~] = GetBestLine(img,lineMat,heights);
        
        disp(imgNum)
    end
else
    for imgNum = imgInds;
        img = IM(:,:,imgNum);
        img(extraArenaInds) = minInt;
        
        [lineMat, ~] = GetMLs(img,fishPos(imgNum,:),dTh,heights);
        
        [midlineInds{imgNum},~] = GetBestLine(img,lineMat,heights);
        
        PlotLineInds(img,fishPos(imgNum,:),midlineInds{imgNum},imgNum)
    end
end

end

%## Helper functions

function lineInds = GetML(im,startPt,varargin)
% lineInds = GetMidline(im,startInd,prevStartPt,dTh,lineLen)

dTh = 1;
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
    dTh  = 1;
end

[rImg,indMat] = RadialFish(im,startPt,dTh,lineLen);
ker = gausswin(round(dTh)*4)*gausswin(round(lineLen/4))'; ker = ker/sum(ker(:));
rImg = conv2(rImg,ker,'same');
allInds = 1:size(indMat,1);
if isempty(prevStartPt)
    farInds = allInds;
    nearInds = [];
else
    [y2,x2] = ind2sub(size(im),indMat(:,end));
    postVecs = [x2, y2] - repmat(startPt,length(x2),1);
    postVecs = postVecs(:,1) + postVecs(:,2).*1i;
    preVec = startPt-prevStartPt;
    preVecs = repmat(preVec(1) + preVec(2)*1i,length(x2),1);
    dAngles = angle(preVecs.*conj(postVecs)) *(180/pi);
    farInds = find(abs(dAngles)<=140);
    nearInds = setdiff(allInds,farInds);
end


%## Finding candidate midlines
nml = rImg2nml(rImg);
thr = 0.1;
maxtab(:,1) = GetPks(nml,'peakThr',thr,'thrType','abs','polarity',1);
maxtab(:,2) = nml(maxtab(:,1));

%## This next bit of code finds peaks that may not have been detected
%## because 'rImg' is linear when it should have been circular
midPt = round(length(nml)/2);
pks_shift = GetPks(circshift(nml(:),midPt),'peakThr',thr,'thrType','abs','polarity',1);
pks_bool = zeros(length(nml),1);
pks_bool(pks_shift) = 1;
pks_bool = circshift(pks_bool(:), -midPt);
pks_shift = find(pks_bool);
pks_shift = union(maxtab(:,1),pks_shift);
maxtab=[];
maxtab(:,1) = pks_shift;
maxtab(:,2) = nml(maxtab(:,1));
%##########

count = 0;
while (isempty(maxtab)) && (count <10)
    thr = thr*0.9;
    [maxtab,~] = peakdet(nml,thr);
    disp('Lowering threshold to find segment...')
    count = count + 1;
end

if isempty(maxtab)
    [maxtab(:,2), maxtab(:,1)] = max(nml);
end
[~,remInds] = intersect(maxtab(:,1),nearInds);
maxtab(remInds,:)=[];
remInds = find(maxtab(:,2) < 0.4*max(maxtab(:,2)));
maxtab(remInds,:)=[];

nPts = round((2/dTh));
probInds = GetPeriPts(maxtab(:,1),nPts);
probInds(probInds<0) = 1;
probInds(probInds > length(nml))= length(nml);
probInds(probInds==0)=1;
probInds = unique(probInds);
keepInds = intersect(probInds,farInds);
if isempty(keepInds)
    blahInds = farInds; % At the moment, not really dealing with a segment not being found!
    nml = nml(farInds);
else
    blahInds = keepInds;
    nml = nml(keepInds);
end

%## Find lines that are not contiguous blocks (i.e. islands) and eliminate
[blockSizes,blockInds] = GetContiguousBlocks(blahInds);
if isempty(blockSizes)
    blockSizes = 1;
    blockInds = 1;
end
blockEndInds = blockInds + blockSizes - 1;

comInds = nan(size(blockSizes));
for blk = 1:numel(blockSizes)
    blkInds = blockInds(blk):blockEndInds(blk);
    if sum(nml(blkInds))~=0
        comInds(blk) = blahInds(round(sum(blkInds(:).*nml(blkInds))/sum(nml(blkInds))));
    end
end
comInds(isnan(comInds))=[];
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

    function periPts = GetPeriPts(pts,nPeri)
        periPts =[];
        for jj = 1:length(pts)
            blah = pts(jj)-nPeri:pts(jj)+nPeri;
            periPts = [periPts;blah(:)];
        end
    end

end


function [lineMat,parentMap] = GetMLs(im,fishPos,dTh,lineLens)
startPt = fishPos;
prevStartPt  =[];
lineInds = cell(numel(lineLens),1);
parentMap = lineInds;
lineInds_first = GetML(im,startPt,prevStartPt, dTh, lineLens(1));
for kk = 1:size(lineInds_first,2);
    parentMap{1}{kk} = num2str(kk);
    lineInds{1}{kk} = lineInds_first(:,kk);
end
for seg = 2:numel(lineLens)
    count = 0;
    for ln = 1:length(lineInds{seg-1})      
        lInds = lineInds{seg-1}{ln};
        startPt = IndToSub(im,lInds(end));
        prevStartPt = IndToSub(im,lInds(length(lInds)-lineLens(seg-1)+1));        
        temp = GetML(im,startPt,prevStartPt,dTh,lineLens(seg));      
        for kk = 1:size(temp,2)
            count = count + 1;
            parentMap{seg}{count} = [parentMap{seg-1}{ln} '_' num2str(kk)];
            lineInds{seg}{count} = [lineInds{seg-1}{ln}; temp(:,kk)];         
        end
    end
    blah = cell2mat(lineInds{seg});
    if size(blah,2)>1
        bestLine = GetBestLine_sub(im,blah);
    else
        bestLine = blah;
    end
    lineInds{seg} = {bestLine};
end
lineMat = cell2mat(lineInds{end});

    function sub  = IndToSub(im,ind)
        [r,c] = ind2sub(size(im),ind);
        x = c; y = r;
        sub = [x,y];
    end
end

function nml = rImg2nml(rImg)
% muPxls = abs((mean(rImg,2).*median(rImg,2).*mode(rImg,2)).^(1/3));
shift = round(size(rImg,1)/3);
triTemp = CreateTailTemplate(11,size(rImg,2),0);
triTemp = triTemp(6:8,:);
blah = conv2(rImg,triTemp,'same');
blah2 = circshift(conv2(circshift(rImg,[shift,0]),triTemp,'same'),[-shift,0]);
blah = blah + blah2;
muPxls = abs((mean(blah,2).* median(blah,2).*mode(blah,2)).^(1/3));
% muPxls2 = median(blah + blah2,2);
% muPxls3 = mode(blah + blah2,2);
rImg2 = sort(rImg,2,'descend');
% temp = abs(rImg2-repmat(mean(rImg2,2)*0.9,1,size(rImg2,2)));
temp = abs(rImg2-repmat(mean(rImg2(:))*0.9,size(rImg2,1),size(rImg2,2)));
% [~,comInds] = min(temp,[],2);
[~,comInds] = min(temp,[],2);
ker = gausswin(max(1,round(size(rImg,1)/30)));
comInds = conv2(comInds(:),ker(:),'same');
nml =  ((Standardize(muPxls(:))).*comInds(:)).^0.5; % Standardization before multiplication is important so as no to convert valleys with negative values into peaks
N = length(nml);
nml = [nml(:); nml(:); nml(:)];
nml = Standardize(conv2(nml, ker(:),'same'));
nml = nml(N+1:2*N);
end

function nml = rImg2nml_bestLine(rImg)
muPxls = mean(rImg,2);
[~,gof]  = GetLineProfileSpread(rImg);
muPxls2 = mean(log2(Standardize(rImg)+1),2);
rImg2 = sort(rImg,2,'descend');
temp = abs(rImg2-repmat(mean(rImg(:))*0.9,size(rImg,1),size(rImg2,2)));
comInds = mean(rImg,2);
nml = abs(Standardize((muPxls(:)).*gof(:).*Standardize(muPxls2(:)).*comInds(:)).^(1/4));
nml = Standardize(nml);
% nml = Standardize(Standardize(muPxls(:)).*comInds(:)); % Standardization before multiplication is important so as no to convert valleys with negative values into peaks

end

function PlotLineInds(img,fishPos,lineInds,imgNum)
cla
inds = [];
for kk = 1:length(lineInds)
    inds = [inds; lineInds{kk}(:)];
end
img(inds) = max(img(:))*1.5;
imagesc(img),axis image
hold on
plot(fishPos(1),fishPos(2),'ko')
title(num2str(imgNum))
shg
% pause()
end

function lineMat = GetBestLine_sub(img,lineMat)
lineProfiles = img(lineMat);
nml = rImg2nml_bestLine(lineProfiles');
[~,ind] = max(nml);
ind = ind(1);
lineMat = lineMat(:,ind);
end

function [midlines, lineMat] = GetBestLine(img,lineMat,heights)
lineProfiles = img(lineMat);
% muPxls = mean(lineProfiles,1);
% [lps,gof] = GetLineProfileSpread(lineProfiles');
% lps = lps';
% nml = muPxls.*lps.*gof;

nml = rImg2nml(lineProfiles');
[~,ind] = max(nml);
ind = ind(1);
lineMat = lineMat(:,ind);
midlines = cell(length(heights),1);
midlines{1} = lineMat(1:heights(1));
startInd = 1;
for seg = 2:numel(heights)
    startInd = startInd + heights(seg-1);
    segInds = startInd:+ startInd + heights(seg)-1;
    midlines{seg} = lineMat(segInds);
end
end


function [lps,gof] = GetLineProfileSpread(rImg)
Y = sort(rImg,2,'ascend');
blah = mean(Y(:,1:round(size(rImg,2)/5)),2);
% lps = blah/sum(rImg(:));
lps = blah;

Y = rImg';
X = [ones(size(Y,1),1), [1:size(Y,1)]'];
B = X\Y;
m = abs(B(2,:));
% b = B(1,:) ;
Y_est = X*B;
gof = sqrt(sum((Y_est-Y).^2,1))./var(Y,[],1);

end
