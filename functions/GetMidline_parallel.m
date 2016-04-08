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

lineLens =  [18 16 14 10 8 8];
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
    lineLens = varargin{2};
elseif nargin ==4
    fishPos = varargin{1};
    lineLens = varargin{2};
    imgExt = varargin{3};
end

if isempty(fishPos)
    fishPos = GetFishPos(ProcessImages(IM),50);
end

if matlabpool('size')==0
    poolSize = 10;
    matlabpool(poolSize);
end
midlineInds = {};
imgNumVec = 1:size(IM,3);
parfor imgNum = imgNumVec
    img = IM(:,:,imgNum);
    img = max(img(:))-img;
    lineInds = {};
    startPt = fishPos(imgNum,:);
    for jj = 1:length(lineLens)
        if jj ==1
            
            lineInds{jj} = GetML(img,startPt,[],[],lineLens(jj));
        else
            lineInds{jj} = GetML(img,startPt,prevStartPt,[],lineLens(jj));
        end
        si = lineInds{jj}(end);
        [r,c] = ind2sub(size(img),si);
        x = c;
        y = r;
        prevStartPt = startPt;
        startPt = [x,y];
        
        midlineInds{imgNum} = lineInds;
    end
    
    cla
    for kk = 1:length(lineInds)
        img(lineInds{kk}) = max(img(:));
        imagesc(img),axis image
    end
    hold on
    plot(fishPos(imgNum,1),fishPos(imgNum,2),'k*')
    title(num2str(imgNum))
    shg
    pause(0.15)
end
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
ker = gausswin(round(1))*gausswin(round(lineLen/4))'; ker = ker/sum(ker(:));
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
    farInds = find(abs(dAngles)<=90);
end


%## Find lines that pass through fish
Standardize = @(x)(x-min(x))/(max(x)-min(x));
muPxls = mean(rImg,2);
backgroundInt = mean(muPxls);
signalMat = rImg > 1*backgroundInt;
rImg(rImg<backgroundInt)=min(rImg(:));
nPxls = sum(signalMat,2);
[lps,mr] = GetLineProfileSpread(rImg);

nml = (nPxls(:).^1.5).*muPxls(:).*lps(:);
nml(isinf(nml))= max(nml);
nml = Standardize(nml);
nml = nml(farInds);
probInds = find(nml>0.4);
mr = Standardize(mr);
nml(probInds) = nml(probInds).*mr(probInds);
probInds= find(nml>0.4);
blahInds = farInds(probInds);
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
for jj = 1:length(blockInds)
    blockMaxes(jj) = max(nml(blockInds(jj):blockEndInds(jj)));
end

%## For 1st line segment choose smaller block because head block will be
%## bigger than tail block

%# Choose a block that has at at least a few lines and high max int, but weight max int more.
[~,bigBlock]= max(blockSizes.*(2*blockMaxes)); 

if numel(blockInds)> bigBlock
    keepInds = blockInds(bigBlock):blockInds(bigBlock+1);
else
    keepInds = blockInds(bigBlock):blockEndInds(bigBlock);
end

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
