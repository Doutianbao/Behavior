function varargout = Karina(varargin)
% Karina - Function for finding fish and fixing position and orientation
% blah = Karina()
% blah = Karina(filePath);

%% Some variables
frameRate  = 500; % For timeseries lowpass filtering
lowPass_ts = 100; % Lowpass cutoff for timeseries filtering
bandPass_img = [5, 10]; % For image spatial filtering
intThr = 0.75; % For image binarization
lineLen = 25; % For orienttion estimation


%% Read Image Stack
if nargin==0
    [fName,fDir] = uigetfile('*');
    fPath = fullfile(fDir,fName);
    dotInd = strfind(fName,'.');
    imgExt = fName(dotInd+1:end);
else
    fPath = varargin{1};
    [fDir,fName,imgExt] = fileparts(fPath);
end

IM = ReadImgStack(fPath,imgExt);
IM = max(IM(:))-IM;

%% Filter image stack
% IM_smooth = SmoothImageStack(IM,10);bwc
disp('Filtering and grayscaling images...')
IM_flt = BandpassStack(IM,bandPass_img(1),bandPass_img(2));
IM_bw = imNormalize999(IM_flt);
IM_bw(IM_bw<intThr) = 0;

nFrames = size(IM,3);
blobs = bwconncomp(IM_bw(:,:,1));
rp = regionprops(blobs,'centroid','Area');
blobPos = cell(size(rp));
for blob = 1:length(blobPos)
    blobPos{blob} = nan(nFrames,3);
    blobPos{blob}(1,1:2) = rp(blob).Centroid;
    blobPos{blob}(1,3) = 1;
end
disp('Matching centroid positions frame-by-frame...')
for frame = 2:nFrames
    %# Centroids in current img
    rp = regionprops(bwconncomp(IM_bw(:,:,frame)),'Centroid');
    cents = nan(length(rp),2);
    for cent = 1:length(rp)
        cents(cent,:) = rp(cent).Centroid;
    end
    %# Centroids in prev img
    cents_prev = nan(length(blobPos),2);
    for cent = 1:size(cents_prev,1)
        cents_prev(cent,:) = blobPos{cent}(frame-1,1:2);
    end
    distMat = GetDistMat(cents,cents_prev);  
    optInds = GetOptimalInds(distMat);  
    sInds = optInds;    
    for rp = 1:size(cents_prev,1)
        blobPos{rp}(frame,1:2)= cents(sInds(rp),:);
        blobPos{rp}(frame,3) = round(frame);
    end
end

%% Getting orientations of fish and tranformed frames
disp('Getting orientations and transformed image stacks...')
orientations = cell(size(blobPos));
IM_adj = orientations;
dOr = orientations;
imgDims = size(IM);
if matlabpool('size')==0
    matlabpool(10)
end
for or = 1:length(orientations)   
    mlInds = GetMidline_parallel(-IM,blobPos{or}(:,1:2),lineLen);
    orientations{or} = chebfilt(GetFishOrientationFromMidlineInds(mlInds,imgDims(1:2)),1/frameRate,lowPass_ts,'low');
    dOrs{or} = DiffOrientation(orientations{or})*180/pi;
    IM_adj{or} = PlayFixedFish(IM,blobPos{or}(:,1:2),orientations{or}+180);
end


if matlabpool('size')>1
    matlabpool close
end

out = struct;
out.fishPos = blobPos;
out.orientation = orientations;
out.dOr = dOrs;
out.IM_adj = IM_adj;
varargout{1} = IM;
varargout{2} = out;



%% Utility functions
    function voteMat = DistMat2VoteMat(distMat)
        voteMat = nan(size(distMat));
        for col = 1:size(distMat,2)
            for row = 1:size(distMat,1)
                blah = distMat;
                blah(row,:)=[];
                blah(:,col)  =[];
                voteMat(row,col) = sum(blah(:));
            end
        end
    end

    function optInds = GetOptimalInds(distMat)
        [sorts,inds] = sort(distMat,'ascend');
        rowInds = []; colInds = [];
        if size(inds,2) > 1
            for c = 1:size(inds,2)-1
                sameInds = find(inds(1,:)==inds(1,c));
                if numel(sameInds)>1
                    addVec = 1:length(sameInds);
                    blah = nan(length(addVec),1);
                    temp = {};
                    for jj = 1:length(addVec)
                        addVec = circshift(addVec,[0,jj-1]);
                        temp{jj} = addVec;
                        blah(jj) = sum(SparseIndex(sorts,addVec,sameInds));
                    end
                    [~, mindInd] = min(blah);
                    rowInds = [rowInds,temp{mindInd}]
                    colInds = [colInds,sameInds];
                    optInds = SparseIndex(inds,rowInds,colInds);
                else
                    optInds = inds(1,:);
                end
            end
        else
            optInds= inds(1,1);
        end
    end

    function out = GetDistMat(m1,m2)
        S = @(x,y) sqrt(sum((x(1)-y(1)).^2 + (x(2)-y(2)).^2));
        out = nan(size(m1,1),size(m2,1));
        for one = 1:size(out,1)
            for two = 1:size(out,2)
                out(one,two) = S(m1(one,:), m2(two,:));
            end
        end
    end

    function IM_flt = BandpassStack(IM,lowerLim,upperLim)
        IM_flt = zeros(size(IM));
        for imNum = 1:size(IM,3)
            IM_flt(:,:,imNum) = Bandpass(IM(:,:,imNum),lowerLim,upperLim);
        end
    end

    function img_flt = Bandpass(img, lowerLim,  upperLim)
        getCtr = @(img)[round(size(img,1)/2+0.599), round(size(img,2)/2+0.599)];
        imgCtr = getCtr(img);
        r_in = lowerLim;
        r_out = upperLim;
        se_in = strel('disk',r_in);
        se_out = strel('disk',r_out);
        se_in_ctr = getCtr(se_in.getnhood);
        se_out_ctr = getCtr(se_out.getnhood);
        mask = zeros(size(img));
        inInds = se_in.getneighbors;
        inInds(:,1) = inInds(:,1) + imgCtr(1);%-se_in_ctr(1);
        inInds(:,2) = inInds(:,2) + imgCtr(2);%-se_in_ctr(2);
        inInds = sub2ind(size(img),inInds(:,1),inInds(:,2));
        outInds = se_out.getneighbors;
        outInds(:,1) = outInds(:,1) + imgCtr(1);%-se_out_ctr(1);
        outInds(:,2) = outInds(:,2) + imgCtr(2);%-se_out_ctr(2);
        outInds = sub2ind(size(img),outInds(:,1),outInds(:,2));
        mask(outInds) = 1;
        mask(inInds) =0;
        img_fft = fftshift(fft2(img));
        img_fft = img_fft.*mask;
        img_flt = ifft2(fftshift(img_fft));
    end

    function IM = ReadImgStack(fPath,imgExt)
        imgInfo = imfinfo(fPath);
        nImages = length(imgInfo);
        IM = zeros(imgInfo(1).Height,imgInfo(1).Width,nImages);
        dispChunk = round(nImages/3);
        disp('Reading images...')
        for imNum = 1:size(IM,3)
            IM(:,:,imNum) = imread(fPath,'tif', imNum);
            if mod(imNum,dispChunk)==0
                disp(num2str(imNum))
            end
        end
    end


    function IM_smooth = SmoothImageStack(IM, kerSize)
        % ker = ones(10,10); ker = ker/sum(ker);
        ker = gausswin(kerSize);
        ker = ker(:)*ker(:)';
        ker = ker/sum(ker);
        IM_smooth = zeros(size(IM));
        for imNum = 1:size(IM_smooth,3)
            IM_smooth(:,:,imNum) = conv2(IM(:,:,imNum),ker,'same');
        end
    end


end
