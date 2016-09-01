
%% Read img stack
[file,path] = uigetfile('*.tif');
imgStack = ReadImgStack(fullfile(path,file),'tif');
imgStack_ab = imgStack;

%% Getting rois
nImages = size(imgStack,3);
imgInds = 9:83;
imgInds = imgInds(:)';
roi = cell(nImages,1);
for jj = imgInds
    img = imgStack(:,:,jj);
    roi{jj} = setEllipticalRois(img,jet(256));
end
img_maxInt = max(imgStack(:,:,imgInds),[],3);
% roi_maxInt = setEllipticalRois(img_maxInt,jet(256));
% for jj = imgInds(:)'
%     roi{jj} = roi_maxInt;roi
% end

%% Find background roi
disp('Background ROI...')
pause(3)
backRoi = setEllipticalRois(img_maxInt,jet(256));
backInds = [];
for jj = 1:length(backRoi)
    if ~isempty(backRoi{jj})
           backInds = [backInds; backRoi{jj}.idx(:)];
    end
 end
rp = randperm(numel(backInds));
backInds = backInds(rp);
n = nextpow2(numel(backInds))-1;
f =2^n;
backInds = backInds(1:f);
backImg = img_maxInt(backInds);
r = floor(n/2);
c  = n-r;
backImg = reshape(backImg,2^r,2^c);
while numel(backImg)< numel(img_maxInt)
    backImg = repmat(backImg,2,2);
end
nPxls = numel(backImg);
rp  = randperm(nPxls);
backImg(1:nPxls) = backImg(rp);
thr = mean(backImg(:)) + 7*std(backImg(:));
break

%% Adjusting
% mn = min(imgStack(:));
disp('Ablating...')
blah = [];
% imgInds = find(ones(size(imgStack,1),size(imgStack,2)));
for jj = 1:length(roi)
    disp(['Slice # ' num2str(jj)])
%     imgInds_shuf = randperm(nPxls);
    if ~isempty(roi{jj})
        img = imgStack_ab(:,:,jj);
        for rr = 1:length(roi{jj});
            if length(roi{jj})>1
                roi_current = roi{jj}{rr};
            else
                roi_current = roi{jj}{1};
            end
            if ~isempty(roi_current)
                disp(['Roi # ' num2str(rr)])
                inds = roi_current.idx;                
                vals = img(inds);
%                 sVals = zscore(vals);
%                 overInds = find(sVals>1.5);
                  overInds = find(vals>thr);
                  if ~isempty(overInds)
                      overInds= inds(overInds);
                  end
                  overInds = inds();
%                 if ~isempty(overInds)
%                 mu = mean(vals(:));
%                 sig = std(vals(:));
%                 mult = Standardize(vals(overInds)) + 0.1;
%                 blah = mu*ones(size(overInds)) + 2.5*(mult.*rand(size(overInds)))*std(vals(:));
%                 blah = blah(randperm(length(blah)));                
%                 vals(overInds) = blah;      
%                 end
%                 img(inds) = vals;
                  img(overInds) = backImg(overInds);
%                   medImg = mean(img(:));
%                   sigImg = std(img(:));
%                   brightInds = find(img(inds)>medImg + 10*sigImg(:));
%                   img(inds(brightInds)) = medImg*ones(size(brightInds)) + sigImg*rand(size(brightInds));
            end
        end
        imgStack_ab(:,:,jj) = img;
    end
end


%% Noise
disp('Noising...')
sig = std(imgStack_ab(:));
sig = repmat(sig,size(imgStack));
imgStack_ab = imgStack_ab + 0.2*sig.*randn(size(imgStack)); 


%% Standardizing and Saving
disp('Standardizing...')
% for jj = 1:size(imgStack,3)
%     imgStack(:,:,jj) = Standardize(imgStack(:,:,jj));
%     imgStack_ab(:,:,jj) = Standardize(imgStack_ab(:,:,jj));
% end
% imgStack = Standardize(imgStack);
imgStack_ab = Standardize(imgStack_ab);

