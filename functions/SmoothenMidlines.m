function tailCurv = SmoothenMidlines(midlineInds,imgStack, varargin)
%SmoothMidlines - Given midline indices generated get GetMidlines.m and the
%   stack of images from whence these indices arise, returns smoothed tail
%   curves (head segment omitted) a la Huang et al, 2013.
%
% tailCurv = SmoothMidlines(midlineInds,imgStack);
% tailCurv = SmoothMidlines(midlineInds,imgStack,nHood);
% tailCurv = SmoothMidlines(midlineInds,imgStack,nHood,'plotBool', plotBool,)
% Inputs:
% midlineInds - Midline segments as generated by GetMidlines.m
% imgStack - Stack of images to which the midline indices correspond
% nHood - Size of local pixel neighborhood used to smooth zig-zag tail
%   curves resulting from integer coordinates of midlineInds. Default = 1.
% 'dsVecs' - A cell array like midlineInds, but containing a set of
%   distances between each pixel of the midline and the head centroid
%   pixel(i.e., fish position). If dsVecs is not empty, will use this to
%   smoothen midlines
%
% Reference:
% Huang, K.-H., Ahrens, M.B., Dunn, T.W., and Engert, F. (2013). Spinal Projection
%   Neurons Control Turning Behaviors in Zebrafish. Current Biology 23, 1566�1573.
%
% Avinash Pujala, Koyama lab/HHMI, 2016

nHood = 2; % Default neighborhood size
plotBool = false;
pauseDur = 0;
smoothFactor = 2;
dsVecs = [];

if nargin < 2
    error('Minimum 2 inputs required!')
elseif nargin > 2
    if ~ischar(varargin{1})
        nHood = varargin{1};
    end
end
for jj = 1:numel(varargin)
    if ischar(varargin{jj})
        switch lower(varargin{jj})
            case 'plotbool'
                plotBool = varargin{jj+1};
            case 'pausedur'
                pauseDur = varargin{jj+1};
            case 'smoothfactor'
                smoothFactor = varargin{jj+1};
            case 'dsvecs'
                dsVecs = varargin{jj+1};
        end
    end
end
if length(midlineInds) ~= size(imgStack,3)
    error('Mismatch in size of inputs, check inputs!')
end
s1= zeros(length(dsVecs),1);
s2 = s1;
ds = s1;
if ~isempty(dsVecs)
    for jj = 1:length(dsVecs)
        if ~isempty(dsVecs{jj})
            s1(jj) = dsVecs{jj}(1);
            s2(jj) = dsVecs{jj}(end);
            ds(jj) = mean(diff(dsVecs{jj}));
        end
    end
end
s1(s1==0) = mode(s1);
s2(s2==0) = mode(s2);
ds(ds==0) = mode(ds);
distVec = median(s1):min(median(ds),1):median(s2);

imgStack = Standardize(imgStack); % Important for smoothing midlines using pxl intensity weighting (avoids -ve values)
if ~isempty(dsVecs)
    bodyLen = length(distVec);
else
    bodyLen = size(cell2mat(midlineInds{1}),1);
    tailLen = bodyLen -length(midlineInds{1}{1});
end

tailCurv = zeros(bodyLen,2,length(midlineInds));
dispChunk = round(length(midlineInds)/10);
% dispChunk = 1;
disp('Smoothing midlines...')
if plotBool
    figure('Name','Smoothed fish tail curvature')
end
tic

imgDims = size(imgStack);
imgDims = imgDims(1:2);

N = length(midlineInds);
for iNum = 1:N
    mlInds = midlineInds{iNum};
    if iscell(mlInds)
        if length(mlInds)>1
            mlInds = cell2mat(mlInds(1:end));
        else
            mlInds= cell2mat(mlInds(1));
        end
    end
    if ~isempty(dsVecs)
        tc = SmoothenMidline(mlInds,imgStack(:,:,iNum),nHood);
        if isempty(tc)
            tc = nan(size(tailCurv,1),size(tailCurv,2));
            tailCurv(:,:,iNum) = tc;
        else
            tc = SplineMidlineInds_ds(tc,dsVecs{iNum},distVec,smoothFactor,imgDims);
            tailCurv(:,:,iNum) = SplineTailCurv(tc,smoothFactor);
        end
    else
        mlInds = SplineMidlineInds(mlInds,smoothFactor,imgDims);
        tc = SmoothenMidline(mlInds,imgStack(:,:,iNum),nHood);
        tailCurv(:,:,iNum) = SplineTailCurv(tc,smoothFactor);
    end
    if mod(iNum,dispChunk)==0
        disp(['Img # ' num2str(iNum) '/' num2str(N)])
    end
    if plotBool
        cla
        img = imgStack(:,:,iNum);
        img(mlInds) = 0;
        imagesc(img),axis image, colormap(gray)
        hold on
        plot(size(imgStack,1)/2+1,size(imgStack,2)/2+1,'b*','markersize',10)
        plot(tailCurv(:,1,iNum), tailCurv(:,2,iNum),'r:','linewidth',1.5)
        %         plot(tc(:,1),tc(:,2),'g.-')
        drawnow
        title(['Img # ' num2str(iNum)])
        shg
        if isempty(pauseDur)
            pause()
        else
            pause(pauseDur)
        end
    end
end
toc

end


function tailCurv = SmoothenMidline(mlInds,img,nHood)
tailCurv = zeros(length(mlInds),2);
[r,c] = ind2sub(size(img),mlInds);
[C,R] = meshgrid(1:size(img,2),1:size(img,1));
imgDims = size(img);
for jj = 1:length(r)
    rInds = r(jj)-nHood:r(jj)+nHood;
    rInds(rInds<=0)=[];
    rInds(rInds>imgDims(1)) = [];
    cInds = c(jj)-nHood:c(jj)+nHood;
    cInds(cInds<=0)=[];
    cInds(cInds>imgDims(2))=[];
    try
        rNeighbors = R(rInds,cInds);
    catch
        rNeighbors = R(rInds,cInds);
    end
    cNeighbors = C(rInds,cInds);
    wts = img(rInds,cInds);
    %     ker = ones(length(wts));
    %     wts = conv2(wts,ker,'same');
    pxlInd(1) = sum(rNeighbors(:).*wts(:))/sum(wts(:));
    pxlInd(2) = sum(cNeighbors(:).*wts(:))/sum(wts(:));
    tailCurv(jj,:)= fliplr(round(pxlInd*10)/10); % Flipping to give in x-y rather than row-col coordinates
end
end


function tailCurv_spline = SplineTailCurv(tailCurv,smoothFactor)
% y = [tailCurv(1,:); tailCurv; tailCurv(end,:)];
y = tailCurv;
t = 1:size(y,1);
ts = linspace(1,size(y,1),size(y,1)/smoothFactor);
ys(:,1) = spline(t,y(:,1),ts);
ys(:,2) = spline(t,y(:,2),ts);
y(:,1) = spline(ts,ys(:,1),t);
y(:,2) = spline(ts,ys(:,2),t);
tailCurv_spline = y;
% tailCurv_spline = [xx(2:end-1); yy(2:end-1)]';

end

function mlInds_spline = SplineMidlineInds(mlInds,smoothFactor,imgDims)
[r,c] = ind2sub(imgDims,mlInds);
t = 1:length(r);
ts= t(1:smoothFactor:end);
rs = r(ts);
cs = c(ts);
r_spline = round(spline(ts,rs,t));
c_spline = round(spline(ts,cs,t));
mlInds_spline = sub2ind(imgDims,r,c);

end

function tc = SplineMidlineInds_ds(tailCurv,dsVec,distVec,smoothFactor,imgDims)
r = tailCurv(:,2);
c = tailCurv(:,1);
ts= distVec;
if length(dsVec) > length(distVec)
    rs = r(1:length(distVec));
    cs = c(1:length(distVec));
else
    t = linspace(distVec(1),distVec(end),min(length(distVec),length(dsVec)));
    rs = interp1(t,r,ts,'cubic');
    cs = interp1(t,c,ts,'cubic');
end

% rs(rs==0) = 1;
% cs(cs==0)=1;
% rs(rs>imgDims(1)) = imgDims(1);
% cs(cs>imgDims(2)) = imgDims(2);

t = 1:length(rs);
temp = t(2:smoothFactor:end-1);
ts = [t(1), temp, t(end)];
rs = rs(ts);
cs = cs(ts);
r_spline = spline(ts,rs,t);
c_spline = spline(ts,cs,t);
% mlInds_spline = sub2ind(imgDims,r_spline,c_spline);
tc = [c_spline(:),r_spline(:)];

end




