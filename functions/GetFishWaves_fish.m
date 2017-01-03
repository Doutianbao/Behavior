function varargout = GetFishWaves_fish(procData,varargin)
%GetFishWaves_fish Given procData.mat (created by FishSwim) or path to it
% returns wavelet transforms for head and tail orienations for all trials
% in addition to plotting these if specified
% W = GetFishWaves_fish([]);
% W = GetFishWaves_fish(procData);
% W = GetFishWaves_fish(pathToProcData);
% W = GetFishWaves_fish(.., 'hr',headRange,'tr',tailRange,'nFramesInTrl',nFramesInTrl,'fps',fps,...
%       'freqRange',freqRange,'dj',dj,'noiseType',noiseType,'stringency',stringency,'freqScale','lin',...
%       'plotOrNot', plotOrNot,'trlList',trlList,'xLim',xLim,'cLim',cLim);
% Inputs:
% 'hr' - Head range for getting head orientation. For instance, inputting
%   [0 25] (default) results in using the beginning 25 of the fish to determine head
%   orientation
% 'tr' - Tail range (default = [75 100]).
% 'nFramesInTrl' - Number of frames in a trial (default = 750)
% 'fps' - Frames per second (default = 500)
% 'freqRange' - Freq range for wavelet transformation (WT) (default = [10 70])
% 'freqScale' - 'Log' or 'lin'; The former results in log2 frequency
%   scales for WT, whereas the latter results in linear scales (default = 'lin')
% 'noiseType' - 'red' or 'white', for statistical significance testing of
%   the WT (Default = 'red').
% 'stringency' - Stringency for significance testing (Default = 0).
% 'sigmaXY' - Divisor for normalization of coefficients in wavelet
%   transform. Can be input for group data. If empty, computes
%   automatically
% 'plotOrNot' - 0 or 1; the latter results in plotting
% 'trlList' - The list of trials to plot. If trList = [], then plots all
%   trials
% 'stimTime' - Stim onset time in milliseconds (default = 100);
% 'xLim' - X limits of plots (default = [-50 650], in milliseconds, customized for escape responses)
% 'cLim' - Color limits for normalized wavelet plots (default = [0.1 3]).
% onsetAlign = 0 or 1. If 1, aligns timeseries of different trials based
%   onset of response (stored in procData.elicitedSwimInfo). If 0, aligns
%   w.r.t stimulus frame
% traceType - 'headTail','curv','both'. Determines the type of timeseries
%   on which to compute WTs.
%
% Avinash Pujala, Koyama lab/HHMI, 2016

headRange = [0 25];
tailRange = [75 100];
nFramesInTrl = 750;
fps = 500;
freqRange = [10 70];
dj = 1/4;
freqScale = 'lin';
noiseType = 'red';
stringency = 0;
plotOrNot = 1;
trlList = [];
xLim = [-50 650]; % Customized for escape
stimTime = 100; % In ms
cLim = [1 100];
onsetAlign = 1;
saveToProc = 1;
sigmaXY = [];
traceType = 'headTail';
diffOrNot = 0;

if nargin ==0 || isempty(procData)
    disp('Getting procData...')
    procData = OpenMatFile();
    [fPath,~] = fileparts(procData.Properties.Source);
    cd(fPath)
elseif nargin > 0
    try
        if isdir(procData)
            procData = OpenMatFile(procData,'nameMatchStr','proc');
        end
    catch
        a = 1;
    end
end

for jj = 1:numel(varargin)   
    if ischar(varargin{jj})      
        switch lower(varargin{jj})
            case 'hr'
                headRange =varargin{jj+1};
            case 'tr'
                tailRange = varargin{jj+1};
            case lower('nFramesInTrl')
                nFramesInTrl = varargin{jj+1};
            case 'fps'
                fps = varargin{jj+1};
            case 'dj'
                dj = varargin{jj+1};
            case 'freqrange'
                freqRange = varargin{jj+1};
            case 'freqscale'
                freqScale = varargin{jj+1};
            case 'noisetype'
                noiseType = varargin{jj+1};
            case 'stringency'
                stringency = varargin{jj+1};
            case 'sigmaxy'
                sigmaXY = varargin{jj+1};
            case 'plotornot'
                plotOrNot = varargin{jj+1};
            case 'trllist'
                trlList = varargin{jj+1};
            case 'xlim'
                xLim = varargin{jj+1};
            case 'clim'
                cLim = varargin{jj+1};
            case 'stimtime'
                stimTime = varargin{jj+1};
            case 'onsetalign'
                onsetAlign = varargin{jj+1};
            case 'savetoproc'
                saveToProc = varargin{jj+1};
            case 'tracetype'
                traceType = varargin{jj+1};
            case 'diffornot'
                diffOrNot = varargin{jj+1};
        end
    end
end

%% Compute head and tail segment orientations
disp('Reading tailCurv...')
tailCurv = procData.tailCurv;
nTrls = size(tailCurv,3)/nFramesInTrl;
if strcmpi(traceType, 'headTail') 
    disp('Getting head orientation...')
    or.head = GetSegOrientationFromTailCurv(tailCurv,headRange);
    disp('Getting tail orientation...')
    or.tail = GetSegOrientationFromTailCurv(tailCurv,tailRange);
    or.head_trl = reshape(or.head,nFramesInTrl,nTrls)';
    or.head_trl = or.head_trl - repmat(or.head_trl(:,1),1,size(or.head_trl,2)); % Zero 1st point in trl
    or.head = reshape(or.head_trl',1,numel(or.head_trl)); % To prevent outrageous std because of jumps at start of trls
    or.tail_trl = reshape(or.tail,nFramesInTrl,nTrls)';
    or.tail_trl = or.tail_trl - repmat(or.tail_trl(:,1),1,size(or.tail_trl,2));
    or.tail = reshape(or.tail_trl',1,numel(or.tail_trl));
    data.or = or;
elseif strcmpi(traceType, 'curv')
    disp('Getting whole body curvature...')
    curv = GetTailTangents(tailCurv,5);
    curv = curv(end,:);
    curv_trl = reshape(curv,nFramesInTrl,nTrls)';
    data.curv = curv;
    data.curv_trl = curv_trl;  
elseif strcmpi(traceType,'both')    
    disp('Getting head, tail orientations and body curvature...')
    or.head = GetSegOrientationFromTailCurv(tailCurv,headRange);
    or.tail = GetSegOrientationFromTailCurv(tailCurv,tailRange);
    or.head_trl = reshape(or.head,nFramesInTrl,nTrls)';
    or.head_trl = or.head_trl - repmat(or.head_trl(:,1),1,size(or.head_trl,2)); % Zero 1st point in trl
    or.head = reshape(or.head_trl',1,numel(or.head_trl)); % To prevent outrageous std because of jumps at start of trls
    or.tail_trl = reshape(or.tail,nFramesInTrl,nTrls)';
    or.tail_trl = or.tail_trl - repmat(or.tail_trl(:,1),1,size(or.tail_trl,2));
    or.tail = reshape(or.tail_trl',1,numel(or.tail_trl));
    data.or = or; 
    curv = GetTailTangents(tailCurv,5);
    curv = curv(end,:);
    curv_trl = reshape(curv,nFramesInTrl,nTrls)';
    data.curv = curv;
    data.curv_trl = curv_trl;
 end

time_trl = (0:nFramesInTrl-1)*1000/fps; % In ms

if onsetAlign
    onsets = nan(nTrls,1);
    blah = procData.elicitedSwimInfo;
    blah = blah.onset;
    nanInds = [];
    for trl = 1:nTrls        
        if ~isempty(blah{trl})
            onsets(trl) = blah{trl}(1);
        else
            onsets(trl)= nan;
            nanInds = [nanInds,trl];
        end
    end
    onsets = onsets + stimTime;
else
    onsets = repmat(stimTime,nTrls,1);
end

data.time = time_trl;
if isempty(trlList)
    trlList = 1:nTrls;
end
onsets = onsets(trlList);
nanInds = find(isnan(onsets));
onsets(nanInds) = [];
trlList(nanInds) = [];

data.trlList = trlList;
data.onsets = onsets;
data.fps = fps;
data.freqRange = freqRange;
data.dj = dj;
data.stringency = stringency;
data.freqScale = freqScale;
data.xLim = xLim;
data.sigmaXY = sigmaXY;
data.cLim = cLim;
data.noiseType = noiseType;
data.traceType = traceType;
data.diffOrNot = diffOrNot;


%% Computing wavelet transforms
W = GetWTs(data);

%% Plotting wavelet transforms

if plotOrNot   
    PlotWTs(W)
end

if saveToProc
    tic
    disp('Appending wavelet data to procData.mat...')
    procData.Properties.Writable = true;
    procData.W = W;
    toc
end

varargout{1} = W;

end

function W = GetWTs(data)
traceType = data.traceType;

if strcmpi(traceType,'headTail')
    W = GetWTs_headTail(data);
elseif strcmpi(traceType,'curv')
    W = GetWTs_curv(data);
elseif strcmpi(traceType,'both')
    W_headTail = GetWTs_headTail(data);
    W_curv = GetWTs_curv(data);
    W = catstruct(W_headTail,W_curv);
end

function W = GetWTs_headTail(data)
trlList = data.trlList;
nTrls = numel(trlList);
xLim = data.xLim;
time = data.time;
noiseType = data.noiseType;
fps = data.fps;
sigmaXY = data.sigmaXY;
tLen_exp = (abs(diff(xLim))/1000)*fps;
freqRange = data.freqRange;
onsets = data.onsets;

disp('Computing wavelet transforms...')
W.head.coeff = cell(nTrls,1);
W.tail.coeff = W.head.coeff;
W.head.ts = W.head.coeff;
W.tail.ts = W.head.coeff;
W.head.avg = [];
W.tail.avg = [];
W.time = [];
if isempty(sigmaXY)
    sigma.hh = std(chebfilt(data.or.head,1/fps,freqRange));
    sigma.tt =  std(chebfilt(data.or.tail,1/fps,freqRange));
    sigma.ht = sigma.hh*sigma.tt;    
else
    sigma.ht = sigmaXY;
end
shortTrls = [];
count = 0;
responseCount = 0;
firstNonZeroFlag = 1;
t = [];
for trl = trlList(:)'
    count = count + 1;
    onInd = find(time>= onsets(count),1,'first');
    time_align = time - time(onInd);
    tInds = find(time_align >= xLim(1) & time_align <= xLim(2));
    if numel(tInds)< tLen_exp
        shortTrls = [shortTrls,trl];
    else
        t  = time_align(tInds);
        x = chebfilt(data.or.head_trl(trl,:),1/fps,freqRange);
        x = x(tInds);
        y  = chebfilt(data.or.tail_trl(trl,:),1/fps,freqRange);
        y = y(tInds);
        [W.head.ts{count},W.tail.ts{count}] = deal(x,y);      
        [W.head.coeff{count},freq, coi] = ComputeXWT(x(:),x(:),t(:)/1000,'freqRange',freqRange,'dj',data.dj,'stringency',data.stringency,...
            'sigmaXY',sigma.ht,'freqScale',data.freqScale,'noiseType',noiseType);     
        [W.tail.coeff{count},~] = ComputeXWT(y(:),y(:),t(:)/1000,'freqRange',freqRange,'dj',data.dj,'stringency',data.stringency,...
            'sigmaXY',sigma.ht,'freqScale',data.freqScale,'noiseType',noiseType);
        if  (~isempty(W.head.coeff{count}) || ~isempty(W.tail.coeff{count})) && firstNonZeroFlag
            W.head.avg = abs(W.head.coeff{count});
            W.tail.avg = abs(W.tail.coeff{count});
            firstNonZeroFlag = 0;
            responseCount = 1;
        elseif ~isempty(W.head.coeff{count}) || ~isempty(W.tail.coeff{count})           
            if any(size(W.head.avg) ~= size(W.head.coeff{count}));
                shortTrls = [shortTrls,trl];
            else
                responseCount = responseCount  + 1;
                W.head.avg = W.head.avg + abs(W.head.coeff{count});
                W.tail.avg = W.tail.avg + abs(W.tail.coeff{count});
            end            
        end
    end
end
W.head.avg = W.head.avg/responseCount;
W.tail.avg = W.tail.avg/responseCount;
[~,delInds] = intersect(trlList,shortTrls);
trlList = setdiff(trlList,shortTrls);
W.head.coeff(delInds) = [];
W.head.ts(delInds) = [];
W.tail.coeff(delInds) = [];
W.tail.ts(delInds) = [];
W_all_head = reshape([W.head.coeff{:}],[size(W.head.avg),numel(trlList)]);
W.head.std = std(W_all_head,[],3);
W_all_tail = reshape([W.tail.coeff{:}],[size(W.tail.avg),numel(trlList)]);
W.tail.std = std(W_all_tail,[],3);
C.head = nan(numel(trlList),1);
C.tail = C.head;
for trl = 1:numel(trlList)  
    C.head(trl) = corr2(W.head.coeff{trl},W.head.avg);
    C.tail(trl) = corr2(W.tail.coeff{trl},W.tail.avg);  
end

if ~isempty(t)
    W.time = t;
else
    error('Not a single complete trial found, please check X limits!')
end
W.freqRange = data.freqRange;
W.stringency = data.stringency;
W.sigma = sigma;
W.freqScale = data.freqScale;
W.dj = data.dj;
W.freq = freq;
W.coi = 1./coi;
W.trlList = trlList;
W.shortTrls = shortTrls;
W.cLim = data.cLim;
W.head.corrVec = C.head;
W.tail.corrVec = C.tail;
W.traceType = data.traceType;
end

function W = GetWTs_curv(data)
trlList = data.trlList;
nTrls = numel(trlList);
xLim = data.xLim;
time = data.time;
noiseType = data.noiseType;
fps = data.fps;
sigmaXY = data.sigmaXY;
tLen_exp = (abs(diff(xLim))/1000)*fps;
freqRange = data.freqRange;
onsets = data.onsets;
diffOrNot = data.diffOrNot;

disp('Computing wavelet transforms...')
W.curv.coeff = cell(nTrls,1);
W.curv.ts = W.curv.coeff;
W.curv.avg = [];
W.time = [];
if isempty(sigmaXY)
    sigma.curv = std(chebfilt(data.curv,1/fps,freqRange))^2;
else
    sigma.curv = sigmaXY;
end
shortTrls = [];
count = 0;
responseCount = 0;
firstNonZeroFlag = 1;
t = [];
for trl = trlList(:)'
    count = count + 1;
    onInd = find(time>= onsets(count),1,'first');
    time_align = time - time(onInd);
    tInds = find(time_align >= xLim(1) & time_align <= xLim(2));
    if numel(tInds)< tLen_exp
        shortTrls = [shortTrls,trl];
    else
        t  = time_align(tInds);
        blah = data.curv_trl(trl,:);
%         dBlah = Standardize(gradient(blah))*2*max(blah);
        if diffOrNot
            dBlah = gradient(blah);
        else
           imf = MyEMD(blah,3);
           dBlah = imf(1).comp + imf(2).comp;
        end   
%         x = chebfilt(data.curv_trl(trl,:),1/fps,freqRange);
        x = blah;        
        x = x(tInds);       
        dBlah = dBlah(tInds);
        W.curv.ts{count} = x;
%         [W.curv.coeff{count},freq] = ComputeXWT(x(:),x(:),t(:)/1000,'freqRange',freqRange,'dj',data.dj,'stringency',data.stringency,...
%             'sigmaXY',sigma.curv,'freqScale',data.freqScale,'noiseType',noiseType);   
        [W.curv.coeff{count},freq,coi] = ComputeXWT(dBlah(:),dBlah(:),t(:)/1000,'freqRange',freqRange,...
            'dj',data.dj,'stringency',data.stringency,...
            'sigmaXY',sigma.curv,'freqScale',data.freqScale,'noiseType',noiseType);  
        if  ~isempty(W.curv.coeff{count}) && firstNonZeroFlag
            W.curv.avg = W.curv.coeff{count};            
            firstNonZeroFlag = 0;
            responseCount = 1;
        elseif ~isempty(W.curv.coeff{count})           
            if any(size(W.curv.avg) ~= size(W.curv.coeff{count}));
                shortTrls = [shortTrls,trl];
            else
                responseCount = responseCount  + 1;
                W.curv.avg = W.curv.avg + W.curv.coeff{count};                
            end            
        end
    end
end
W.curv.avg = W.curv.avg/responseCount;
[~,delInds] = intersect(trlList,shortTrls);
trlList = setdiff(trlList,shortTrls);
W.curv.coeff(delInds) = [];
W.curv.ts(delInds) = [];
W_all_curv = reshape([W.curv.coeff{:}],[size(W.curv.avg),numel(trlList)]);
W.curv.std = std(W_all_curv,[],3);
C.curv = nan(numel(trlList),1);
for trl = 1:numel(trlList)
    C.curv(trl) = corr2(W.curv.coeff{trl},W.curv.avg); 
end

if ~isempty(t)
    W.time = t;
else
    error('Not a single complete trial found, please check X limits!')
end
W.freqRange = data.freqRange;
W.stringency = data.stringency;
W.sigma = sigma;
W.freqScale = data.freqScale;
W.dj = data.dj;
W.freq = freq;
W.coi = 1./coi;
W.trlList = trlList;
W.shortTrls = shortTrls;
W.cLim = data.cLim;
W.curv.corrVec = C.curv;
W.traceType = data.traceType;
end

end



