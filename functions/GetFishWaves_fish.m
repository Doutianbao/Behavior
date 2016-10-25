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
% 'plotOrNot' - 0 or 1; the latter results in plotting
% 'trlList' - The list of trials to plot. If trList = [], then plots all
%   trials
% 'stimTime' - Stim onset time in milliseconds (default = 100);
% 'xLim' - X limits of plots (default = [-50 650], in milliseconds, customized for escape responses)
% 'cLim' - Color limits for normalized wavelet plots (default = [0.1 3]).
% 'onsetAlign' = 0 or 1. If 1, aligns timeseries of different trials based
%   onset of response (stored in procData.elicitedSwimInfo). If 0, aligns
%   w.r.t stimulus frame
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

if nargin ==0 || isempty(procData)
    disp('Getting procData...')
    procData = OpenMatFile();
    [fPath,~] = fileparts(procData.Properties.Source);
    cd(fPath)
elseif nargin > 0
    try
        if isdir(procData)
            procData = OpenMatFile(procData);
        end
    catch
        a = 1;
    end
end

for jj = 1:numel(varargin)
    if ischar(varargin{jj})
        switch lower(varargin{jj})
            case 'hr'
                headRange = varargin{jj+1};
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
        end
    end
end

%% Compute head and tail segment orientations
disp('Reading tailCurv...')
tailCurv = procData.tailCurv;
disp('Getting head orientation...')
or.head = GetSegOrientationFromTailCurv(tailCurv,headRange);
disp('Getting tail orientation')
or.tail = GetSegOrientationFromTailCurv(tailCurv,tailRange);
nTrls = size(tailCurv,3)/nFramesInTrl;
or.head_trl = reshape(or.head,nFramesInTrl,nTrls)';
or.head_trl = or.head_trl - repmat(or.head_trl(:,1),1,size(or.head_trl,2)); % Zero 1st point in trl
or.head = reshape(or.head_trl',1,numel(or.head_trl)); % To prevent outrageous std because of jumps at start of trls
or.tail_trl = reshape(or.tail,nFramesInTrl,nTrls)';
or.tail_trl = or.tail_trl - repmat(or.tail_trl(:,1),1,size(or.tail_trl,2));
or.tail = reshape(or.tail_trl',1,numel(or.tail_trl));
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

data.or = or;
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



%% Computing wavelet transforms
W = GetWTs(data);

%% Plotting wavelet transforms
if plotOrNot
    W.cLim = cLim;
    PlotWTs(W,data)
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
trlList = data.trlList;
nTrls = numel(trlList);
xLim = data.xLim;
time = data.time;
fps = data.fps;
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
sigma.hh = std(chebfilt(data.or.head,1/fps,freqRange));
sigma.tt =  std(chebfilt(data.or.tail,1/fps,freqRange));
sigma.ht = sigma.hh*sigma.tt;
shortTrls = [];
count = 0;
for trl = trlList(:)'
    count = count + 1;
    onInd = find(time>= onsets(count),1,'first');
    time_align = time - time(onInd);
    tInds = find(time_align >= xLim(1) & time_align <= xLim(2));
    if numel(tInds)< tLen_exp
        count = count-1;
        shortTrls = [shortTrls,trl];
    else
        t  = time_align(tInds);
        x = chebfilt(data.or.head_trl(trl,:),1/fps,freqRange);
        x = x(tInds);
        y  = chebfilt(data.or.tail_trl(trl,:),1/fps,freqRange);
        y = y(tInds);
        [W.head.ts{trl},W.tail.ts{trl}] = deal(x,y);
        [W.head.coeff{trl},freq] = ComputeXWT(x(:),x(:),t(:)/1000,'freqRange',freqRange,'dj',data.dj,'stringency',data.stringency,...
            'sigmaXY',sigma.ht,'freqScale',data.freqScale);
        [W.tail.coeff{trl},freq] = ComputeXWT(y(:),y(:),t(:)/1000,'freqRange',freqRange,'dj',data.dj,'stringency',data.stringency,...
            'sigmaXY',sigma.ht,'freqScale',data.freqScale);
        if count ==1
            W.head.avg = W.head.coeff{count};
            W.tail.avg = W.tail.coeff{count};
        else           
            W.head.avg = W.head.avg + W.head.coeff{trl};
            W.tail.avg = W.tail.avg + W.tail.coeff{trl};          
        end
    end
end
[~,delInds] = intersect(trlList,shortTrls);
trlList = setdiff(trlList,shortTrls);
W.head.coeff(delInds) = [];
W.head.ts(delInds) = [];
W.tail.coeff(delInds) = [];
W.tail.ts(delInds) = [];
W.head.avg = W.head.avg/count;
W.tail.avg = W.tail.avg/count;
if ~isempty(t)
    W.time = t;
else
    error('Not a single complete trial found, please check X limits!')
end
W.freqRange = data.freqRange;
W.stringency = data.stringency;
W.freqScale = data.freqScale;
W.dj = data.dj;
W.freq = freq;
W.trlList = trlList;
W.shortTrls = shortTrls;
end

function PlotWTs(W,data)
yShift = 300;
trlList = W.trlList;
t = W.time;
freq = W.freq;
cLim = W.cLim;
ax = {};
ax{1} = [1 0.39 0 0.61];
ax{2} = [1 0.39 0 0.21];
ax{3} = [1 0.2 0 0];
count = 0;
for trl = trlList(:)'
    count = count + 1;
    fh = figure('Name','Wavelet transforms of orientation timeseries');
    axH = CreateSubaxes(fh,ax{1},ax{2},ax{3});
    % -- Head wavelet --
    axes(axH(1));
    imagesc(t,freq,abs(W.head.coeff{trl}))
    set(gca,'ydir','normal','xtick',[],'clim',cLim);
    ylabel({'Head';' Freq (Hz)'})
    xlim([t(1) t(end)])
    title(['Head and tail orientation, trl = ' num2str(trl)])
    box off
    
    % -- Tail wavelet --
    axes(axH(2));
    imagesc(t, freq,abs(W.tail.coeff{trl}))
    set(gca,'ydir','normal','xtick',[],'clim',cLim);
    ylabel({'Tail' ; 'Freq (Hz)'})
    xlim([t(1) t(end)])
    
    % -- Head and tail orientation timeseries
    axes(axH(3))
    try
        plot(t,W.head.ts{trl},'g.')
    catch
        a = 1;
    end
    hold on
    plot(t,W.tail.ts{trl},'m.')
    xlim([t(1) t(end)])
    ylim([-200 200])
    box off
    ylabel({'Orientation'; '(deg)'})
    xTick = get(gca,'xtick');
    xTick(mod(xTick,100)~=0)=[];
    set(gca,'tickdir','out','xtick',xTick,'ytick',[-100 0 100],'color','k')
    xlabel('Time (ms)')
    shg
    linkaxes(axH,'x');
end
% -- Avg head and tail --
fh = figure('Name','Avg WT for head and tail');
ax{1} =[0.8 0.4 0 0.6];
ax{2} = [0.2 0.4 0.8 0.6];
ax{3} = [0.8 0.4 0 0.2];
ax{4} = [0.2 0.4 0.8 0.2];
ax{5} = [0.8 0.2 0 0];
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3},ax{4},ax{5});
axes(axH(1))
imagesc(t, freq,abs(W.head.avg))
set(gca,'ydir','normal','xtick',[],'clim',[cLim(1) cLim(2)*0.9]);
ylabel({'Head' ; 'Freq (Hz)'})
xlim([t(1) t(end)])
title('Avg WT for head and tail orientation')

axes(axH(2))
plot(mean(abs(W.head.avg),2),freq,'g')
ylim([freq(end) freq(1)])
box off
xlim([-inf inf])
set(gca,'ytick',[],'xaxislocation','top','color','k')

axes(axH(3))
imagesc(t, freq, abs(W.tail.avg))
set(gca,'ydir','normal','xtick',[],'clim',cLim);
ylabel({'Tail' ; 'Freq (Hz)'})
xlim([t(1) t(end)])

axes(axH(4))
plot(mean(abs(W.tail.avg),2),freq,'m')
box off
ylim([freq(end) freq(1)])
xlim([-inf inf])
set(gca,'ytick',[],'color','k')
xlabel('$\Sigma$ power','interpreter','latex')

axes(axH(5))
y = Standardize(mean(abs(W.tail.avg),1)) - Standardize(mean(abs(W.head.avg),1));
plot(t,y,'r')
hold on
plot(t, zeros(size(t)),'y--')
box off
set(gca,'tickdir','out','color','k','xtick',xTick)
xlim([t(1) t(end)])
ylim([-inf inf])
xlabel('Time(ms)')
ylabel('$\Sigma |T| - \Sigma |H|$','interpreter','latex')
end
