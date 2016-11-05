function varargout = GetFishWaves_group(pathList,varargin)
%GetFishWaves_fish Given procData.mat (created by FishSwim) or path to it
% returns wavelet transforms for head and tail orienations for all trials
% in addition to plotting these if specified
% W = GetFishWaves_group([]);
% W = GetFishWaves_group(pathList);
% W = GetFishWaves_group(pathToXlSheetWithPaths);
% W = GetFishWaves_group(.., 'hr',headRange,'tr',tailRange,'nFramesInTrl',nFramesInTrl,'fps',fps,...
%       'freqRange',freqRange,'dj',dj,'noiseType',noiseType,'stringency',stringency,'freqScale','lin',...
%       'plotOrNot', plotOrNot,'trlList',trlList,'xLim',xLim,'cLim',cLim);
% [W,procData] = GetFishWaves_group(...);
% Inputs:
% pathList - Cell array of directories where procData.mat files created by
%   FishSwim are created
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
% 'sigmaXY' - Signal cross variance to normalize WT coefficients by. If
%   sigmaXY =[], then auto determines for each fish, if sigmaXY = NaN, then
%   determines sigma for the entire group and normalizes individual values
%   by this factor.
% 'plotOrNot' - 0 or 1; the latter results in plotting
% 'trlList' - The list of trials to plot. If trList = [], then plots all
%   trials
% 'stimTime' - Stim onset time in milliseconds (default = 100);
% 'xLim' - X limits of plots (default = [-50 650], in milliseconds, customized for escape responses)
% 'cLim' - Color limits for normalized wavelet plots (default = [0.1 3]).
% 'onsetAlign' = 0 or 1. If 1, aligns timeseries of different trials based
%   onset of response (stored in procData.elicitedSwimInfo). If 0, aligns
%   w.r.t stimulus frame
% Outputs:
% W - Cell array containig wavelet transform data for each fish
% procData - Cell array containing pointers to procData.mat files
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
sigmaXY = [];
plotOrNot = 1;
trlList = [];
xLim = [-50 650]; % Customized for escape
stimTime = 100; % In ms
cLim = [1 100];
onsetAlign = 1;
saveToProc = 1;

currPath = cd;
if nargin ==0 || isempty(pathList)
    disp('Getting paths...')
    cd('S:\Avinash\Ablations and behavior')
    pathList = GetFilteredPathsFromXLS();
 elseif nargin > 0
    try
        if isdir(pathList)
            pathList = GetFilteredPathsFromXLS(pathList);
        else
            pathList = varargin{1};
        end
    catch
        a = 1;
    end
end
if ~iscell(pathList)
    pathList = {pathList};
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
        end
    end
end

nPaths = length(pathList);
grpFlag = 0;
sigmaXY_grp = nan(nPaths,1);
W = cell(nPaths,1);
procData = cell(nPaths,1);
for pp = 1:nPaths
    cd(pathList{pp})
    disp(pathList{pp})
    procData{pp} = OpenMatFile(pathList{pp});
    if isnan(sigmaXY) % Use group sigmaXY
        plotOrNot_new = 0;
        sigmaXY_new = [];
        grpFlag = 1;
        saveToProc_new = 0;
    else
        plotOrNot_new = plotOrNot;
    end
    W{pp}  = GetFishWaves_fish(procData{pp},'hr',headRange,'tr',tailRange,...
        'nFramesInTrl',nFramesInTrl,'fps',fps,'freqRange',freqRange,'dj',dj,...
        'plotOrNot',plotOrNot_new,'xLim',xLim,'cLim',cLim,'stimTime',stimTime,...
        'onsetAlign',onsetAlign,'saveToProc',saveToProc_new,'noiseType',noiseType,...
        'stringency',stringency,'sigmaXY',sigmaXY_new);
    sigmaXY_grp(pp) = W{pp}.sigma.ht;
    disp(['Completed for ' pathList{pp}])
end

disp('Adjusting WT coeffs based on group sigmaXY and appending to procData...')
if grpFlag
    sigmaXY = mean(sigmaXY_grp);
    for pp = 1:nPaths
        disp(pathList{pp})
        for trl = 1:length(W{pp}.tail.coeff)
            mult = sigmaXY/W{pp}.sigma.ht;
            W{pp}.head.coeff{trl} = W{pp}.head.coeff{trl} * mult; 
            W{pp}.tail.coeff{trl} = W{pp}.tail.coeff{trl} * mult;
        end
        W{pp}.head.avg = W{pp}.head.avg * mult;
        W{pp}.tail.avg = W{pp}.tail.avg * mult;
        if saveToProc
            disp('Appending WT to procData...')
            procData{pp}.Properties.Writable = true;
            procData{pp}.W = W{pp};
        end
        if plotOrNot
            PlotWTs(W{pp})
        end
    end
end

varargout{1} = W;
varargout{2} = procData;

end

function PlotWTs(W)
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
    imagesc(t,freq,abs(W.head.coeff{count}))
    set(gca,'ydir','normal','xtick',[],'clim',cLim);
    ylabel({'Head';' Freq (Hz)'})
    xlim([t(1) t(end)])
    title(['Head and tail orientation, trl = ' num2str(trl)])
    box off
    
    % -- Tail wavelet --
    axes(axH(2));
    imagesc(t, freq,abs(W.tail.coeff{count}))
    set(gca,'ydir','normal','xtick',[],'clim',cLim);
    ylabel({'Tail' ; 'Freq (Hz)'})
    xlim([t(1) t(end)])
    
    % -- Head and tail orientation timeseries
    axes(axH(3))
    plot(t,W.head.ts{count},'g.')
    hold on
    plot(t,W.tail.ts{count},'m.')
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

