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

W = cell(length(pathList),1);
for pp = 1:length(pathList)
    cd(pathList{pp})
    procData = OpenMatFile();
    W{pp}  = GetFishWaves_fish(procData,'hr',headRange,'tr',tailRange,...
        'nFramesInTrl',nFramesInTrl,'fps',fps,'freqRange',freqRange,'dj',dj,...
        'plotOrNot',plotOrNot,'xLim',xLim,'cLim',cLim,'stimTime',stimTime,...
        'onsetAlign',onsetAlign,'saveToProc',saveToProc,'noiseType',noiseType,...
        'stringency',stringency,'sigmaXY',sigmaXY); 
    disp(['Completed for ' pathList{pp}])
end

end


