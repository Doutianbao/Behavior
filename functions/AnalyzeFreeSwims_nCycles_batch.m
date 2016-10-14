function out = AnalyzeFreeSwims_nCycles_batch(pathList,varargin)
%AnalyzeFreeSwims_nCycles_bath Given a lit of paths that point to
%   procData.mat files created by SlowSwim and related scripts, appends
%   elicited swim info to procData and also returns a cell array containing
%   all this info.
% out = AnalyzeFreeSwims_nCycles_batch(pathList);
% out = AnalyzeFreeSwims_ncycles_batch(procData,'fps',fps,'nFramesInTrl',nFramesInTrl,'paramList', paramList,'xLim',xLim);
% Inputs:
% pathList - A list (cell array) of paths to procData.m
%   scripts
% 'fps'  - Frames per second (default: 500)
% 'nFramesInTrl' - Number frames in a single trial (default: 750)
% 'paramList' - List of params to extract in this function (default:
%   {'bodyAmp', 'bodyPer', 'headAmp'}).
%   'bodyAmp' - Total body bending amplitudes
%   'angVel' - Period for total body bends
%   'headAmp' - Amplitudes for only head segment bends
% 'xLim' - Limits of x-axes for interative selection of peaks
% 
% Avinash Pujala, Koyama lab/HHMI, 2016

fps = 500;
nFramesInTrl = 750;
preStimPeriod = 0.1;
% xLim = [0 8-00];  % For vibration
xLim = [0 1500]; % For dark flash
stringency = 1.5;
paramList_all = {'bodyAmp','angVel','headAmp'};
paramList = paramList_all;

if isempty(pathList)
    error('Must input a list of paths')
end

if ~iscell(pathList)
    pathList = {pathList};
end

for jj = 1:numel(varargin)
    if ischar(varargin{jj})
        switch lower(varargin{jj})
            case 'fps'
                fps = varargin{jj+1};
            case lower('nFramesInTrl')
                nFramesInTrl = varargin{jj+1};
            case lower('preStimPeriod')
                preStimPeriod = varargin{jj+1};
            case lower('paramList')
                paramList = varargin{jj+1};
            case lower('xLim')
                xLim = varargin{jj+1};
        end
    end
end


out = cell(length(pathList),1);
for pp = 1:length(pathList)
    cp = pathList{pp};
    [~,lastDir] = fileparts(cp);
    if ~strcmpi(lastDir,'proc')
        cp = fullfile(cp,'proc');
    end
    cd(cp)
    out{jj} = AnalyzeFreeSwims_nCycles([],'fps',fps,'nFramesInTrl',nFramesInTrl, 'preStimPeriod', preStimPeriod,...
        'paramList',paramList,'xLim',xLim);
end

end

