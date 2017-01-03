function varargout = AnalyzeFreeSwims_nCycles(varargin)
%AnalyzeFreeSwims Given procData.mat generated by SlowSwim, returns some
%   useful info after analyzing swims
% out = AnalyzeFreeSwims_nCycles();
% out = AnalyzeFreeSwims_nCycles(procData,'fps',fps,'nFramesInTrl',nFramesInTrl,'paramList', paramList);
% out = AnalyzeFreeSwim_nCycles([],'fps',fps,...);
% Inputs:
% procData - Processed data .mat file created by SlowSwim.m and/or related
%   scripts
% 'fps'  - Frames per second (default: 500)
% 'nFramesInTrl' - Number frames in a single trial (default: 750)
% 'paramList' - List of params to extract in this function (default:
%   {'bodyAmp', 'angVel', 'headAmp'}).
%   'bodyAmp' - Total body bending amplitudes
%   'angVel' - Period for total body bends
%   'headAmp' - Amplitudes for only head segment bends
% 'freqRange' - Frequency range for crosswavelet transform
% Avinash Pujala, Koyama lab/HHMI, 2016

fps = 500;
nFramesInTrl = 750;
preStimPeriod = 0.1;
xLim = [0 0.8*1000];  % For vibration
% xLim = [0 1.5*1000]; % For dark flash
stringency = 1.5;
paramList_all = {'bodyAmp','angVel','headAmp'};
paramList = paramList_all;
freqRange = [15 60];
freqChannels = linspace(freqRange(1),freqRange(2),4);

% cd('S:\Avinash\Ablations and behavior\Intermediate RS\20160715')
if nargin ==0
    procData = OpenMatFile();
elseif nargin == 1
    procData = varargin{1};
    if isempty(procData)
        procData = OpenMatFile();
    elseif isdir(procData)
        procData = OpenMatFile(procData,'nameMatchStr','proc');
    end
    fNames = fieldnames(procData);
    if  any(strcmpi(fNames,'fps'))
        fps = procData.fps;
    end
    if any(strcmpi(fNames,'nFramesInTrl'))
        nFramesInTrl = procData.nFramesInTrl;
    end
else
    procData = varargin{1};
    if isempty(procData)
        procData = OpenMatFile();
    elseif isdir(procData)
        procData = OpenMatFile(procData,'nameMatchStr','proc');
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
                case lower('freqRange')
                    freqRange = varargin{jj+1};
            end
        end
    end
end

% cd(pathName)
if ~iscell(paramList)
    paramList ={paramList};
end

paramInds = zeros(length(paramList),1);
for jj = 1:length(paramList)
    ind = find(strcmpi(paramList_all,paramList{jj}));
    if ~isempty(ind)
        paramInds(jj) = ind;
    end
end

disp('Reading data from procData...')
tic
disp('Tail curvature...')
tailCurv = procData.tailCurv;
procData.Properties.Writable = true;
toc

time = (0:size(tailCurv,3)-1)*(1/fps);
nTrls = length(time)/nFramesInTrl;

tA_5 = GetTailTangents(tailCurv,5);
curv = tA_5(end,:)';
curv_seg1 = tA_5(1,:)';
% curv_head = tA_5(1,:)';
disp('Getting head orientation from body curve..')
curv_head = GetSegOrientationFromTailCurv(tailCurv);
curv_head = chebfilt(curv_head,1/fps,50,'low');
tA_trl = reshape(curv,nFramesInTrl,nTrls);
tA_trl_seg1 = reshape(curv_seg1,nFramesInTrl,nTrls);
tA_trl_head = reshape(curv_head,nFramesInTrl,nTrls);
tA_trl_head = tA_trl_head - repmat(tA_trl_head(1,:),size(tA_trl_head,1),1);
curv_head = reshape(tA_trl_head,1,prod(size(tA_trl_head)));
time_trl = time(1:nFramesInTrl);
pkThr1 = stringency*std(tA_trl(:));
pkThr3 = stringency*std(tA_trl_head(:));
maxFreq = 80;
minIntPts = max(floor((0.5/maxFreq)*fps),1);

blah = chebfilt(tA_trl,1/fps,30,'low');
dTrace_all = gradient(blah')';
pkThr2 = stringency*std(dTrace_all(:));

disp('Getting peak info...')
out = struct;
figure('Name','Pk info')
out.bendAmp = cell(nTrls,1);
out.bendPer = out.bendAmp;
out.onset = out.bendAmp;
out.bendAngVel = out.bendAmp;
for trl = 1:nTrls
    tr = tA_trl(:,trl);
    tr_seg1 = tA_trl_seg1(:,trl);
    tr_head = tA_trl_head(:,trl);
    pks1 = GetPks(tr,'polarity',0, 'peakThr',pkThr1,'thrType','rel','minPkDist',minIntPts);
    mT = max(abs(tr));
    dTrace = dTrace_all(:,trl);
    pks2 = GetPks(dTrace,'polarity',0, 'peakThr',pkThr2,'thrType','rel','minPkDist',minIntPts);
    try
    pks3 = GetPks(tr_head,'polarity',0, 'peakThr',20,'thrType','rel','minPkDist',minIntPts);
    catch
        pks3 = GetPks(tr_head,'polarity',0, 'peakThr',20,'thrType','rel','minPkDist',minIntPts);
    end
    mDT = max(dTrace);
    sf = 0.5*mT/mDT;
    x_trace1 = 0;
    paramInds(paramInds==0)=[];
    for traceType = paramInds(:)'
        cla
        maxY = max(tA_trl(:,trl));
        minY = min(tA_trl(:,trl));
        if traceType ==1
            plot(time_trl*1000,tr)
            hold on
            plot(time_trl(pks1)*1000,tr(pks1),'ko')
            plot(time_trl*1000,dTrace*sf,'r:')
            yl = [min([minY,-250]) max([maxY,250])];
            ylim(yl)
            plot([preStimPeriod, preStimPeriod]*1000,yl,'--','color',[0 0.5 0])
        elseif traceType ==2
            plot(time_trl*1000,tr,'b:')
            hold on
            plot(time_trl*1000,dTrace*sf,'r')
            plot(time_trl(pks2)*1000,dTrace(pks2)*sf,'ko')
            yl = [min([minY,-250]) max([maxY,250])];
            ylim(yl)
            plot([preStimPeriod, preStimPeriod]*1000,yl,'--','color',[0 0.5 0])
        else
            cla
            plot(time_trl*1000,tr_head)
            hold on
            plot(time_trl(pks3)*1000,tr_head(pks3),'ro')
            plot(time_trl*1000,tr/2 + tr_head(1),'m:')
            yl = [min(tr_head),max(tr_head)];
            ylim([-inf inf])
            plot([preStimPeriod preStimPeriod]*1000, yl,'--','color',[0 0.5 0])
        end
        ylabel(paramList_all{traceType})
        box off
        xlim(xLim)
        set(gca,'xtick',[100 500 1000 15000])
        title(['Click on 5 pts to get onset, 1st and 3rd undulation info, Trl # ' num2str(trl)])
        shg
        [x,y,~] = ginput_plot();
        if ~isempty(x) && traceType ==1
            x_trace1 = x;
            preInds = find((time_trl(pks1)*1000) < min(x));
            postInds = find((time_trl(pks1)*1000) > max(x));
            prePostInds = union(preInds,postInds);
            pks1(prePostInds)=[];
            time_trl = time_trl(:);
            tr = tr(:);
            x = [x(:);time_trl(pks1)*1000];
            y = [y(:); tr(pks1)];
            [x,inds] = sort(x);
            y = y(inds);
            dx = (diff(x)/1000)*fps;
            inds = find(dx < minIntPts);
            x(inds) = [];
            y(inds)=[];
            tr = chebfilt(tr,1/fps,10,'high');
            tr_seg1 = chebfilt(tr_seg1,1/fps,10,'high');
            [Wxy,freq] =  ComputeXWT(tr(:),tr_seg1(:),time_trl(:),freqRange,1/64,0,'all');
            [~, freq_pow] = GetWaveFreqTimeseries(Wxy,freq);
            inds = find((time_trl*1000)>=x(1) & (time_trl*1000) <= x(end));
            per = (1000./freq_pow.freq(inds));
            ph = freq_pow.phase(inds);
            %             w = freq_pow.pow(inds);  % Not using weighted linear fit at the moment, which I can do using lscov, if need be
            tau = time_trl(inds);
            try
                tau = (tau-tau(1));
            catch
                error('Must click at two locations at least')
            end
            B_per = polyfit(tau(:),per(:),1);
            B_ph = polyfit(tau(:),ph(:),1);
            onset = x(1)- (preStimPeriod*1000);
            for bend = 2:numel(x)
                out.onset{trl}(bend-1) = onset;
                out.bendAmp{trl}(bend-1) = y(bend)-y(bend-1);
                out.bendPer{trl}(bend-1) = x(bend)-x(bend-1);
                out.perSlope{trl}(bend-1) = B_per(1);
                out.perOff{trl}(bend-1) = B_per(2);
                out.phSlope{trl}(bend-1) = B_ph(1);
                out.phOff{trl}(bend-1) = B_ph(2);
                [~, ind] = min(abs((time_trl*1000)-x(bend)));
                out.xwPer{trl}(bend-1) = (1000/freq_pow.freq(ind));
                out.xWPhase{trl}(bend-1) = freq_pow.phase(ind);
                for freqChan = 1:numel(freqChannels)-1
                    fInds = find(freq > freqChannels(freqChan) & freq <= freqChannels(freqChan+1));
                    suffix = [num2str(freqChannels(freqChan)) '_' num2str(freqChannels(freqChan+1))];
                    fldName = ['xw_pow_' suffix];
                    out.(fldName){trl}(bend-1) = mean(abs(Wxy(fInds,ind)));
                    fldName = ['xw_phase_' suffix];
                    out.(fldName){trl}(bend-1) = angle(mean(Wxy(fInds,ind)))*180/pi;
                end
            end
        elseif ~isempty(x) && traceType ==2
            [~, inds] = sort(x);
            y = y(inds);
            preInds = find((time_trl(pks2)*1000) < min(x_trace1));
            postInds = find((time_trl(pks2)*1000) > max(x_trace1));
            prePostInds = union(preInds,postInds);
            pks2(prePostInds)=[];
            time_trl = time_trl(:);
            dTrace = dTrace(:);
            x = [x(:);time_trl(pks2)*1000];
            y = [y(:); dTrace(pks2)];
            [x,inds] = sort(x);
            y = y(inds);
            dx = (diff(x)/1000)*fps;
            inds = find(dx < minIntPts);
            x(inds) = [];
            y(inds)=[];
            for bend = 1:numel(x)
                out.bendAngVel{trl}(bend) = y(bend)/sf;
            end
        elseif ~isempty(x) && traceType ==3
            x_trace1 = x;
            preInds = find((time_trl(pks3)*1000) < min(x));
            postInds = find((time_trl(pks3)*1000) > max(x));
            prePostInds = union(preInds,postInds);
            pks3(prePostInds)=[];
            time_trl = time_trl(:);
            tr = tr(:);
            x = [x(:);time_trl(pks3)*1000];
            y = [y(:); tr_head(pks3)];
            [x,inds] = sort(x);
            y = y(inds);
            dx = (diff(x)/1000)*fps;
            inds = find(dx < minIntPts);
            x(inds) = [];
            y(inds)=[];
            
            tr_head = chebfilt(tr_head,1/fps,10,'high');
            [Wxy,freq] =  ComputeXWT(tr_head(:),tr_head(:),time_trl(:),freqRange,1/64,0,'all');
            [~, freq_pow] = GetWaveFreqTimeseries(Wxy,freq);
            inds = find((time_trl*1000)>=x(1) & (time_trl*1000) <= x(end));
            per = (1000./freq_pow.freq(inds));
            %             w = freq_pow.pow(inds);  % Not using weighted linear fit at the moment, which I can do using lscov, if need be
            tau = time_trl(inds);
            tau = (tau-tau(1));
            B_per = polyfit(tau(:),per(:),1);
            onset_head = x(1)- (preStimPeriod*1000);
            for bend = 2:numel(x)               
                out.headAmp{trl}(bend-1) = y(bend)-y(bend-1);
                out.headPer{trl}(bend-1) = x(bend)-x(bend-1);
                out.headPerSlope{trl}(bend-1) = B_per(1);
                out.headPerOff{trl}(bend-1) = B_per(2);
                for freqChan = 1:numel(freqChannels)-1
                    [~, ind] = min(abs((time_trl*1000)-x(bend)));
                    fInds = size(Wxy,1)-find(freq > freqChannels(freqChan) & freq <= freqChannels(freqChan+1));
                    suffix = [num2str(freqChannels(freqChan)) '_' num2str(freqChannels(freqChan+1))];
                    fldName = ['head_xw_pow_' suffix];
                    out.(fldName){trl}(bend-1) = mean(abs(Wxy(fInds,ind)));
                end
            end
        else
            for bend = 1:numel(x)
                out.bendAngVel{trl}(bend) = y(bend)/sf;
            end
        end
    end
    
end

count = 0;
repeat = true;
while count < 5 && repeat
    saveOrNot = input('Append pk data to procData? (y/n)','s');
    count = count + 1;
    if strcmpi(saveOrNot,'y')
        procData.Properties.Writable = true;
        procData.hOr = curv_head;
        if sum(strcmpi(fieldnames(procData),'elicitedSwimInfo'))
            esi = procData.elicitedSwimInfo;
            esi = catstruct(esi,out);
        end
        procData.elicitedSwimInfo = esi;
        repeat = false;
    elseif strcmpi(saveOrNot,'n')
        disp('Data not appended to procData!')
        repeat = false;
    else
        disp('Please enter "y" or "n"')
        repeat = true;
    end
end

varargout{1} = out;
varargout{2} = procData;
end

