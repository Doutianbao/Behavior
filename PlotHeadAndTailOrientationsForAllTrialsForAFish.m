
%% Inputs
segRange.head = [0 25];
segRange.tail = [75 100];
nFramesInTrl = 750;
fps = 500;
dj = 1/512;
freqRange = [8 64];
outPath = 'S:\Avinash\Ablations and behavior\PlayPen';

%% Get procData
tic
disp('Getting procData...')
procData = OpenMatFile();
[fPath,~] = fileparts(procData.Properties.Source);
cd(fPath)

%% Compute head and tail segment orientations
disp('Reading tailCurv...')
tailCurv = procData.tailCurv;
disp('Getting head orientation...')
or.head = GetSegOrientationFromTailCurv(tailCurv,segRange.head);
disp('Getting tail orientation')
or.tail = GetSegOrientationFromTailCurv(tailCurv,segRange.tail);
nTrls = size(tailCurv,3)/nFramesInTrl;
or.head_trl = reshape(or.head,nFramesInTrl,nTrls)';
or.head_trl = or.head_trl - repmat(or.head_trl(:,1),1,size(or.head_trl,2)); % Zero 1st point in trl
or.tail_trl = reshape(or.tail,nFramesInTrl,nTrls)';
or.tail_trl = or.tail_trl - repmat(or.tail_trl(:,1),1,size(or.tail_trl,2));
time_trl = (0:nFramesInTrl-1)*1000/fps; % In ms

data.or = or;
data.time = time_trl;




%% Plotting data
trlList = 1:nTrls;
% trlList = 1:15;
% trlList = [1 5 6 7 8 9 10 11 14 15]; % ctrl
% trlList = [1 4 6 7 9 10 17 18 19 20]; % abl
yShift = 300;
xLim = [50 650];
stimTime = 100; % In ms
% cLim = [0.01 0.4]; % Ctrl
cLim = [0.02 1]; % Abl

%### Head and tail orientation
yOff = zeros(length(trlList),1);
figure('Name','Head and tail orientations')
count = 0;
for trl = trlList(:)'
    count = count + 1;
    yOff(count) = (count-1)*yShift;
    plot(time_trl,or.head_trl(trl,:)-yOff(count),'g')
    hold on
    plot(time_trl,or.tail_trl(trl,:)-yOff(count),'m')    
end
yLim = [min(yOff)-yShift, max(yOff)+yShift];
yOff = sort(-yOff,'ascend');
ytl = sort(trlList,'descend');
lh = legend('Head','Tail');
set(lh,'color','k','textcolor',[0.9 0.9 0.9])
set(gca,'color','k','ytick',yOff,'yticklabel',ytl,'tickdir','out')
box off
plot([stimTime, stimTime],-yLim,'y--')
xlim(xLim)
xlabel('Time (ms)')
ylabel('Trl #')
title('Head and tail orientations for different trials')

%###### Wavelet plots for tail
W.head = cell(nTrls,1);
W.tail = W.head;
W.avg.head = [];
W.avg.tail = [];
ax = cell(2,1);
xPow = zeros(nTrls,1);
sigma.hh = std(or.head)^2;
sigma.tt = std(or.tail)^2;
count = 0;
for trl = trlList
    count = count + 1;
    fh = figure('Name','Wavelet transforms of orientation timeseries');
    tInds = find(time_trl >= xLim(1) & time_trl <= xLim(2));
    t  = time_trl(tInds);
    x = chebfilt(or.head_trl(trl,:),1/fps,freqRange);
    x = x(tInds);
    y  = chebfilt(or.tail_trl(trl,:),1/fps,freqRange);
    y = y(tInds);
    [W.head{trl},freq] = ComputeXWT(x(:),x(:),t(:)/1000,freqRange, dj,0,'all',sigma.hh);
    [W.tail{trl},freq] = ComputeXWT(y(:),y(:),t(:)/1000,freqRange, dj,0,'all',sigma.tt);
    if count ==1
        W.avg.head = W.head{trl};
        W.avg.tail = W.tail{trl};
    else
        W.avg.head = W.avg.head + W.head{trl};
        W.avg.tail = W.avg.tail + W.tail{trl};
    end
    %     blah = W.tail{trl};
    %     xPow(trl) = mean(abs(blah(:)));
    ax{1} = [1 0.39 0 0.61];
    ax{2} = [1 0.39 0 0.21];
    ax{3} = [1 0.2 0 0];
    axH = CreateSubaxes(fh,ax{1},ax{2},ax{3});
    
    % -- Head wavelet --
    axes(axH(1));
    imagesc(t, log2(freq),abs(W.head{trl}))
    ytl = str2num(get(gca,'yticklabel'));
    ytl(mod(ytl,1)~=0) =[];
    yTick = ytl;
    ytl = 2.^ytl;
    set(gca,'ytick',yTick,'yticklabel',ytl,'ydir','normal','xtick',[],'clim',cLim);
    ylabel({'Head';' Freq (Hz)'})
    xlim([t(1) t(end)])
    title(['Head and tail orientation, trl = ' num2str(trl)])
    
    % -- Tail wavelet --
    axes(axH(2));
    imagesc(t, log2(freq),abs(W.tail{trl}))
    set(gca,'ytick',yTick(1:end-1),'yticklabel',ytl(1:end-1),'ydir','normal','xtick',[],'clim',cLim);
    ylabel({'Tail' ; 'Freq (Hz)'})
    xlim([t(1) t(end)])
    
    % -- Head and tail orientation timeserie     
    axes(axH(3))
    plot(t,x,'g.')
    hold on
    plot(t,y,'m.')
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
W.avg.head = W.avg.head/count;
W.avg.tail = W.avg.tail/count;

% -- Avg head and tail --
fh = figure('Name','Avg WT for head and tail');
ax{1} =[0.8 0.4 0 0.6];
ax{2} = [0.2 0.4 0.8 0.6];
ax{3} = [0.8 0.4 0 0.2];
ax{4} = [0.2 0.4 0.8 0.2];
ax{5} = [0.8 0.2 0 0];
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3},ax{4},ax{5});

axes(axH(1))
imagesc(t, log2(freq),abs(W.avg.head))
set(gca,'ytick',yTick(1:end-1),'yticklabel',ytl(1:end-1),'ydir','normal','xtick',[],'clim',cLim);
ylabel({'Head' ; 'Freq (Hz)'})
xlim([t(1) t(end)])
title('Avg WT for head and tail orientation')

axes(axH(2))
plot(sum(abs(W.avg.head),2),log2(freq),'g')
box off
ylim(log2([freq(end) freq(1)]))
xlim([-inf inf])
set(gca,'ytick',[],'xaxislocation','top','color','k')

axes(axH(3))
imagesc(t, log2(freq),abs(W.avg.tail))
set(gca,'ytick',yTick(1:end-1),'yticklabel',ytl(1:end-1),'ydir','normal','xtick',[],'clim',cLim);
ylabel({'Tail' ; 'Freq (Hz)'})
xlim([t(1) t(end)])

axes(axH(4))
plot(sum(abs(W.avg.tail),2),log2(freq),'m')
box off
ylim(log2([freq(end) freq(1)]))
xlim([-inf inf])
set(gca,'ytick',[],'color','k')
xlabel('$\Sigma$ power','interpreter','latex')

axes(axH(5))
y = zscore(sum(abs(W.avg.tail) - abs(W.avg.head),1));
plot(t,y,'r')
hold on
plot(t, zeros(size(t)),'y--')
box off
set(gca,'tickdir','out','color','k','xtick',xTick)
xlim([t(1) t(end)])
ylim([-inf inf])
xlabel('Time(ms)')
ylabel('$\Sigma |T| - \Sigma |H|$','interpreter','latex')


break

%% Saving data
% saveOrNot = input('Save orientation data? (y/n)','s');
% if strcmpi(saveOrNot,'y')
%     dataMat = [time_trl(:),or.head_trl', or.tail_trl'];
%     disp('Saving data...')
%     fileName = fullfile(outPath,'HeadAndTailOrientation.mat');
%     save(fileName,'dataMat')
% end
% toc

