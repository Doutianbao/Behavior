function PlotWTs(W)
trlList = W.trlList;
t = W.time;
freq = W.freq;
cLim = W.cLim;
cLim_avg = [cLim(1) cLim(2)*0.7];
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
set(gca,'ydir','normal','xtick',[],'clim',cLim_avg);
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
set(gca,'ydir','normal','xtick',[],'clim',cLim_avg);
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