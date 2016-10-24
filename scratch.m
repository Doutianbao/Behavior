
trl =1;
nFramesInTrl = 750;
fps = 500;
freqRange = [20 50];
dj = 1/64;
time = (0:size(tA_5,2)-1)*(1/fps);
inds = (trl-1)*nFramesInTrl + 1 : (trl-1)*nFramesInTrl + 500;
x = chebfilt(tA_5(end,inds),1/fps,10,'high');
y =  chebfilt(tA_5(1,inds),1/fps,10,'high');
t = time(inds);
% xx = chebfilt(hOr,1/fps,50,'low');
% xx  = xx(inds); 
[Wxy,freq] = ComputeXWT(x(:),y(:),t(:),freqRange,1/64,0,'all');
% [Wxy,~] = ComputeXWT(xx(:),xx(:),t(:),freqRange,1/64,0,'all');

figure
imagesc(t*1000, log2(freq),abs(Wxy))
ytl = 2.^str2num(get(gca,'yticklabel'));
set(gca,'yticklabel',ytl,'ydir','normal');
shg


[freq_mean, freq_pow] = GetWaveFreqTimeseries(Wxy,freq);

PlotWxy(Wxy, t(:), freq, [zscore(x(:)) zscore(y(:))])

figure
ax =[];
% thr  = mean(freq_pow.pow) + 0*std(freq_pow.pow);
thr = max(freq_pow.pow)*0.1;
inds  = find(freq_pow.pow >= thr);
ax(1) = subplot(3,1,1);
plot(t*1000,[zscore(x(:)), zscore(y(:))])
ylabel('Zscore amp')

ax(2) = subplot(3,1,2);
plot(t(inds)*1000,1000./freq_mean.freq(inds),'.')
hold on
plot(t(inds)*1000,1000./freq_pow.freq(inds),'r.')
ylabel('Period (ms)')
ylim([10 70])


ax(3) = subplot(3,1,3);
plot(t(inds)*1000,abs(freq_mean.phase(inds)),'.');
hold on
plot(t(inds)*1000,abs(freq_pow.phase(inds)),'r.');
ylabel('Phase shift (deg)')
ylim([0 180])

% ax(4) = subplot(4,1,4);
% var = (abs(freq_mean.phase(inds))/360).*(1000./freq_mean.freq(inds));
% plot(t(inds)*1000,var,'.');
% ylabel('Phase shift (ms)')
% ylim([0 inf])

linkaxes(ax,'x')



%%
% trl  = 1;
figure('Name', 'Dynamic pow in freq channels')
for trl = 1:20
blah = [out.xw_pow_15_30{trl}; out.xw_pow_30_45{trl}; out.xw_pow_45_60{trl}];  
blah = blah./repmat(max(blah,[],2),1,size(blah,2));
plot(blah(1,:),'.-')
hold on
% plot(blah(2,:),'.-','color',[0 0.4 0])
% plot(blah(3,:),'r.-')
box off
xlim([-inf inf])
end
legend('15-30 Hz', '30-45 Hz', '45-60 Hz' )





