
% ---------------------------------------------------------
% Attempting a version of wavelet and xwt with linear scales, never mind
% the redundant computations that result from this
% ---------------------------------------------------------

%% Creating a signal
freqRange = [10 80];
dj = 1/64;
freqScale = 'log';
% dj = 1/2;
% freqScale = 'lin';
f0 = 15;
f1 = 55;


t = linspace(0,2,1000);
dt = t(2)-t(1);
y = chirp(t,f0,t(end),f1,'q',[],'concave');
y = cat(1, y(:),flipud(y(:)));
y = y + randn(size(y))/6;
t = linspace(0,2*t(end),length(y));

[Wxy,freq] = ComputeXWT(y,y,t,'dj',dj, 'freqRange',freqRange,'stringency',0,'freqScale',freqScale);


ax = [];
figure('Name','WT')
ax(1) = subplot(3,1,1);
if strcmpi(freqScale,'lin')
    imagesc(t,freq,abs(Wxy))
    set(gca,'tickdir','out','ydir','normal')
else
    imagesc(t,log2(freq),abs(Wxy))
    yt = get(gca,'ytick');
    yt(mod(yt,1)~=0) = [];
    ytl = 2.^yt;
    set(gca,'tickdir','out','ydir','normal','ytick',yt,'yticklabel',ytl)
end
    box off


ax(2) = subplot(3,1,2);
plot(t,y)
axis([-inf inf -inf inf])
set(gca,'tickdir','out')
box off

ax(3) = subplot(3,1,3);
[~, pf] = GetWaveFreqTimeseries(Wxy,freq);

plot(t,pf.freq,'r.')
axis([-inf inf 0 inf])
box off
set(gca,'tickdir','out')

linkaxes(ax,'x')
