
% A_pos = data.abl.vib.mean >= (data.ctrl.vib.mean + 0.5*data.ctrl.vib.std);
% A_neg = data.abl.vib.mean <= (data.ctrl.vib.mean - 0.5*data.ctrl.vib.std);
% 
% A_pn = (A_pos | A_neg);

% w = data.ctrl.vib.procData{1}.W; 
% time = w.time; freq = w.freq;
% A = (data.abl.vib.mean - data.ctrl.vib.mean)./(data.ctrl.vib.std);
% M = data.ctrl.vib.mean > mean(data.ctrl.vib.mean(:));
% A = A.*M;
% A(abs(A) <0.75) = 0;
% figure
% imagesc(time,freq,A), set(gca,'ydir','normal')
% box off
% xlabel('Time (ms)')
% ylabel('Freq (Hz)')
% colorbar


%% Working with grpData that contains ablation data from all ablation groups to make some intra- and inter-group comparisons

[intra,inter] = deal(struct);

% --- Intra group diff maps for vib stim ---
% --- mHom group
intra.vib.mHom_inter = grpData.inter.ctrl.vib.mean - grpData.mHom.ctrl.vib.mean;
intra.vib.mHom_vent = grpData.vent.ctrl.vib.mean - grpData.mHom.ctrl.vib.mean;

% --- Intermediate group
intra.vib.inter_vent = grpData.vent.ctrl.vib.mean - grpData.inter.ctrl.vib.mean;
intra.vib.inter_mHom = grpData.mHom.ctrl.vib.mean - grpData.inter.ctrl.vib.mean;

% --- Ventral group
intra.vib.vent_mHom = grpData.mHom.ctrl.vib.mean - grpData.vent.ctrl.vib.mean;
intra.vib.vent_inter = grpData.inter.ctrl.vib.mean - grpData.vent.ctrl.vib.mean;

% --- Inter group diff maps for vib stim ---
% --- mHom group
inter.vib.mHom_mHom = grpData.mHom.abl.vib.mean - grpData.mHom.ctrl.vib.mean;
inter.vib.mHom_inter = grpData.mHom.abl.vib.mean - grpData.inter.ctrl.vib.mean;
inter.vib.mHom_vent = grpData.mHom.abl.vib.mean - grpData.vent.ctrl.vib.mean;

% --- Intermediate group
inter.vib.inter_inter = grpData.inter.abl.vib.mean - grpData.inter.ctrl.vib.mean;
inter.vib.inter_vent = grpData.inter.abl.vib.mean - grpData.vent.ctrl.vib.mean;
inter.vib.inter_mHom = grpData.inter.abl.vib.mean - grpData.mHom.ctrl.vib.mean;

% --- Ventral group
inter.vib.vent_vent = grpData.vent.abl.vib.mean - grpData.vent.ctrl.vib.mean;
inter.vib.vent_mHom = grpData.vent.abl.vib.mean - grpData.mHom.ctrl.vib.mean;
inter.vib.vent_inter = grpData.vent.abl.vib.mean - grpData.inter.ctrl.vib.mean;


% --- Intra group diff maps for dark stim ---
% --- mHom group
intra.dark.mHom_inter = grpData.inter.ctrl.dark.mean - grpData.mHom.ctrl.dark.mean;
intra.dark.mHom_vent = grpData.vent.ctrl.dark.mean - grpData.mHom.ctrl.dark.mean;

% --- Intermediate group
intra.dark.inter_vent = grpData.vent.ctrl.dark.mean - grpData.inter.ctrl.dark.mean;
intra.dark.inter_mHom = grpData.mHom.ctrl.dark.mean - grpData.inter.ctrl.dark.mean;

% --- Ventral group
intra.dark.vent_mHom = grpData.mHom.ctrl.dark.mean - grpData.vent.ctrl.dark.mean;
intra.dark.vent_inter = grpData.inter.ctrl.dark.mean - grpData.vent.ctrl.dark.mean;

% --- Inter group diff maps for vib stim ---
% --- mHom group
inter.dark.mHom_mHom = grpData.mHom.abl.dark.mean - grpData.mHom.ctrl.dark.mean;
inter.dark.mHom_inter = grpData.mHom.abl.dark.mean - grpData.inter.ctrl.dark.mean;
inter.dark.mHom_vent = grpData.mHom.abl.dark.mean - grpData.vent.ctrl.dark.mean;

% --- Intermediate group
inter.dark.inter_inter = grpData.inter.abl.dark.mean - grpData.inter.ctrl.dark.mean;
inter.dark.inter_vent = grpData.inter.abl.dark.mean - grpData.vent.ctrl.dark.mean;
inter.dark.inter_mHom = grpData.inter.abl.dark.mean - grpData.mHom.ctrl.dark.mean;

% --- Ventral group
inter.dark.vent_vent = grpData.vent.abl.dark.mean - grpData.vent.ctrl.dark.mean;
inter.dark.vent_mHom = grpData.vent.abl.dark.mean - grpData.mHom.ctrl.dark.mean;
inter.dark.vent_inter = grpData.vent.abl.dark.mean - grpData.inter.ctrl.dark.mean;


%% Plotting

%---- ---- VIB STIM --- ---
stimType = 'vib';
cl = [-100 100];
xl = [-50 300];

W = grpData.mHom.ctrl.(stimType).procData{1}.W;
freq = W.freq;
time = W.time;

%--- mHom group variation
W_inter_mHom = cat(3,inter.(stimType).mHom_mHom, inter.(stimType).mHom_inter, inter.(stimType).mHom_vent);
M = mean(W_inter_mHom,3);
[r_max,c_max] = find(abs(M)==max(abs(M(:))));
ax = {};
ax{1} =[1 0.2 0 0.8];
ax{2} = [1 0.2 0 0.6];
ax{3} = [1 0.2 0 0.4];
ax{4} = [1 0.2 0 0.2];
ax{5} = [1 0.2 0 0];

fh = figure('Name','Variation, mHom');
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3},ax{4},ax{5});
axes(axH(1))
imagesc(time,freq,W_inter_mHom(:,:,1))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'mHom abl - mHom ctrl','Freq (Hz)'})
xlim(xl)
title('mHom ablation group')

axes(axH(2))
imagesc(time,freq,W_inter_mHom(:,:,2))
hold on
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'mHom abl - Inter ctrl','Freq (Hz)'})
xlim(xl)

axes(axH(3))
imagesc(time,freq,W_inter_mHom(:,:,3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'mHom abl - Vent ctrl','Freq (Hz)'})
xlim(xl)

axes(axH(4))
imagesc(time,freq,mean(W_inter_mHom,3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Avg','Freq (Hz)'})
xlim(xl)

axes(axH(5))
imagesc(time,freq,std(W_inter_mHom,[],3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out')
box off
ylabel({'Std','Freq (Hz)'})
xlabel('Time (ms)')
xlim(xl)
cbh = colorbar;
set(cbh,'color','k','ytick',[-200 0 200],'location','East')


%--- Intermediate group variation
W_inter_inter = cat(3,inter.(stimType).inter_mHom, inter.(stimType).inter_inter, inter.(stimType).inter_vent);
M = mean(W_inter_inter,3);
[r_max,c_max] = find(abs(M)==max(abs(M(:))));
ax = {};
ax{1} =[1 0.2 0 0.8];
ax{2} = [1 0.2 0 0.6];
ax{3} = [1 0.2 0 0.4];
ax{4} = [1 0.2 0 0.2];
ax{5} = [1 0.2 0 0];

fh = figure('Name','Variation, Intermediate');
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3},ax{4},ax{5});
axes(axH(1))
imagesc(time,freq,W_inter_inter(:,:,1))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Inter abl - mHom ctrl','Freq (Hz)'})
xlim(xl)
title('Intermediate ablation group')


axes(axH(2))
imagesc(time,freq,W_inter_inter(:,:,2))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Inter abl - Inter ctrl','Freq (Hz)'})
xlim(xl)

axes(axH(3))
imagesc(time,freq,W_inter_inter(:,:,3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Inter abl - Vent ctrl','Freq (Hz)'})
xlim(xl)

axes(axH(4))
imagesc(time,freq,mean(W_inter_inter,3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Avg','Freq (Hz)'})
xlim(xl)

axes(axH(5))
imagesc(time,freq,std(W_inter_inter,[],3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out')
box off
ylabel({'Std','Freq (Hz)'})
xlabel('Time (ms)')
xlim(xl)
cbh = colorbar;
set(cbh,'color','k','ytick',[-200 0 200],'location','East')

%--- Ventral group variation
W_inter_vent = cat(3,inter.(stimType).vent_mHom, inter.(stimType).vent_inter, inter.(stimType).vent_vent);
M = mean(W_inter_vent,3);
[r_max,c_max] = find(abs(M)==max(abs(M(:))));
ax = {};
ax{1} =[1 0.2 0 0.8];
ax{2} = [1 0.2 0 0.6];
ax{3} = [1 0.2 0 0.4];
ax{4} = [1 0.2 0 0.2];
ax{5} = [1 0.2 0 0];

fh = figure('Name','Variation, Ventral');
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3},ax{4},ax{5});
axes(axH(1))
imagesc(time,freq,W_inter_vent(:,:,1))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Vent abl - mHom ctrl','Freq (Hz)'})
xlim(xl)
title('Ventral ablation group')

axes(axH(2))
imagesc(time,freq,W_inter_vent(:,:,2))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Vent abl - Inter ctrl','Freq (Hz)'})
xlim(xl)

axes(axH(3))
imagesc(time,freq,W_inter_vent(:,:,3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Vent abl - Vent ctrl','Freq (Hz)'})
xlim(xl)

axes(axH(4))
imagesc(time,freq,mean(W_inter_vent,3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Avg','Freq (Hz)'})
xlim(xl)

axes(axH(5))
imagesc(time,freq,std(W_inter_vent,[],3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out')
box off
ylabel({'Std','Freq (Hz)'})
xlabel('Time (ms)')
xlim(xl)
cbh = colorbar;
set(cbh,'color','k','ytick',[-200 0 200],'location','East')

%---- ---- DARK FLASH STIM --- ---
stimType = 'dark';
mult = 1;
xl = [-50 300];

W = grpData.mHom.ctrl.(stimType).procData{1}.W;
freq = W.freq;
time = W.time;

%--- mHom group variation
W_inter_mHom = cat(3,inter.(stimType).mHom_mHom, inter.(stimType).mHom_inter, inter.(stimType).mHom_vent)*mult;
M = mean(W_inter_mHom,3);
[r_max,c_max] = find(abs(M)==max(abs(M(:))));
ax = {};
ax{1} =[1 0.2 0 0.8];
ax{2} = [1 0.2 0 0.6];
ax{3} = [1 0.2 0 0.4];
ax{4} = [1 0.2 0 0.2];
ax{5} = [1 0.2 0 0];

fh = figure('Name','Variation, mHom');
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3},ax{4},ax{5});
axes(axH(1))
imagesc(time,freq,W_inter_mHom(:,:,1))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'mHom abl - mHom ctrl','Freq (Hz)'})
xlim(xl)
title('mHom ablation group')

axes(axH(2))
imagesc(time,freq,W_inter_mHom(:,:,2))
hold on
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'mHom abl - Inter ctrl','Freq (Hz)'})
xlim(xl)

axes(axH(3))
imagesc(time,freq,W_inter_mHom(:,:,3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'mHom abl - Vent ctrl','Freq (Hz)'})
xlim(xl)

axes(axH(4))
imagesc(time,freq,mean(W_inter_mHom,3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Avg','Freq (Hz)'})
xlim(xl)

axes(axH(5))
imagesc(time,freq,std(W_inter_mHom,[],3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out')
box off
ylabel({'Std','Freq (Hz)'})
xlabel('Time (ms)')
xlim(xl)
cbh = colorbar;
set(cbh,'color','k','ytick',[-200 0 200],'location','East')


%--- Intermediate group variation
W_inter_inter = cat(3,inter.(stimType).inter_mHom, inter.(stimType).inter_inter, inter.(stimType).inter_vent)*mult;
M = mean(W_inter_inter,3);
[r_max,c_max] = find(abs(M)==max(abs(M(:))));
ax = {};
ax{1} =[1 0.2 0 0.8];
ax{2} = [1 0.2 0 0.6];
ax{3} = [1 0.2 0 0.4];
ax{4} = [1 0.2 0 0.2];
ax{5} = [1 0.2 0 0];

fh = figure('Name','Variation, Intermediate');
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3},ax{4},ax{5});
axes(axH(1))
imagesc(time,freq,W_inter_inter(:,:,1))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Inter abl - mHom ctrl','Freq (Hz)'})
xlim(xl)
title('Intermediate ablation group')


axes(axH(2))
imagesc(time,freq,W_inter_inter(:,:,2))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Inter abl - Inter ctrl','Freq (Hz)'})
xlim(xl)

axes(axH(3))
imagesc(time,freq,W_inter_inter(:,:,3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Inter abl - Vent ctrl','Freq (Hz)'})
xlim(xl)

axes(axH(4))
imagesc(time,freq,mean(W_inter_inter,3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Avg','Freq (Hz)'})
xlim(xl)

axes(axH(5))
imagesc(time,freq,std(W_inter_inter,[],3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out')
box off
ylabel({'Std','Freq (Hz)'})
xlabel('Time (ms)')
xlim(xl)
cbh = colorbar;
set(cbh,'color','k','ytick',[-200 0 200],'location','East')

%--- Ventral group variation
W_inter_vent = cat(3,inter.(stimType).vent_mHom, inter.(stimType).vent_inter, inter.(stimType).vent_vent)*mult;
M = mean(W_inter_vent,3);
[r_max,c_max] = find(abs(M)==max(abs(M(:))));
ax = {};
ax{1} =[1 0.2 0 0.8];
ax{2} = [1 0.2 0 0.6];
ax{3} = [1 0.2 0 0.4];
ax{4} = [1 0.2 0 0.2];
ax{5} = [1 0.2 0 0];

fh = figure('Name','Variation, Ventral');
axH = CreateSubaxes(fh,ax{1},ax{2},ax{3},ax{4},ax{5});
axes(axH(1))
imagesc(time,freq,W_inter_vent(:,:,1))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Vent abl - mHom ctrl','Freq (Hz)'})
xlim(xl)
title('Ventral ablation group')

axes(axH(2))
imagesc(time,freq,W_inter_vent(:,:,2))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Vent abl - Inter ctrl','Freq (Hz)'})
xlim(xl)

axes(axH(3))
imagesc(time,freq,W_inter_vent(:,:,3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Vent abl - Vent ctrl','Freq (Hz)'})
xlim(xl)

axes(axH(4))
imagesc(time,freq,mean(W_inter_vent,3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out','xtick',[])
box off
ylabel({'Avg','Freq (Hz)'})
xlim(xl)

axes(axH(5))
imagesc(time,freq,std(W_inter_vent,[],3))
hold on
plot([time(1) time(end)],freq([r_max r_max]),'w:')
plot(time([c_max c_max]),[freq(1), freq(end)],'w:')
set(gca,'clim',cl,'ydir','normal','tickdir','out')
box off
ylabel({'Std','Freq (Hz)'})
xlabel('Time (ms)')
xlim(xl)
cbh = colorbar;
set(cbh,'color','k','ytick',[-200 0 200],'location','East')


