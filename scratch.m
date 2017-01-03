
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


out = AnalyzeFreeSwims_nCycles_batch(paths(1), 'xLim',[0 750],'paramList',...
    {'headAmp','bodyAmp'});
