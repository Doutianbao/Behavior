M = struct;
S = struct;

C = cat(3, W.ctrl.mHom.vib, W.ctrl.inter.vib, W.ctrl.vent.vib);
M.ctrl.vib = mean(C,3);
S.ctrl.vib = std(C,[],3);

M.abl.mHom.vib  = mean(W.abl.mHom.vib,3);
S.abl.mHom.vib = std(W.abl.mHom.vib,[],3);

thr_ctrl = M.ctrl.vib + 1 *S.ctrl.vib;
thr_abl = M.abl.mHom.vib + 1*S.abl.mHom.vib;

inds1 = find(M.abl.mHom.vib < thr_ctrl);
inds2 = find(M.ctrl.vib < thr_abl);

B1 = M.ctrl.vib; B1(inds1)=0;
B2 = M.abl.mHom.vib; B2(inds2)=0;



