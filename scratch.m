
var1 = data.ctrl.vib;
var2 = data.abl.vib;

nFish = length(var1);
W = cell(nFish,1);
for fn = 2
    disp(['Fish # ' num2str(fn)])
    disp('Reading W... ')
    W = var1{fn}.W;
    PlotWTs(W)
    a = 1;
end


% nFish = length(var2);
% W = cell(nFish,1);
% for fn = 4
%     disp(['Fish # ' num2str(fn)])
%     disp('Reading W... ')
%     W = var2{fn}.W;
%     PlotWTs(W)
%     a = 1;
% end