%% The average WT across fish from average WTs for fish across trials
tic
W_all = cell(2,1);

%--- Ctrl ---
disp('Avg WTs from all control fish...')
var = data.ctrl.vib;
nFish = length(var);
for fishNum = 1:nFish
    disp(['Fish # ' num2str(fishNum)])
    W = var{fishNum}.W;
    W = cat(1,W.head.avg,W.tail.avg);
    if fishNum ==1
        W_all{1} = W;
    else
        W_all{1} = cat(3,W_all{1},W);
    end    
end
data.ctrl.vib_mean = mean(W_all{1},3);
data.ctrl.vib_std = std(W_all{1},[],3);
data.ctrl.vib_cv = data.ctrl.vib_std./(data.ctrl.vib_mean + 1);
corrVec = nan(1,nFish);
for fishNum = 1:nFish
    corrVec(fishNum) = corr2(data.ctrl.vib_mean,W_all{1}(:,:,fishNum));
end
data.ctrl.vib_corrVec = corrVec;

%--- Abl ---
disp('Avg WTs from all ablated fish...')
var = data.abl.vib;
nFish = length(var);
for fishNum = 1:nFish
    disp(['Fish # ' num2str(fishNum)])
    W = var{fishNum}.W;
    W = cat(1,W.head.avg,W.tail.avg);
    if fishNum ==1
        W_all{2} = W;
    else
        W_all{2} = cat(3,W_all{2},W);
    end    
end

data.abl.vib_mean = mean(W_all{2},3);
data.abl.vib_std = std(W_all{2},[],3);
data.abl.vib_cv = data.abl.vib_std./(data.abl.vib_mean + 1);
corrVec = nan(1,nFish);
for fishNum = 1:nFish
    corrVec(fishNum) = corr2(data.abl.vib_mean,W_all{2}(:,:,fishNum));
end
data.abl.vib_corrVec = corrVec;
toc
