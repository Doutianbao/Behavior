
%% Vib Stim
% --- Unraveling trls ---
stimType = 'vib';
disp(stimType)
var = {};
var{1} = data.ctrl.(stimType).procData;
var{2} = data.abl.(stimType).procData;
W_all = cell(size(var));
wDims = [];
for trtmnt = 1:2
    disp(['Stim # ' num2str(trtmnt)])
    for fishNum = 1:length(var{trtmnt})
        disp(['Fish # ' num2str(fishNum)])
        W = var{trtmnt}{fishNum}.W;
        if fishNum == 1
            wDims = size(W.curv.coeff{1});
            nTrls = numel(W.trlList);
            W = [W.curv.coeff{:}];
            W = reshape(W,[wDims nTrls]);
            W_all{trtmnt} =  W;
        else
            nTrls = numel(W.trlList);
            W = [W.curv.coeff{:}];
            W = reshape(W,[wDims nTrls]);
            W_all{trtmnt} =  cat(3,W_all{trtmnt},W);
        end
    end
end

% --- Correlations ---
tic
disp('Getting correlations...')
disp('Ctrl')
D = struct;
nTrls = size(W_all{1},3);
blah = nan(nTrls,2);
trlList = 1:nTrls;
for trl = trlList(:)'
    w1 = abs(W_all{1}(:,:,trl));
    otherTrls = setdiff(trlList,trl);
    w_mu1 = mean(abs(W_all{1}(:,:,otherTrls)),3);
    w_mu2 = mean(abs(W_all{2}),3);
    blah(trl,1) = corr2(w1,w_mu1);
    blah(trl,2) = corr2(w1,w_mu2);
end
data.ctrl.(stimType).corrMat = blah;
toc

disp('Getting correlations...')
disp('Abl')
nTrls = size(W_all{2},3);
blah = nan(nTrls,2);
trlList = 1:nTrls;
for trl = trlList(:)'
    w1 = abs(W_all{2}(:,:,trl));
    otherTrls = setdiff(trlList,trl);
    w_mu1 = mean(abs(W_all{2}(:,:,otherTrls)),3);
    w_mu2 = mean(abs(W_all{1}),3);
    blah(trl,1) = corr2(w1,w_mu1);
    blah(trl,2) = corr2(w1,w_mu2);
end
data.abl.(stimType).corrMat = blah;
toc


%% Dark flash stim
% --- Unraveling trls ---
stimType = 'dark';
var = {};
var{1} = data.ctrl.(stimType).procData;
var{2} = data.abl.(stimType).procData;

W_all = cell(size(var));
wDims = [];
for trtmnt = 1:2
    disp(['Stim # ' num2str(trtmnt)])
    for fishNum = 1:length(var{trtmnt})
        disp(['Fish # ' num2str(fishNum)])
        W = var{trtmnt}{fishNum}.W;
        if fishNum == 1
            wDims = size(W.curv.coeff{1});
            nTrls = numel(W.trlList);
            W = [W.curv.coeff{:}];
            W = reshape(W,[wDims nTrls]);
            W_all{trtmnt} =  W;
        else
            nTrls = numel(W.trlList);
            W = [W.curv.coeff{:}];
            W = reshape(W,[wDims nTrls]);
            W_all{trtmnt} =  cat(3,W_all{trtmnt},W);
        end
    end
end

% --- Correlations ---
tic
disp('Getting correlations...')
disp('Ctrl')
D = struct;
nTrls = size(W_all{1},3);
blah = nan(nTrls,2);
trlList = 1:nTrls;
for trl = trlList(:)'
    w1 = abs(W_all{1}(:,:,trl));
    otherTrls = setdiff(trlList,trl);
    w_mu1 = mean(abs(W_all{1}(:,:,otherTrls)),3);
    w_mu2 = mean(abs(W_all{2}),3);
    blah(trl,1) = corr2(w1,w_mu1);
    blah(trl,2) = corr2(w1,w_mu2);
end
data.ctrl.(stimType).corrMat = blah;
toc

disp('Getting correlations...')
disp('Abl')
nTrls = size(W_all{2},3);
blah = nan(nTrls,2);
trlList = 1:nTrls;
for trl = trlList(:)'
    w1 = abs(W_all{2}(:,:,trl));
    otherTrls = setdiff(trlList,trl);
    w_mu1 = mean(abs(W_all{2}(:,:,otherTrls)),3);
    w_mu2 = mean(abs(W_all{1}),3);
    blah(trl,1) = corr2(w1,w_mu1);
    blah(trl,2) = corr2(w1,w_mu2);
end
data.abl.(stimType).corrMat = blah;
toc

