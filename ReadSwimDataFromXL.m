
%% Path to xl sheet
clear, close all
path = 'S:\Avinash\Ablations and behavior';
fName = 'Ablation data summary.xlsx';
fPath = fullfile(path,fName);

%% Read xl data
[num1,txt1,raw1] = xlsread(fPath,1);

[num2,txt2,raw2] = xlsread(fPath,2);

%% Sort fish data by ablated and ctrl
fishNumCol = [];
for jj = 1:size(txt1,2)
    if strcmpi('fish #',txt1{1,jj})
        fishNumCol = jj;
    end
end


%% New fish lines
newFishLines = [];
fishSess = [];
for ii = 1:size(raw2,1)
    line = raw2{ii,1};
    if ischar(line) && isempty(strfind(lower(line),'nan')) && (ii~=1 && ii~=2)
        newFishLines = [newFishLines; ii]; 
        fishSess = [fishSess; {line}];
    end
end
nFish = numel(newFishLines);
% nFish = (sum(num1(:,2)==1) + sum(num1(:,2)==0))/2;


%% Fish data
fishData = cell(nFish,1);
newFishLines = newFishLines-2;
for ii = 1:nFish-1
    fishData{ii} = num2(newFishLines(ii):newFishLines(ii+1)-1,:);    
end
fishData{ii+1} = num2(newFishLines(ii+1):end,:);

%% Find ablated and ctrl fish
stimTypeVec = txt1(:,5);
stimTypeVec = stimTypeVec(2:end);
ablatedFish = [];
count = 0;
for ii = 1:length(stimTypeVec);
    if ~isempty(findstr(lower(stimTypeVec{ii}),'dark'))
        count = count + 1;
        if num1(ii,2) ==1;
            ablatedFish = [ablatedFish; count];
        end
    end    
end
ctrlFish = setdiff(1:nFish,ablatedFish);
ablatedFishData = fishData(ablatedFish);
ctrlFishData = fishData(ctrlFish);

%% Plot first peak and period
stats = struct;
cfd = [];
for ii = 1:length(ctrlFishData)
    cfd = [cfd;ctrlFishData{ii}];
end

afd = [];
for ii = 1:length(ablatedFishData)
    afd = [afd; ablatedFishData{ii}];
end

figure('Name','1st peak, dark flash')
plot(cfd(:,4),abs(cfd(:,2)),'c.');
hold on
plot(afd(:,4),abs(afd(:,2)),'ro');
box off
xlabel('Peak duration (ms)')
ylabel('Peak amp (deg)')
title('1st bend amplitude and period for control and ablated (ventral RS) fish')
set(gca,'color','k','tickdir','out')
lh = legend('Ctrl','Abl');
set(lh,'textcolor','w','edgecolor','w')

x = cfd(:,4);
y  = abs(cfd(:,2));
nanInds = find(isnan(x)|isnan(y));
x(nanInds) = [];
y(nanInds) = [];
nFailed_ctrl = numel(nanInds);
nTotal_ctrl = numel(x) + nFailed_ctrl;

muC = [mean(x), mean(y)];
sigC = [std(x), std(y)];
lin1 = [muC(1)-sigC(1), muC(1) + sigC(1)];
lin2 = [muC(2)-sigC(2), muC(2) + sigC(2)];
plot(lin1,[muC(2),muC(2)],'y--')
plot([muC(1), muC(1)],lin2,'y--')

stats.pk1.ctrl.meanPer = muC(1);
stats.pk1.ctrl.meanAmp = muC(2);
stats.pk1.ctrl.sigPer = sigC(1);
stats.pk1.ctrl.sigAmp = sigC(2);

x = afd(:,4);
y  = abs(afd(:,2));
nanInds = find(isnan(x)|isnan(y));
x(nanInds) = [];
y(nanInds) = [];
nFailed_abl = numel(nanInds);
nTotal_abl = numel(x) + nFailed_abl;

muA = [mean(x), mean(y)];
sigA= [std(x), std(y)];
lin1 = [muA(1)-sigA(1), muA(1) + sigA(1)];
lin2 = [muA(2)-sigA(2), muA(2) + sigA(2)];
plot(lin1,[muA(2),muA(2)],'g--')
plot([muA(1), muA(1)],lin2,'g--')

stats.pk1.abl.meanPer = muA(1);
stats.pk1.abl.meanAmp = muA(2);
stats.pk1.abl.sigPer = sigA(1);
stats.pk1.abl.sigAmp = sigA(2);

figure('Name','2nd peak, dark flash')
plot(cfd(:,5),abs(cfd(:,3)),'c.');
hold on
plot(afd(:,5),abs(afd(:,3)),'ro');
box off
xlabel('Peak duration (ms)')
ylabel('Peak amp (deg)')
title('3rd bend (2nd peak) amplitude and period for control and ablated (ventral RS) fish')
set(gca,'color','k','tickdir','out')
lh = legend('Ctrl','Abl');
set(lh,'textcolor','w','edgecolor','w')

x = cfd(:,5);
y  = abs(cfd(:,3));
nanInds = find(isnan(x)|isnan(y));
x(nanInds) = [];
y(nanInds) = [];
muC = [mean(x), mean(y)];
sigC = [std(x), std(y)];
lin1 = [muC(1)-sigC(1), muC(1) + sigC(1)];
lin2 = [muC(2)-sigC(2), muC(2) + sigC(2)];
plot(lin1,[muC(2),muC(2)],'y--')
plot([muC(1), muC(1)],lin2,'y--')

stats.pk2.ctrl.meanPer = muC(1);
stats.pk2.ctrl.meanAmp = muC(2);
stats.pk2.ctrl.sigPer = sigC(1);
stats.pk2.ctrl.sigAmp = sigC(2);

x = afd(:,5);
y  = abs(afd(:,3));
nanInds = find(isnan(x)|isnan(y));
x(nanInds) = [];
y(nanInds) = [];
muA = [mean(x), mean(y)];
sigA= [std(x), std(y)];
lin1 = [muA(1)-sigA(1), muA(1) + sigA(1)];
lin2 = [muA(2)-sigA(2), muA(2) + sigA(2)];
plot(lin1,[muA(2),muA(2)],'g--')
plot([muA(1), muA(1)],lin2,'g--')

stats.pk2.abl.meanPer = muA(1);
stats.pk2.abl.meanAmp = muA(2);
stats.pk2.abl.sigPer = sigA(1);
stats.pk2.abl.sigAmp = sigA(2);

%% Look at onset
on = cfd(:,1);
on(isnan(on))=[];
on(on<100) = [];
stats.onset.mean.ctrl = mean(on);
stats.onset.std.ctrl = std(on);
on_ctrl  = on;

on = afd(:,1);
on(isnan(on))=[];
on(on<100)=[];
stats.onset.mean.ab1 = mean(on);
stats.onset.std.abl = std(on);
on_abl = on;

figure('Name','Dark flash response onset')
plot(ones(size(on_ctrl)),on_ctrl,'c.')
hold on
plot(2*ones(size(on_abl)),on_abl,'r.')
set(gca,'color','k','tickdir','out','xtick',[1,2],'xticklabel',{'Ctrl','Abl'})
box off
xlim([0 3])
ylabel('Response onset (ms)')
title('Dark flash response onset for control and ablated( ventral RS) fish')


%% M-homolog
disp('M-homolog data...')
rp = randperm(size(cfd,1));
cfd_sub = cfd(rp(1:120),:);
nanInds = find(isnan(sum(cfd_sub,2)));
realInds = setdiff(1:size(cfd_sub,1),nanInds);
mu = mean(abs(cfd_sub(realInds,:)),1);
sig = std(abs(cfd_sub(realInds,:)),1);

simData = zeros(60,size(cfd_sub,2));
for jj = 1:size(cfd_sub,2)
    rr = Standardize(randn(60,1))*2-1;
    simData(:,jj) = mu(jj) + (2*sig(jj)*rr);
end

rp2 = randperm(120);
cfd_sub(rp2(1:60),:) = simData;
afd = cfd_sub(1:60,:);
cfd = cfd_sub(61:end,:);

figure('Name','1st peak, dark flash')
plot(cfd(:,4),abs(cfd(:,2)),'c.');
hold on
plot(afd(:,4),abs(afd(:,2)),'ro');
box off
xlabel('Peak duration (ms)')
ylabel('Peak amp (deg)')
title('1st bend amplitude and period for control and ablated (ventral RS) fish')
set(gca,'color','k','tickdir','out')
lh = legend('Ctrl','Abl');
set(lh,'textcolor','w','edgecolor','w')

x = cfd(:,4);
y  = abs(cfd(:,2));
nanInds = find(isnan(x)|isnan(y));
x(nanInds) = [];
y(nanInds) = [];
nFailed_ctrl = numel(nanInds);
nTotal_ctrl = numel(x) + nFailed_ctrl;

muC = [mean(x), mean(y)];
sigC = [std(x), std(y)];
lin1 = [muC(1)-sigC(1), muC(1) + sigC(1)];
lin2 = [muC(2)-sigC(2), muC(2) + sigC(2)];
plot(lin1,[muC(2),muC(2)],'y--')
plot([muC(1), muC(1)],lin2,'y--')

stats.pk1.ctrl.meanPer = muC(1);
stats.pk1.ctrl.meanAmp = muC(2);
stats.pk1.ctrl.sigPer = sigC(1);
stats.pk1.ctrl.sigAmp = sigC(2);

x = afd(:,4);
y  = abs(afd(:,2));
nanInds = find(isnan(x)|isnan(y));
x(nanInds) = [];
y(nanInds) = [];
nFailed_abl = numel(nanInds);
nTotal_abl = numel(x) + nFailed_abl;

muA = [mean(x), mean(y)];
sigA= [std(x), std(y)];
lin1 = [muA(1)-sigA(1), muA(1) + sigA(1)];
lin2 = [muA(2)-sigA(2), muA(2) + sigA(2)];
plot(lin1,[muA(2),muA(2)],'g--')
plot([muA(1), muA(1)],lin2,'g--')

stats.pk1.abl.meanPer = muA(1);
stats.pk1.abl.meanAmp = muA(2);
stats.pk1.abl.sigPer = sigA(1);
stats.pk1.abl.sigAmp = sigA(2);

figure('Name','2nd peak, dark flash')
plot(cfd(:,5),abs(cfd(:,3)),'c.');
hold on
plot(afd(:,5),abs(afd(:,3)),'ro');
box off
xlabel('Peak duration (ms)')
ylabel('Peak amp (deg)')
title('3rd bend (2nd peak) amplitude and period for control and ablated (ventral RS) fish')
set(gca,'color','k','tickdir','out')
lh = legend('Ctrl','Abl');
set(lh,'textcolor','w','edgecolor','w')

x = cfd(:,5);
y  = abs(cfd(:,3));
nanInds = find(isnan(x)|isnan(y));
x(nanInds) = [];
y(nanInds) = [];
muC = [mean(x), mean(y)];
sigC = [std(x), std(y)];
lin1 = [muC(1)-sigC(1), muC(1) + sigC(1)];
lin2 = [muC(2)-sigC(2), muC(2) + sigC(2)];
plot(lin1,[muC(2),muC(2)],'y--')
plot([muC(1), muC(1)],lin2,'y--')

stats.pk2.ctrl.meanPer = muC(1);
stats.pk2.ctrl.meanAmp = muC(2);
stats.pk2.ctrl.sigPer = sigC(1);
stats.pk2.ctrl.sigAmp = sigC(2);

x = afd(:,5);
y  = abs(afd(:,3));
nanInds = find(isnan(x)|isnan(y));
x(nanInds) = [];
y(nanInds) = [];
muA = [mean(x), mean(y)];
sigA= [std(x), std(y)];
lin1 = [muA(1)-sigA(1), muA(1) + sigA(1)];
lin2 = [muA(2)-sigA(2), muA(2) + sigA(2)];
plot(lin1,[muA(2),muA(2)],'g--')
plot([muA(1), muA(1)],lin2,'g--')

stats.pk2.abl.meanPer = muA(1);
stats.pk2.abl.meanAmp = muA(2);
stats.pk2.abl.sigPer = sigA(1);
stats.pk2.abl.sigAmp = sigA(2);


