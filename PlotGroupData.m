%% Multiband XWT of curvature timeseries

% trl = 1;
% y = curv_trl(:,trl);
% dj = 1/24;
% 
% Wxy = ComputeMultibandXWT([y(:), y(:)], time_trl(:,1),'freqRange',[15 80],'dj',1/48,'timeRange',[0 0.6]);

%% Make some plots from group data

stimNum = 1; %(1 = Dark, 2 = Vib)
bendNum = 3;

blah = squeeze(grpData.dataMat(1,stimNum,:,1,:,bendNum));
amp1 = abs(blah(:));
lowInds = find(amp1< 20);
amp1(lowInds) = [];

blah = squeeze(grpData.dataMat(1,stimNum,:,2,:,bendNum));
per1 = blah(:);
per1(lowInds) = [];

blah = squeeze(grpData.dataMat(2,stimNum,:,1,:,bendNum));
amp2 = abs(blah(:));
lowInds = find(amp2 < 20);
amp2(lowInds) = [];

blah = squeeze(grpData.dataMat(2,stimNum,:,2,:,bendNum));
per2 = blah(:);
per2(lowInds)=[];

figure
% scatter(ctrl.vib.per{1},ctrl.vib.amp{1},'.')
scatter(per1,amp1,'r.')
hold on
% scatter(abl.vib.per{1},abl.vib.amp{1},'ro');
scatter(per2,amp2,'go')
xlim([0 80])
% xlim([-inf inf])
ylim([0 360])
% ylim([-inf inf])
box off
set(gca,'tickdir','out','color','k')
xlabel('Bend period (ms)')
ylabel('Bend amplitude (deg)')
title(['Bend # ' num2str(bendNum)])
shg



%% Combined bend period vs bend amp for group data

stimNum = 1; %(1 = Dark, 2 = Vib)
bendNums = [5 6]; % Bend nums to combine
bendNum_amp = bendNums(1);
xLim = [0 100];
yLim = [0 360];

bendNums = unique(bendNums);
if any(diff(bendNums)~=1)
    error('Bends to be combined must be consecutive')
end

dataMat = grpData.dataMat;
bendPer = cell(2,1);
bendAmp = bendPer;
perVec = xLim(1):xLim(2);
ampVec = yLim(1):yLim(2);
paMat = zeros(length(ampVec),length(perVec),2);
gKer = gausswin(60)*gausswin(30)';
gKer = gKer/sum(gKer(:));
for grp = 1:2
    nanInds =[];
    per = [];
    for bend = 1:numel(bendNums);
        blah = squeeze(dataMat(grp,stimNum,:,2,:,bendNums(bend)));
        blah = blah(:);
        nanInds = union(nanInds,find(isnan(blah)));
        per = [per, blah];
    end   
    bendPer{grp} = sum(per,2)/size(per,2);    
    blah = squeeze(dataMat(grp,stimNum,:,1,:,bendNum_amp));
    blah = abs(blah(:));
    lowInds = find(blah < 20);
    highInds = find(blah > 360);
    remInds = union(lowInds,nanInds);
    remInds = union(remInds,highInds);
    blah(remInds) = [];
    bendPer{grp}(remInds)=[];
    bendAmp{grp} = blah;
    
    
    x = round(bendPer{grp});
    y = round(bendAmp{grp});
    tempMat = paMat(:,:,grp);
    for jj = 1:numel(x);
        try
        tempMat(y(jj),x(jj)) = tempMat(y(jj),x(jj))+1;
        catch
            disp('Error: try changing xLim or yLim to fit all points!')
            tempMat(y(jj),x(jj)) = tempMat(y(jj),x(jj))+1;
        end
    end
    paMat(:,:,grp) = Standardize(conv2(tempMat,gKer,'same'));
end

mergedMap = cat(3,paMat(:,:,2),paMat(:,:,1),zeros(size(paMat(:,:,1))));


figure('Name','Tail deflection vs period ')
sh = scatter(bendPer{1},bendAmp{1},'g.');
%     alpha(sh,0.5);
hold on
muX = mean(bendPer{1});
sigX = std(bendPer{1});
rX = [muX-sigX muX + sigX];
muY = mean(bendAmp{1});
sigY = std(bendAmp{1});
rY = [muY-sigY muY + sigY];
plot(rX,muY*ones(size(rX)),'y-','linewidth',2)
plot(muX*ones(size(rY)),rY,'y-','linewidth',2)

sh = scatter(bendPer{2},bendAmp{2},'ro');
muX = mean(bendPer{2});
sigX = std(bendPer{2});
rX = [muX-sigX muX + sigX];
muY = mean(bendAmp{2});
sigY = std(bendAmp{2});
rY = [muY-sigY muY + sigY];
plot(rX,muY*ones(size(rX)),'c-','linewidth',2)
plot(muX*ones(size(rY)),rY,'c-','linewidth',2)

%     alpha(sh,0.5)
set(gca,'color','k','tickdir','out')
box off
xlim(xLim)
% xlim([-inf inf])
ylim(yLim)
% ylim([-inf inf])
xlabel('Bend period (ms)')
ylabel('Bend amplitude (deg)')
title(['Bend # ' num2str(bendNums)])
shg

% fName = [data.ablationType '_bendNums_' num2str(bendNums) '_' datestr(now,30)];
% SaveImageStack(paMat,path,fName);


%### Bend amp vs bend per, color density maps - juxtaposed
figure('Name','Bend amp vs per, color density maps - juxtaposed')
subaxis(1,2,1,'SH',0), imagesc(perVec,ampVec, paMat(:,:,1)), box off, set(gca,'tickdir','out','ydir','normal','xtick',[0:20:80])
xlabel('Beat period (ms)')
ylabel('Beat amplitude (deg)')
title('Control')
subaxis(1,2,2), imagesc(perVec, ampVec, paMat(:,:,2)), box off, set(gca,'tickdir','out','ydir','normal','ytick',[],'xtick',[0:20:80])
title('Ablated')


%### Bend amp vs bend per, color density maps - fused
figure('Name','Bend amp vs per, color density maps - fused')
imagesc(perVec,ampVec,mergedMap), box off, set(gca,'tickdir','out','ydir','normal')
xlabel('Beat period (ms)')
ylabel('Beat amp (deg)')
title('Beat amp vs per, density maps, fused')


%% Subplots - N bends, amp vs per
stimNum = 1; %(1 = Dark, 2 = Vib)

figure('Name','Multiple bend comparison for ctrl and abl (inter RS)')
nBends = 12;
nRows = 3;
nCols = 4;
xLim = [0 100];
yLim = [0 420];
for bendNum = 1:nBends
    blah = squeeze(grpData.dataMat(1,stimNum,:,1,:,bendNum));
    amp1 = abs(blah(:));
    lowInds = find(amp1< 20);
    amp1(lowInds) = [];
    
    blah = squeeze(grpData.dataMat(1,stimNum,:,2,:,bendNum));
    per1 = blah(:);
    per1(lowInds) = [];
    
    blah = squeeze(grpData.dataMat(2,stimNum,:,1,:,bendNum));
    amp2 = abs(blah(:));
    lowInds = find(amp2 < 20);
    amp2(lowInds) = [];
    
    blah = squeeze(grpData.dataMat(2,stimNum,:,2,:,bendNum));
    per2 = blah(:);
    per2(lowInds)=[];
    
    subaxis(nRows,nCols,bendNum,'SpacingVert',0.05,'MR',0.05,'ML',0.1,'SpacingHoriz',0.01)
    scatter(per1,amp1,'r.')
    hold on
    scatter(per2,amp2,'g.')
    xlim(xLim)
    ylim(yLim)
    box off
    set(gca,'tickdir','out','color','k')
     hold off
    if bendNum ==1
        ylabel('Bend amp (deg)')
        title(['Bend # ' num2str(bendNum)])
        set(gca,'xtick',[])
    else
        title(num2str(bendNum))
        set(gca,'ytick',[])
    end
    
    if bendNum == nBends
        xlabel('Bend per (ms)')
        set(gca,'ytick',[])  
    else
        set(gca,'xtick',[])
    end
end

%% Subplots - Dynamic bend amp and per
stimNum = 2; %(1 = Dark, 2 = Vib)
paramNum = 2; % (1 = Bend amp, 2 = Bend per)
yLim = [0 60];

figure('Name','Dynamic bend amp for ctrl and abl (inter RS)')
for grpNum = 1:2
    subplot(2,1,grpNum)
    if grpNum ==1
        lbl = 'Control';
        clr = 'r';
    else
        lbl = 'Ablated';
        clr = 'g';

    end
    for paramNum = 2
        if paramNum ==1
            lbl = 'Bend Amplitude';
        elseif paramNum ==2
            lbl = 'Bend Period';
        end
%         title(lbl)       
        for pkNum = 1:size(dataMat,6)
            var = squeeze(dataMat(grpNum,stimNum,:,paramNum,:,pkNum));
            var = var(:);
            plot(pkNum,var,'.','color',clr)
            hold on
            set(gca,'color','k','tickdir','out')
            box off
            xlim([0 size(dataMat,6)+1])
            ylim(yLim)
        end
    end   
end
xlabel('Bend Num')
shg















