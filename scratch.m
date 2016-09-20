%% Multiband XWT of curvature timeseries

% trl = 1;
% y = curv_trl(:,trl);
% dj = 1/24;
% 
% Wxy = ComputeMultibandXWT([y(:), y(:)], time_trl(:,1),'freqRange',[15 80],'dj',1/48,'timeRange',[0 0.6]);

%% Make some plots from group data

stimNum = 2; %(1 = Dark, 2 = Vib)
bendNum = 1;

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
ylim([0 360])
box off
set(gca,'tickdir','out','color','k')
xlabel('Bend period (ms)')
ylabel('Bend amplitude (deg)')
shg


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
paramNum = 1; % (1 = Bend amp, 2 = Bend per)

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
            ylim([-inf inf])
        end
    end   
end
xlabel('Bend Num')
shg















