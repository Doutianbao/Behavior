
tic
nFramesInTrl = 750;
fps = 500;

time = (0:size(tailCurv,3)-1)*(1/fps);
nTrls = size(tailCurv,3)/nFramesInTrl;

% tc = zeros(size(tailCurv));
% x = gradient(squeeze(tailCurv(:,1,:)),2);
% tc(:,1,:) = x;
% y = gradient(squeeze(tailCurv(:,2,:)),2);
% tc(:,2,:) = y;


disp('Calculating tail angles...')
% tailCurv = tailCurv_uncorrected;
tA = GetTailTangents(tailCurv);
tA_trl = reshape(tA,size(tA,1),nFramesInTrl,nTrls);
time_trl = reshape(time,nFramesInTrl,nTrls);

tA_5 = GetTailTangents(tailCurv,5);
curv = tA_5(end,:)';
curv_trl = reshape(curv,nFramesInTrl,nTrls);
toc


%%
% trls = [1:17 19];
trls = 1:nTrls;
xLim = [-inf 750];
cLim = [-225 225];
yShift = 300;
% blah = zeros(size(tA_trl,1),size(tA_trl,2)*numel(trls));
blah = tA_trl(:,:,trls(1));
nRows = size(tA_trl,1);
yTick = cumsum(repmat(nRows,numel(trls),1))+round(nRows/2)-nRows;
yTL = cell(numel(yTick),1);
yTL{1} = num2str(trls(1));
for jj = 2:numel(trls)
    temp = squeeze(tA_trl(:,:,trls(jj)));
    blah = cat(1,blah,temp);
    yTL{jj} = num2str(trls(jj));
end
blah(abs(blah)<=10)=0;
% blah = round(blah/10)*10;
figure
cMap = jet(64*4);
imagesc(blah),colormap(cMap)
box off
colorbar
xTick = [50 250:250:nFramesInTrl];
xtl = (xTick/fps)*1000;
set(gca,'ytick',yTick,'yticklabel',yTL,'xtick',xTick,'xTickLabel',xtl,'tickdir','out','clim',cLim)
xlim(xLim)
xlabel('Time (ms)')
ylabel('Trl #')
title('Tail undulations for different trials')


figure('Name','Trlzed curv ts')
yTick = zeros(numel(trls),1);
count = 0;
for trl = trls
    count = count + 1;
    yTick(count) = yShift*(count-1);
    if mod(count,2)==0
        plot(time_trl*1000,curv_trl(:,trl)-yTick(count),'r')
        hold on
    else
        plot(time_trl*1000,curv_trl(:,trl)-yTick(count),'b')
        hold on
    end
end
xlim((xLim/fps)*1000)
ylim([-yTick(end)-yShift yShift])
xlabel('Time (ms)')
set(gca,'tickdir','out','ytick',flipud(-yTick),'yticklabel',fliplr(trls),'xtick',[100 500:500:time_trl(end)*1000])
box off
title(['Timeseries of total tail bend angle for different trials (yShift = ' num2str(yShift) '^o)'])