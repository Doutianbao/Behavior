
tic
nFramesInTrl = 750;
fps = 500;
stimTime = 100; % In ms

tailCurv = procData.tailCurv;


time = (0:size(tailCurv,3)-1)*(1/fps);
% time = time-stimTime/1000;
nTrls = size(tailCurv,3)/nFramesInTrl;

% tc = zeros(size(tailCurv));
% x = gradient(squeeze(tailCurv(:,1,:)),2);
% tc(:,1,:) = x;
% y = gradient(squeeze(tailCurv(:,2,:)),2);
% tc(:,2,:) = y;


disp('Calculating tail angles...')
% tailCurv = tailCurv_uncorrected;
% tA = GetTailTangents(tailCurv);
if exist('dsVecs')==1
    tA = GetTailTangents(tailCurv,[],dsVecs);
else
    tA = GetTailTangents(tailCurv);
end

tA_trl = reshape(tA,size(tA,1),nFramesInTrl,nTrls);
time_trl = reshape(time,nFramesInTrl,nTrls);

tA_5 = GetTailTangents(tailCurv,5);
curv = tA_5(end,:)';
% curv = (tA_5(end,:)-tA_5(1,:))';
% curv = max(tA_5)';
curv_trl = reshape(curv,nFramesInTrl,nTrls);

toc


%% Head Orientation
% OrientationFromVec = @(V)(angle(V(:,1) + V(:,2)*1i));
% orVec = squeeze(diff(tailCurv([1 10],:,:),[],1))';
% orVec2 = nan(size(orVec));
% for tt = 1:length(hOr);
%     if ~isempty(hOr{tt})
%          orVec2(tt,:) = diff([hOr{tt}(1,:); hOr{tt}(end,:)],[],1);
%     end
%
% end
% or = OrientationFromVec(orVec)*180/pi;
% or2 = OrientationFromVec(orVec2)*180/pi;
% ker = gausswin(10); ker = ker/sum(ker);
% % or  = chebfilt(or,1/fps,50,'low');
% or = conv2(or(:),ker(:),'same');
% or2 = conv2(or2(:),ker(:),'same');
% % or = cumsum(gradient(or));
% or_trl = reshape(or,nFramesInTrl,nTrls);
% or_trl2 = reshape(or2,nFramesInTrl,nTrls);
% for trl = 1:size(or_trl,2)
%     or_trl(:,trl) = or_trl(:,trl)-or_trl(1,trl);
%     or_trl2(:,trl) = or_trl2(:,trl)-or_trl2(1,trl);
% end
% yShift = 90;
% count = 0;
% figure('Name','Trialized Head Orientation')
% % trlList = [2 4 5 7 8 9 11];
% trlList = [1 3 11 16 18 20]
% for trl = trlList
%     count = count + 1;
%     if mod(count,2)==0
%         clr = 'r';
%     else
%         clr = 'b';
%     end
%     plot(time_trl(:,1),or_trl(:,trl)-(yShift*(trl-1)),'color',clr)
%     hold on
% %     pause()
% end


%%
% trls = [2 5 7 9 13];  % ctrl
% trls = [3 4 6 11 12] ; % abl
trls = 1:nTrls;
xLim = [-inf inf]; % In frames
cLim = [-225 225];
yShift = 350;
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
        plot((time_trl(:,1)*1000),curv_trl(:,trl)-yTick(count),'r')
        hold on
    else
        plot((time_trl(:,1)*1000),curv_trl(:,trl)-yTick(count),'b')
        hold on
    end
end
xlim((xLim/fps)*1000)
yLim =[-yTick(end)-yShift yShift];
ylim(yLim)
plot([stimTime, stimTime],yLim,'k:')
xlabel('Time (ms)')
set(gca,'tickdir','out','ytick',flipud(-yTick),'yticklabel',fliplr(trls),'xtick',[100 500:500:time_trl(end)*1000])
box off
title(['Timeseries of total tail bend angle for different trials (yShift = ' num2str(yShift) '^o)'])