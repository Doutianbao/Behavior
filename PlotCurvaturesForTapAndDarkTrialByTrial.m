%% Inputs
fps = 500;
nFramesInTrl = 750;
preStimPer = 0.1; % In seconds
firstTapTrl = 1;
yOff = [220 100];
seg = 3; % [1 - head, 2 = tail, 3 - combined]

%% Trializing 
curv = motionInfo.curv;
nTrls = size(curv,1)/nFramesInTrl;
if mod(nTrls,1) ~= 0
    error('Check value of nFramesInTrl, does not evenly divide into total # of frames')
else
    disp([num2str(nTrls) ' trials detected'])
end

curv_trl = permute(reshape(curv,[nFramesInTrl nTrls, size(curv,2)]),[1 3 2]);

time_trl = ((0:nFramesInTrl-1)/fps)-0.1;

%% Separating tap and dark flash trls
tap = curv_trl(:,:,firstTapTrl:2:end);
dark = curv_trl(:,:,firstTapTrl+1:2:end);
curv_trl = cat(4,tap,dark);

%% Plotting tap and dark trls
postStimPer = (nFramesInTrl*(1/500))-0.1;
for stim = 1:size(curv_trl,4)
    figure('Name',['Stim ' num2str(stim)])
    count = 0;
    for trl = 1:size(curv_trl,3)
        yShift = (trl-1)*yOff(stim);
        if mod(count,2)==0
            clr = 'r'
        else
            clr = 'g'
        end
        plot(time_trl*1000,curv_trl(:,seg,trl,stim)-yShift, 'color',clr)
        hold on
        drawnow
        count = count + 1;
    end    
    box off
    set(gca,'tickdir','out','color','k')
    xlabel('Time (ms)')
    ylim([-yShift-yOff(stim) yOff(stim)])
    xlim([-preStimPer*1000 postStimPer*1000])
end