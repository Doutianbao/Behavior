function PlayAnnotatedFishFrames(IM,varargin)
% PlayAnnotatedFishFrames - Plays fish frames with annotation
% PlayAnnotatedFishFrames(IM, 'frames',frames,'fishPos',fishPos,'fishOrientation',fishOrientation,'pauseDur',[])
% Inputs:
% IM - Image stack of size M x N x T, where M = image height, N = image
%   width, T = # of time points
% fishPosVec - T x 2 vec where 1st & 2nd cols contain row & col positions
%   of the fish in the images
% startFrame - Frame from which to display video
% endFrame - Last Frame of video display
% pauseDur - Duration of pause between subsequent frames

frames = [];
fishpos = [];
fishorientation = [];
pausedur = 0;
pathlength = 500;

nArgs = nargin;
args = varargin;
% str = cellfun(@isstr,args);
% strInds = find(str);
% args_str = args(strInds);


for jj = 1:2:length(args)
    varName = lower(args{jj});
    eval([varName '= args{jj+1};']);
end

if isnan(pausedur)
    pausedur = 0.1;
end

if isempty(frames)
    frames = 1:size(IM,3);
end
turnInds = getTurnInds(fishpos,frames);

figure('Name', 'Annotated Fish Tracking')
tic
count = 0;
pos = [nan nan];
frames(frames==1)=[];
for imgNum = frames(:)'
    count = count + 1;
    cla
    imagesc(IM(:,:,imgNum)),axis image, axis off, colormap(gray), drawnow
    hold on
    pos = [pos; [fishpos(imgNum,1), fishpos(imgNum,2)]];
    if mod(count,pathlength)==0
        pos = pos(end-2:end,:);
    end
    plot(pos(:,1),pos(:,2),'r.-'), drawnow
    
    if count > 1
        prevFrame = frames(count);
    else
        prevFrame = nan;
    end
    
    turnDir = 'N/A';
    if ~isempty(intersect(prevFrame,turnInds{1}))
        turnDir =  'Left';
    elseif ~isempty(intersect(prevFrame,turnInds{2}))
        turnDir = 'Right';
    elseif ~isempty(intersect(prevFrame,turnInds{3}))
        turnDir = 'Straight';
    end
    
    title(['Frame: ' num2str(imgNum) ', Turn: ' turnDir ', Angle: ' num2str(abs(round(turnInds{4}(count))))])
    shg
    if isempty(pausedur)
        pause()
    else
        pause(pausedur)
    end
end
end


function turnInds = getTurnInds(fishPos,motionFrames)
turnAngleThr = 5;

fishPos_motion = fishPos(motionFrames,:);
motionVecs = diff(fishPos_motion);
motionVecs = motionVecs(:,1) + motionVecs(:,2)*1i;
dTh = angle(motionVecs(1:end-1).*conj(motionVecs(2:end)));
dTh_motion = [0; dTh(:)]*180/pi;

rightInds = motionFrames(find(dTh_motion < -turnAngleThr)+1);
leftInds = motionFrames(find(dTh_motion > turnAngleThr)+1);
straightInds = motionFrames(find(abs(dTh_motion) <=5)+1);

turnInds = {};
turnInds{1} = leftInds;
turnInds{2} = rightInds;
turnInds{3} = straightInds;
turnInds{4} = dTh_motion;

end

