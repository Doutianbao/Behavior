function data = StatsOnData(dataDir,varargin)
%StatsOnData Reads relevant fish data files(.mat) in the input dir
%  as long as the files are labeled in a certain way.
% data = StatsOnData(dataDir, sideStems);
% Inputs:
% dataDir - Directory where group data is stored (requires fish pos, orientation,
%   and ref image)
% sideStems - File stem which indicates if ablation was to the left, right
%   or if the data is from an control unablated fish; {'Left','Right'}, or
%   {'Ctrl'}

clear data
sideStems = {'Left','Right'};
% sideStems = {'Ctrl'};
if nargin ==2
    sideStems = varargin{1};
end
data = GetSideData(dataDir,sideStems);

disp('Computing and appending motion info...')

data = AppendMotionInfo(data);


disp('Flipping trajectory angles on the Y-axis for left ablations...')
data = FlipOnYForLeft(data);

disp('Appending trajectory angle histogram info...')
[data,fldNames_orig] = AppendTrajAngleHists(data);

disp('Appending traj speed histogram info...')
data = AppendTrajSpeedHists(data,fldNames_orig);

disp('Appending 1st order turn angle histogram info...')
data = Append1stOrderTurnAngleHist(data,fldNames_orig);

disp('Appending group curvature info...')
data = AppendGroupCurv(data,fldNames_orig);

disp('Appending motion vec info...')
data = AppendMotionVecInfo(data,fldNames_orig);

%# Saving data
saveOrNot = input('Save data (y/n)?  :', 's');
if strcmpi(saveOrNot,'y')
    procDir = fullfile(dataDir,'Proc');
    if exist(procDir)~=7
        mkdir(procDir)
    end
    save(fullfile(procDir,['groupData_' [sideStems{:}] '.mat'] ),'data');
    disp(['Saved to ' procDir])
end

%# Plot traj angle histograms for group data
PlotGroupAngleHists(data)

%# Plot swim velocity histograms for group data
PlotGroupSpeedHists(data)

%# Plot 1st order turn angle histograms
Plot1stOrderTurnAngleHists(data)

% PlotTrajAngles(data,fldNames_orig)

%# Plot head curvature histograms
PlotCurvHist(data)

%# Plot motion vec histograms
PlotMotionVecHists(data)
end

function sideData = GetSideData(dataDir,sideStems)
sideData = struct;
fileNames = sort(GetFileNames(dataDir,'.mat'));
for side = 1:length(sideStems)
    disp(['Reading ' sideStems{side} ' side...'])
    sideLen = length(sideStems{side});
    for fn = 1:length(fileNames)
        sInd = strfind(lower(fileNames{fn}),lower(sideStems{side}));
        if ~isempty(sInd)
            disp(['Reading file ' fileNames{fn}])
            fldName = lower(fileNames{fn}(sInd:sInd+ sideLen));
            blah = load(fullfile(dataDir,fileNames{fn}));
            fldName_sub = fieldnames(blah);
            if strcmpi(fldName_sub{1},'tracexy');
                fldName_sub_mod = 'fishpos';
            else
                fldName_sub_mod = lower(fldName_sub{1});
            end
            sideData.(fldName).(fldName_sub_mod) = blah.(fldName_sub{1});
        end
    end
end
end

function fileNames = GetFileNames(fileDir,fileExt)
if nargin < 2
    error('At least 2 inputs required!')
end
filesInDir = dir(fileDir);
filesInDir= {filesInDir.name};
remInds = [];
for fn = 1:length(filesInDir)
    if isempty(strfind(filesInDir{fn},fileExt))
        remInds = [remInds; fn];
    end
end
filesInDir(remInds) = [];
fileNames = filesInDir;
end

function data_motion = AppendMotionInfo(data)
fldNames = fieldnames(data);
for fldNum = 1:length(fldNames)
    fldName = fldNames{fldNum};
    var = data.(fldName);
    imgDims = size(var.ref);
    motionInfo = GetMotionInfo(var.fishpos,var.orientation,imgDims(1));
    data.(fldName).motionInfo = motionInfo;
    data_motion = data;
end

end

function [data,fldNames] = AppendTrajAngleHists(data)
binVec = union(0:10:80, 100:10:360);
fldNames = fieldnames(data);
tah = nan(length(fldNames),length(binVec)+1);
bv = nan(1,length(binVec)+1);
nTrls = nan(length(fldNames),1);
for fldNum = 1:length(fldNames)
    fldName = fldNames{fldNum};
    angles = data.(fldName).motionInfo.traj_angle;
    nTrls(fldNum) = numel(angles);
    [rho_prob,theta_prob] = hist(angles,binVec);
    rho_prob = [rho_prob(:); rho_prob(1)];
    theta_prob = [theta_prob(:); theta_prob(1)];
    rho_prob = rho_prob/sum(rho_prob);
    tah(fldNum,:) = rho_prob;
end
bv(1,:) = theta_prob;
data.hist.traj_angle_prob = tah;
data.hist.traj_angle_vals = bv;
data.hist.traj_angle_nTrls = nTrls;
end

function data = AppendTrajSpeedHists(data,fldNames)
binVec = 1:0.75:30;
prob_ipsi = nan(length(fldNames),length(binVec));
vals = nan(1,length(binVec));
prob_contra = prob_ipsi;
for fldNum  = 1:length(fldNames)
    fldName = fldNames{fldNum};
    speed = data.(fldName).motionInfo.traj_speed;
    angles = data.(fldName).motionInfo.traj_angle-90;
    nanInds = find(isnan(speed));
    zerInds = find(speed==0);
    remInds = union(nanInds,zerInds);
    speed(remInds)=[];
    angles(remInds) = [];
    leftInds = find(angles > 5);
    rightInds = find(angles <5);
    [prob_left,vals] = hist(speed(leftInds),binVec);
%     prob_left = prob_left/sum(prob_left);
    
    [prob_right,~] = hist(speed(rightInds),binVec);
%     prob_right = prob_right/sum(prob_right);

%     N = sum(prob_right) + sum(prob_left);
%    prob_left = prob_left/N;
%    prob_right = prob_right/N;
    if ~isempty(strfind(lower(fldName),'left'))
        prob_contra(fldNum,:) = prob_right;
        prob_ipsi(fldNum,:) = prob_left;
    elseif ~isempty(strfind(lower(fldName),'right'))
        prob_contra(fldNum,:) = prob_left;
        prob_ipsi(fldNum,:) = prob_right;
    else
        toss = rand(1);
        if toss < 0.5
            prob_contra(fldNum,:) = prob_right;
            prob_ipsi(fldNum,:) = prob_left;
        else
            prob_contra(fldNum,:) = prob_right;
            prob_ipsi(fldNum,:) = prob_left;
        end
    end
end
data.hist.traj_speed_prob_ipsi = prob_ipsi;
data.hist.traj_speed_prob_contra = prob_contra;
data.hist.traj_speed_vals(1,:) = vals;
end

function data = Append1stOrderTurnAngleHist(data,fldNames)
binVec = 3:3:120;
prob_ipsi = nan(length(fldNames),length(binVec));
vals = nan(1,length(binVec));
prob_contra = prob_ipsi;
for fldNum  = 1:length(fldNames)
    fldName = fldNames{fldNum};
    dOr = data.(fldName).motionInfo.dOr;
    dOr(abs(dOr)<4)=[]; dOr(abs(dOr)>120)=[];
    y = -dOr(dOr<0);
    [blah,~] = hist(y,binVec); % Sign already flipped based on 'left' or 'right' in filename
%     blah = blah/sum(blah);
    prob_contra(fldNum,:) = blah;
    [blah,vals] = hist(dOr(dOr>0),binVec);
%     blah = blah/sum(blah);
    prob_ipsi(fldNum,:) = blah;
%     N = sum(prob_contra(fldNum,:)) + sum(prob_ipsi(fldNum,:));
%     prob_contra(fldNum,:) = prob_contra(fldNum,:)/N;
%     prob_ipsi(fldNum,:) = prob_ipsi(fldNum,:)/N;
end
data.hist.turnAngle_prob_ipsi = prob_ipsi;
data.hist.turnAngle_prob_contra = prob_contra;
data.hist.turnAngle_vals(1,:) = vals;

end

function data = AppendGroupCurv(data,fldNames)
binVec = 3:3:120;
prob_ipsi = nan(length(fldNames),length(binVec));
vals = nan(1,length(binVec));
prob_contra = prob_ipsi;
for fldNum  = 1:length(fldNames)
    fldName = fldNames{fldNum};
    curv = data.(fldName).motionInfo.curv;
    if strfind(lower(fldName),'left')
        curv = -curv;
    end
    curv(abs(curv)<4)=[]; curv(abs(curv)>120)=[];
    y = -curv(curv<0);
    [blah,~] = hist(y,binVec); % Sign already flipped based on 'left' or 'right' in filename
%     blah = blah/sum(blah);
    prob_contra(fldNum,:) = blah;
    [blah,vals] = hist(curv(curv>0),binVec);
%     blah = blah/sum(blah);
    prob_ipsi(fldNum,:) = blah;
%     N = sum(prob_contra(fldNum,:)) + sum(prob_ipsi(fldNum,:));
%     prob_contra(fldNum,:) = prob_contra(fldNum,:)/N;
%     prob_ipsi(fldNum,:) = prob_ipsi(fldNum,:)/N;
end
data.hist.curv_prob_ipsi = prob_ipsi;
data.hist.curv_prob_contra = prob_contra;
data.hist.curv_vals(1,:) = vals;

end

function data = AppendMotionVecInfo(data,fldNames)
motionThr = 5;
binVec1 = 2:100;  
prob_ipsi_mag = nan(length(fldNames),length(binVec1));
vals_mag = nan(1,length(binVec1));
prob_contra_mag = prob_ipsi_mag;

binVec2 = 4:4:120;    
prob_ipsi_ang = nan(length(fldNames),length(binVec2));
prob_contra_ang = prob_ipsi_ang;
vals_ang =nan(1,length(binVec2));

for fldNum  = 1:length(fldNames)
    fldName = fldNames{fldNum};
    mv = data.(fldName).motionInfo.motionVecs;
    %## Remove straight and noise trajectories
    remInds = find(mv(:,2) <=motionThr | mv(:,2) > 100); % No motion and artifactual jump frames
    mv(remInds,:) = [];
    remInds = find(abs(mv(:,1))<=4 | abs(mv(:,1))>=120); % Straight movt and artifactual turn frames
    mv(remInds,:) = [];
    
    if strfind(lower(fldName),'left')
        mv(:,1) = -mv(:,1);
    end
    %## Magnitudes      
    [blah,~] = hist(mv(mv(:,1)<0,2),binVec1); 
%     blah = blah/sum(blah);
    prob_contra_mag(fldNum,:) = blah;
    [blah,vals_mag] = hist(mv(mv(:,1)>0,2),binVec1);
%     blah = blah/sum(blah);
    prob_ipsi_mag(fldNum,:) = blah;
    
%     N = sum(prob_contra_mag(fldNum,:)) + sum(prob_ipsi_mag(fldNum,:)); 
%     prob_contra_mag(fldNum,:) = prob_contra_mag(fldNum,:)/N;
%     prob_ipsi_mag(fldNum,:) = prob_ipsi_mag(fldNum,:)/N;
    
    %## Angles
    [blah,~] = hist(-mv(mv(:,1)<0,1),binVec2); 
%     blah = blah/sum(blah);
    prob_contra_ang(fldNum,:) = blah;
    [blah,vals_ang] = hist(mv(mv(:,1)>0),binVec2);
%     blah = blah/sum(blah);
    prob_ipsi_ang(fldNum,:) = blah;
    
%     N  = sum(prob_contra_ang(fldNum,:)) + sum(prob_ipsi_ang(fldNum,:));
%     prob_contra_ang(fldNum,:) = prob_contra_ang(fldNum,:)/N;
%     prob_ipsi_ang(fldNum,:) = prob_ipsi_ang(fldNum,:)/N;
end
data.hist.motionVecMag_prob_ipsi = prob_ipsi_mag;
data.hist.motionVecMag_prob_contra = prob_contra_mag;
data.hist.motionVecMag_vals(1,:) = vals_mag;

data.hist.motionVecAng_prob_ipsi = prob_ipsi_ang;
data.hist.motionVecAng_prob_contra = prob_contra_ang;
data.hist.motionVecAng_vals(1,:) = vals_ang;
end

function data_flipped = FlipOnYForLeft(data)
fldNames = fieldnames(data);
for fldNum = 1:length(fldNames)
    fldName = fldNames{fldNum};
    if ~isempty(strfind(lower(fldName),lower('left')))
        blah =data.(fldName).motionInfo.traj_angle;
        data.(fldName).motionInfo.traj_angle = ...
            ReflectAngleOnY(blah);
        data.(fldName).motionInfo.dOr = -data.(fldName).motionInfo.dOr;
        disp(['Flipped ' fldName])
    end
end
data_flipped = data;
end

function PlotGroupAngleHists(data)
var = data.hist.traj_angle_prob;
muVec = mean(var,1);
sigVec = std(var,[],1);
thetas = data.hist.traj_angle_vals*pi/180;

figure('Name','Group traj angle histogram_mu and sigma')
ph = polar(thetas,muVec+sigVec,'r--');
set(ph,'linewidth',2)
hold on
ph = polar(thetas,muVec-sigVec,'r--');
set(ph,'linewidth',2)
ph = polar(thetas,muVec,'k');
set(ph,'linewidth',2)
title([num2str(size(var,1)) ' fish, ' num2str(sum(data.hist.traj_angle_nTrls)) ' swim episodes'])

figure('Name','Group traj angle histogram_all hists')
mRows = max(var,[],2);
[~,mRows] = sort(mRows,'descend');
for jj = mRows(:)'
    polar(thetas,var(jj,:))
    hold on
end
title([num2str(size(var,1)) ' fish, ' num2str(sum(data.hist.traj_angle_nTrls)) ' swim episodes'])

end

function PlotGroupSpeedHists(data)
prob_contra = data.hist.traj_speed_prob_contra;
vals = data.hist.traj_speed_vals;
prob_ipsi = data.hist.traj_speed_prob_ipsi;

xLabel = 'Swim speed';
yLabel = 'Prob';
ttl = 'Swim speed histogram';
xLim = [0 25];
PlotGroupData(prob_contra, prob_ipsi, vals,xLabel,yLabel,ttl, xLim)
PlotFishData(prob_contra,prob_ipsi,vals, xLabel, yLabel, ttl, xLim)
end

function Plot1stOrderTurnAngleHists(data)
prob_contra = data.hist.turnAngle_prob_contra;
vals = data.hist.turnAngle_vals;
prob_ipsi = data.hist.turnAngle_prob_ipsi;

xLabel = '1st order turn angle';
yLabel = 'Prob';
ttl = '1st order turn angle hist';
xLim = [0 70];
PlotGroupData(prob_contra, prob_ipsi, vals,xLabel,yLabel,ttl, xLim)
PlotFishData(prob_contra,prob_ipsi,vals, xLabel, yLabel, ttl, xLim)

end

function PlotTrajAngles(data,fldNames)
for fldNum = 1:length(fldNames)
    fldName = fldNames{fldNum};
    figure
    var =data.(fldName).motionInfo.traj_adj;
    for jj = 2:length(var)
        x = var{jj}(:,1);
        y = var{jj}(:,2);
        if length(x)>3
            a = angle(x(3) + y(3)*1i)*180/pi;
            if abs(a)< 180;
                patchline(x,y,'facealpha',0.1);
                %             plot(x,y)
                drawnow
                axis image
                hold on
            end
        end
        
    end
    title('Trajectories')
end

end

function PlotCurvHist(data)
prob_contra = data.hist.curv_prob_contra;
vals = data.hist.curv_vals;
prob_ipsi = data.hist.curv_prob_ipsi;

xLabel = 'Heading angles';
yLabel = 'Prob';
ttl = 'Head curvature histogram';
xLim = [0 70];
PlotGroupData(prob_contra, prob_ipsi, vals,xLabel,yLabel,ttl, xLim)
PlotFishData(prob_contra,prob_ipsi,vals, xLabel, yLabel, ttl, xLim)

end

function PlotMotionVecHists(data)
prob_contra = data.hist.motionVecAng_prob_contra;
vals = data.hist.motionVecAng_vals;
prob_ipsi = data.hist.motionVecAng_prob_ipsi;

xLabel = 'Motion vec angle';
yLabel = 'Prob';
ttl = 'Histogram of motion vec angles';
xLim = [0 70];
PlotGroupData(prob_contra, prob_ipsi, vals,xLabel,yLabel,ttl, xLim)
PlotFishData(prob_contra,prob_ipsi,vals, xLabel, yLabel, ttl, xLim)


prob_contra = data.hist.motionVecMag_prob_contra;
vals = data.hist.motionVecMag_vals;
prob_ipsi = data.hist.motionVecMag_prob_ipsi;

xLabel = 'Motion vec length';
yLabel = 'Prob';
ttl = 'Histogram of motion vec lengths';
xLim = [0 50];

PlotGroupData(prob_contra, prob_ipsi, vals,xLabel,yLabel,ttl, xLim)
PlotFishData(prob_contra,prob_ipsi,vals, xLabel, yLabel, ttl, xLim)

end


    function PlotGroupData(prob_contra, prob_ipsi,vals,xLabel,yLabel,ttl,xLim) 
        N_vec = sum(prob_ipsi,2) + sum(prob_contra,2);
        N = repmat(N_vec,1,size(prob_ipsi,2));
        prob_ipsi = prob_ipsi./N;
        prob_contra = prob_contra./N;
        
        
        prob_contra_mu = mean(prob_contra,1);
        prob_contra_sig = std(prob_contra,[],1);
        
        prob_ipsi_mu = mean(prob_ipsi,1);
        prob_ipsi_sig = std(prob_ipsi,[],1);
        
        s{1} = [prob_contra_mu - prob_contra_sig, ...
            fliplr(prob_contra_mu + prob_contra_sig)];
        
        s{2} = [prob_ipsi_mu - prob_ipsi_sig, ...
            fliplr(prob_ipsi_mu + prob_ipsi_sig)];
        x =  [vals, fliplr(vals)];
        
        figure
        fh{1} = fill(x,s{1},'b');
        set(fh{1},'FaceAlpha',0.4)
        hold on
        fh{2} = fill(x,s{2},'r');
        set(fh{2},'FaceAlpha',0.4)
        plot(vals, prob_contra_mu,'k','linewidth',2)
        plot(vals, prob_ipsi_mu,'k--','linewidth',2)
        
        y = max([prob_contra_mu(:); prob_ipsi_mu(:)]);
        y = y(1)*0.9;
        x = ceil(max(vals)/3);
        nTrls = round(2*sum(vals)*size(prob_contra,1));
        text(x,y,num2str(nTrls))
        
        legend('Contra +/- std','Ipsi +/ std','Contra mean','Ipsi mean')
        xlim(xLim)
        ylim([-inf inf])
        box off
        set(gca,'tickdir','out')
        xlabel(xLabel)
        ylabel(yLabel)
        title(ttl)
        
    end

    function PlotFishData(prob_contra,prob_ipsi,vals,xLabel,yLabel,ttl, xLim)
        N_vec = sum(prob_ipsi,2) + sum(prob_contra,2);
        N = repmat(N_vec,1,size(prob_ipsi,2));
        prob_ipsi = prob_ipsi./N;
        prob_contra = prob_contra./N;
        plotAllFish = input('Plot motion vecs for all fish? (y/n): ','s');
        if strcmpi(plotAllFish,'y')
            for jj= 1:size(prob_contra,1)
                figure('Name',['Fish # ' num2str(jj)])
                title(['Fish # ', num2str(jj)])
                plot(vals,prob_contra(jj,:),'b.-')
                hold on;
                plot(vals,prob_ipsi(jj,:),'r.-')
                legend('Contra', 'Ipsi')
                y = max([prob_contra(jj,:), prob_ipsi(jj,:)]);
                y = y(1);
                x = ceil(max(vals)/3);
                nTrls = N_vec(jj);
                text(x,y,num2str(nTrls))
                xlim(xLim)
                xlabel(xLabel)
                ylabel(yLabel)
                title(ttl)
            end
        end
    end


