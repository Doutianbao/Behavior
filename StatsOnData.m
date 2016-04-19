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
% sideStems = {'Left','Right'};
sideStems = {'Ctrl'};
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

function data_flipped = FlipOnYForLeft(data)
fldNames = fieldnames(data);
for fldNum = 1:length(fldNames)
    fldName = fldNames{fldNum};
    if ~isempty(strfind(lower(fldName),lower('left')))
            blah =data.(fldName).motionInfo.traj_angle;
        data.(fldName).motionInfo.traj_angle = ...
            ReflectAngleOnY(blah);
        disp(['Flipped ' fldName])
    end
end
data_flipped = data;
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
    prob_left = prob_left/sum(prob_left);
    
    [prob_right,~] = hist(speed(rightInds),binVec);
    prob_right = prob_right/sum(prob_right);
    
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

function PlotGroupSpeedHists(data)
prob_contra = data.hist.traj_speed_prob_contra;
vals = data.hist.traj_speed_vals;
prob_ipsi = data.hist.traj_speed_prob_ipsi;

prob_contra_mu = mean(prob_contra,1);
prob_contra_sig = std(prob_contra,[],1);

prob_ipsi_mu = mean(prob_ipsi,1);
prob_ipsi_sig = std(prob_ipsi,[],1);

s{1} = [prob_contra_mu - prob_contra_sig, ...
    fliplr(prob_contra_mu + prob_contra_sig)];

s{2} = [prob_ipsi_mu - prob_ipsi_sig, ...
    fliplr(prob_ipsi_mu + prob_ipsi_sig)];
x =  [vals, fliplr(vals)];

figure('Name','Ipsi and Contra Swim Speed Hist')
fh{1} = fill(x,s{1},'b');
set(fh{1},'FaceAlpha',0.4)
hold on
fh{2} = fill(x,s{2},'r');
set(fh{2},'FaceAlpha',0.4)
plot(vals, prob_contra_mu,'k','linewidth',2)
plot(vals, prob_ipsi_mu,'k--','linewidth',2)
legend('Contra +/- std','Ipsi +/ std','Contra mean','Ipsi mean')
xlim([-inf 25])
ylim([-inf inf])
box off
set(gca,'tickdir','out')
xlabel('Swim speed (pxls/frame)')
ylabel('Probability')
title('Average swim speed for 1st 3 frames towards contra and ipsilateral to the side of ablation')

end
