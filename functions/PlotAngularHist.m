function varargout = PlotAngularHist(varargin)
%PlotAngleHist Given a vector of angles (in radians), plots angular 
%   histogram as polar plot and rose diagram
% 
% PlotAngleHist(angles)
% figHandle = PlotAngleHist(angles, nBins)

tTickDelta = 45;
rTickVal = [0.5 1];
ttl = {' 90^{\circ} ','0^{\circ}','-90^{\circ} ','\pm 180^{\circ} '};
angles = varargin{1};
nBins = 50;
if nargin ==2
    nBins = varargin{2};
end


[rho,theta] = hist(angles,nBins);
rho = [rho(:); rho(1)];
nTrials = sum(rho);
theta = [theta(:); theta(1)]; 
rho_max = rho/max(rho);

binVec = union(0:10:80, 100:10:360)*pi/180;
[rho_prob,theta_prob] = hist(angles,binVec);
rho_prob = [rho_prob(:); rho_prob(1)];
theta_prob = [theta_prob(:); theta_prob(1)];
rho_prob = rho_prob/sum(rho_prob);

figure('Name','Polar histogram plot')
mmpolar(theta,rho_max,'TTickDelta',tTickDelta,'RTickValue',rTickVal,'TGridColor','r','Border','off','TTickLabel',ttl);

figure('Name','Polar histogram_probability, zero offset')
mmpolar(theta_prob,rho_prob,'TTickDelta',tTickDelta,'RTickValue',rTickVal,'TGridColor','r','Border','off','TTickLabel',ttl);
title(['nTrials = ' num2str(nTrials)])

figure('Name','Rose histogram plot')
ph{1} = rose(angles,nBins);
varargout{1} = ph;

%## Zero offset histogram
figure('Name','Rose plot_zero offset')
ph{2} = rose(angles,binVec);
title(['nTrials = ' num2str(nTrials)])
end

