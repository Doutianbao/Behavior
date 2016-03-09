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
rho = rho/max(rho);
rho = [rho(:); rho(1)];
theta = [theta(:); theta(1)]; 

figure('Name','Polar histogram plot')
ph{1} = mmpolar(theta,rho,'TTickDelta',tTickDelta,'RTickValue',rTickVal,'TGridColor','r','Border','off','TTickLabel',ttl);


figure('Name','Rose histogram plot')
ph{2} = rose(angles,nBins);
varargout{1} = ph;


end

