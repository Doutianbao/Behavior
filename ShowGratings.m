function varargout = ShowGratings(varargin)
%ShowGratings Displays square- or sine wave gratings based on the input
%   parameters
% parameterVec = ShowGratings(gratingDims,gratingFrequencies,gratingAngles,gratingSpeeds, gratingDurations, gratingType);
%  Inputs:
% gratingDims - Size of grating matrix
% gratingFrequencies - Vector of spatial frequencies at which to generate
%   grating
% gratingAngles - Vector of angles at which to generate grating (degree units)
% gratingSpeeds - Vector os speeds at which to move grating (given in units of phase offset in degrees)
% gratingDurations - Vector of durations (number of frames) at which to
%   display each grating
% gratingType - 'square' or 'sine';
% Outputs:
% pVec - Parameter vec that shows the timecourse of the different grating
%   parameters. T X 3 matrix, where T = length(gratingFrequencies)*
%   length(gratingAngles)*length(gratingSpeeds)*D. The 1st, 2nd, and 2rd
%   col are the timecourses of the grating frequency, angle, and speed
%   respectively.
%
% Avinash Pujala, HHMI, 2016

gDims = [600, 600]; % Grating dimensions
F = 10; % Grating spatial frequency
A = 0; % Grating angle in degrees
V = 15; % Grating velocity (Inverse of pause duration before showing gratings circshifted by 1 pixel)
D = 10; % Duration in seconds
gType = 'square'; % Grating type
testMode = 0;

if nargin ==1
    gDims = varargin{1};
elseif nargin == 2
    gDims = varargin{1};
    F = varargin{2};
elseif nargin  == 3;
    gDims = varargin{1};
    F = varargin{2};
    A = varargin{3};
elseif nargin  ==4
    gDims = varargin{1};
    F = varargin{2};
    A = varargin{3};
    V = varargin{4};
elseif nargin ==5
    gDims = varargin{1};
    F = varargin{2};
    A = varargin{3};
    V = varargin{4};
    D = varargin{5};
elseif nargin ==6
    gDims = varargin{1};
    F = varargin{2};
    A = varargin{3};
    V = varargin{4};
    D = varargin{5};
    gType = varargin{6};
elseif nargin == 7
    gDims = varargin{1};
    F = varargin{2};
    A = varargin{3};
    V = varargin{4};
    D = varargin{5};
    gType = varargin{6};
    testMode = varargin{7};
elseif nargin > 7
    error('Too many inputs!')
end

if isempty(gDims)
    gDims =[600, 600];
end

if isempty(gType)
    gType = 'square';
end

if isempty(testMode)
    testMode = 0;
end
if mod(D,1)~=0
    error('Duration must be specified as number of frames, and must be an integer!')
end

V = round(V*(gDims(1)/100));
count = 0;
pVec = zeros(length(F)*length(A)*length(V)*D,3);
fh = figure;
set(fh,'units','normalized','menubar','none','NumberTitle','off');
if testMode == 0;
    set(fh,'position',[0.9 0 1 1])
end
for f = F(:)'
    for a = A(:)'
        for v = round(V(:))'
            for d = 1:round(D);           
            G = CreateGrating(gDims,f,0,a,gType);
            G = circshift(G,[count 0]);
            count = count + 1;
            pVec(count,1) = f;
            pVec(count,2)  = a;
            pVec(count,3)= v;
            cla
            imagesc(G),colormap(gray), %axis image
            axis off
            %                 title(['Frame # ' num2str(count) ', freq = ' num2str(f) ', angle = ' num2str(a) ', vel = ' num2str(v)])
            drawnow
            shg
            pause(1/v)
            end
        end
    end
end

varargout{1} = pVec;
end

function G = CreateGrating(varargin)
% CreateGrating - Create a square- or sinewave grating based on input
%   parameters
% G = CreateGrating(gratingDims,spatialFrequency,phaseOffset,gratingAngle,gratingType)
% Outputs:
% G - Grating with defined input parameters

gDims = varargin{1}; % Dimensions of the grating
f = varargin{2}; % Spatial frequency of the grating
phi = varargin{3}; % Phase offset of the grating
gAngle = varargin{4}; % Angle of the grating
gType = varargin{5}; % Type of grating, square or sine
phi = phi*pi/180; % Convert to radians

lenDiag = ceil(sqrt(sum(gDims.^2)));
x = linspace(0,1,lenDiag);
if strcmpi(gType,'square')
    y = square(2*pi*f*x + phi);
elseif strcmpi(gType,'sine')
    y = sin(2*pi*f*x + phi);
end

[~,G] = meshgrid(x,y);
G = imrotate(G,gAngle,'crop');
xInds = ceil([(size(G,1)-gDims(1))/2, gDims(1)+(size(G,1)-gDims(1))/2]);
yInds = ceil([(size(G,2)-gDims(2))/2, gDims(2)+(size(G,2)-gDims(2))/2]);
G = G(xInds(1):xInds(2)-1,:);
G = G(:,yInds(1):yInds(2)-1);


end