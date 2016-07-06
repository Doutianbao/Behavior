function ipd = GetPxlSpacing(ref,varargin)
%GetPxlsSpacing Given a reference image wherein the circular edge of the fish arena
%   is clearly visible, uses the specified diameter of the arena edge to
%   determine the inter-pixel spacing of the image in the same units as the
%   specified diameter
%   First, the edge of the arena is detected and fit with a circle. The
%   dimater of this circle is then estimated and matched to the diameter
%   input by the user. This allows for conversion of distance in pixels to
%   distances in the units in which the arena diameter is entered
% 
% pxlDist = GetPxlSpacing(refImg);
% pxlDist = GetPxlSpacing(refImg, 'diam, arenaDiameter,'detThr',detThr,'nIter',nIter,'tol',tol,'plotBool',plotBool)
% 
% Inputs:
% refImg - Reference image in which to detect arena edge
% arenaDiameter - Actual arena diameter specified by user in units of
%   choice
% detThr - Threshold for detecting a few points on an the arena's edge
% nIter - Max # of iterations to step through to fit a circle to the
%   arena's edge
% tol - Tolerance for the reduction in circle fit to the edge from one iteration to the next.
%   If reduction in error of fit from one iteration to the next does not
%   falls at or below this tolerance, then does not go through all
%   iterations.
% plotBool - 0 or 1. 1 results in plotting the fit of the circle to the
%   arena's edge
% 
% Avinash Pujala, Koyama lab/HHMI, 2016

diam = 25; % This is the diameter of the arena that I am currently using in mm.
detThr = 0.7;
nIter = 100;
tol = 0.1;
plotBool = 1;

for jj = 1:numel(varargin)
    switch lower(varargin{jj}) 
        case 'diam'
            diam = varargin{jj+1};
        case lower('detThr')
            detThr = varargin{jj+1};
        case lower('nIter')
            nIter = varargin{jj+1};
        case 'tol'
            tol = varargin{jj+1};
        case lower('plotBool')
            plotBool = varargin{jj+1};
    end
end

edgeInds = GetArenaEdge(ref,'nIter',nIter,'detThr',detThr,'nIter',nIter,'tol',tol,'plotBool',plotBool);

diam_pxl  = 2* EstimateCircleDiam(edgeInds);

ipd = diam/diam_pxl;

end

function diam = EstimateCircleDiam(cInds)
% Given a set of inds which outline a cicle in 2D Cartesian space, returns the estimated
%   diameter of the circle in pxl units.
% circleDiam = EstimateCircleDiam(circleInds);
% Inputs:
% circleInds - N x 2 matrix, where the 1st and 2nd cols specify the x and y positions of the circle in Cartesian space  

ctr = mean(cInds,1);
cInds = cInds - repmat(ctr,size(cInds,1),1);

[~, rho] = cart2pol(cInds(:,1), cInds(:,2));

diam = mean(rho);

end