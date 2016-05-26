
function varargout = GetArenaEdge(varargin)
%GetArenaEdge Given a reference image, returns the indices of a circle
%   that has been fit to the edge of the fish arena
% edgeInds = GetArenaEdge(refImg)
% edgeInds = GetArenaEdge(refImg,'nIter', nIter)
% edgeInds = GetArenaEdge(refImg,nIter,'detThr', detThr);
% Inputs:
% refImg - Reference image in which to find the arena edge (usually a single image that is the mean
%   of all the images in the image stack)
% detThr - Threshold for detecting some edge points in the original image.
% nIter - Max number of iterations to carry out in order fit the circle to
%   the arena edge
% tol - Error tolerance from one iteration to the next for accepting the circle fit to edge as
%   good (Default =  0.25). If this level of tolerance is reached before the number of
%   iterations then, breaks out of the loop. 
% plotBool - 0 or 1, determines whether or not to display fit

% Outputs:
% edgeInds - T x 2 matrix where 1st and 2nd col are the x, and y
%   coordinates of the perimeter of a circle fit to the arena edge
%
% Avinash Pujala, HHMI, 2016

detThr = 0.8;
nIter = 10;
tol = 0.25;
plotBool = 1;

for jj = 2:numel(varargin)
    switch lower(varargin{jj})
        case lower('nIter')
            nIter = varargin{jj+1};
        case lower('detThr')
            detThr = varargin{jj+1};
        case 'tol'
            tol = varargin{jj+1};
        case lower('plotBool')
            plotBool = varargin{jj+1};
    end
end
refImg = varargin{1};

Standardize = @(x)(x-min(x(:)))/(max(x(:))-min(x(:)));
[G,~]  = imgradient(refImg);
G = Standardize(G);
E = zeros(size(G));
E(G>detThr)=1;
[eX,eY] = find(E);

[edgeInds,error] = FitCircle([eX, eY],nIter,tol);
edgeInds = fliplr(edgeInds);

if plotBool
    figure
    imagesc(refImg), axis image
    hold on
    plot(eY,eX,'g.')
    plot(edgeInds(:,1),edgeInds(:,2),'m.')
    legend('Detected points','Fit points','Location','best')
    title(['Circle fit to the edge of the arena, with error: ' num2str(error) '%'])
end

varargout{1} = edgeInds;

end

function varargout = FitCircle(cInds,varargin)
%FitCircle Given a set of indices of roughly circular shape, returns
%   indices of a circle that fits the shape
%
% circleInds = FitCircle(circleLikeInds,nIter,tolerance)
% [circleInds,error] = FitCircle(...)
% Inputs:
% cirleLikeInds - x,y indices of a circle like shape (such as the edge indices for an arena in which a
%    fish swims in during behavior expts)
% nIter - Number of iterations, through which the program must go in order to converge onto a circle
% tolerance - Error tolerance (default: 0.5). Iterations stop when
%   difference in error for previous and currrent iterations falls below
%   the tolerance
% Outputs:
% circleInds
%
% Avinash Pujala, Koyamalab/HHMI, 2016


tol = 0.25; % Tolerance for percentage error improvement between one iteration and the next
shift = [mean(cInds(:,1)), mean(cInds(:,2))];
shift_orig = shift;
x = cInds(:,1)-shift(1);
y = cInds(:,2)-shift(2);

if nargin ==1
    nIter = 1;
elseif nargin ==2
    nIter = varargin{1};
elseif nargin ==3;
    nIter = varargin{1};
    tol = varargin{2};
end
%## Iteratively find a circle that fits indices well

for iter = 1:nIter
    [theta,rho] = cart2pol(x,y);
    rho_fit = mean(rho)*ones(size(rho));
    error = sqrt(sum((rho-rho_fit).^2))/numel(rho);
    disp(['Iter# ' num2str(iter) ', Error = ' num2str(error*100) '%'])
    [x_fit,y_fit] = pol2cart(theta,rho_fit);
    if iter > 1 && ((error_prev-error)*100) < tol
        disp(['Error within tolerance of ' num2str(tol) ' %, quitting!'])
        break
    end
    offset = [mean(x)-mean(x_fit), mean(y)-mean(y_fit)];
    shift = shift + offset;
    x = x - offset(1);
    y = y - offset(2);
    error_prev = error;
end
%## Fill in missing points in the circle and uniformly downsample
[theta,rho] = cart2pol(x_fit,y_fit);
tt = linspace(-pi,pi,numel(theta));
rr = interp1(theta,rho,tt,'spline');
[x_fit,y_fit]= pol2cart(tt,rr);

% figure
% plot(x_fit,y_fit,'r.'), hold on
% plot(x,y,'b.'),axis image
% legend('Fit','Original','Location','best')
% xlim([-inf inf])
% title(['Iter# ' num2str(iter) ', Error = ' num2str(error*100) '%'])

varargout{1} = [x_fit(:)+shift(1), y_fit(:)+shift(2)];
varargout{2} = error;

end


