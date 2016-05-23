
function varargout = GetArenaEdge(varargin)
%GetArenaEdge Given a reference image, returns the indices of a circle
%   that has been fit to the edge of the fish arena
% edgeInds = GetArenaEdge(refImg)
% edgeInds = GetArenaEdge(refImg,nIter)
% Inputs: 
% refImg - Reference image in which to find the arena edge (usually a single image that is the mean 
%   of all the images in the image stack)
% nIter - Max number of iterations to carry out in order fit the circle to
%   the arena edge
% Outputs:
% edgeInds - T x 2 matrix where 1st and 2nd col are the x, and y
%   coordinates of the perimeter of a circle fit to the arena edge
% 
% Avinash Pujala, HHMI, 2016

nIter = 10;
if nargin == 1
    refImg = varargin{1};
elseif nargin > 1
    refImg = varargin{1};
    nIter = varargin{2};
end

Standardize = @(x)(x-min(x(:)))/(max(x(:))-min(x(:)));
[G,~]  = imgradient(refImg);
G = Standardize(G);
E = zeros(size(G));
E(G>0.8)=1;
[eX,eY] = find(E);

edgeInds = FitCircle([eX, eY],nIter);

varargout{1} = fliplr(edgeInds);

end

