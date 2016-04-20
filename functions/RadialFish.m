
function varargout = RadialFish(im,varargin)
% RadialFish - Rotates a line of certain length around fish head centroid
%   and reconstructs fish pixels that way
% radialImg = RadialFish(im)
% radialImg = RadialFish(im,fishPos)
% radialImg = RadialFish(...,dTh)
% radialImg = RadialFish(...,lineLen)
% [radialImg,lineIndMat] = ...
% [..., lineThetas] = ....
% Inputs:
% im - Image containing fish
% fishPos - x,y coordinates of head centroid. If not specified or [], then
%   automatically determines fish pos
% dTh - Angles by which to rotate line (default = 5)
% lineLen  - Length of line to rotate (default = 70)
% Outputs:
% radialImg - Radial img constructed from rotating line
% lineIndMat - Pixel indices corresponding to radial img
% lineThetas - Thetas by which line is rotated


dTh = 5;
lineLen = 70;
if nargin ==1
    fishPos = GetFishPos(im,30);
elseif nargin == 2
    fishPos = varargin{1};    
elseif nargin ==3
    fishPos = varargin{1};
    dTh = varargin{2};
elseif nargin ==4
    fishPos = varargin{1};
    dTh = varargin{2};
    lineLen = varargin{3};
end

if isempty(fishPos)
    fishPos = GetFishPos(im,30);
end
if isempty(dTh)
    dTh  = 5;
end

x = fishPos(:,1);
y = fishPos(:,2);
lineThetas = (0:dTh:360)*pi/180;
th = {};
cartCoords = {};
orientation = zeros(size(x));
orLine = {};
orLine{1} = 0;

thetaInds = 1:numel(lineThetas);
lineProf = {};
radialImg = zeros(numel(thetaInds),lineLen);
lineIndMat = radialImg;
for jj = thetaInds(:)'
    th{jj} = repmat(lineThetas(jj),1,lineLen);
    [blahX,blahY] = pol2cart(th{jj},0:lineLen-1);
    blahX = blahX + x;
    yy = y;  
    blahY = blahY  + yy;
    
    outInds = find(blahX > size(im,2));
    blahX(outInds) = 1:numel(outInds);
    %     blahX(blahX>size(im,2)) = size(im,2);
    
    outInds = find(blahX < 1);
    blahX(outInds) = size(im,2)-numel(outInds)+1:size(im,2);    
%     blahX(blahX<1) = 1;
    
    outInds = find(blahY > size(im,1));
    blahY(outInds) = 1:numel(outInds);
%     blahY(blahY>size(im,1)) = size(im,1);
    
    outInds = find(blahY < 1);
    blahY(outInds) = size(im,1)-numel(outInds)+1:size(im,1);   
%     blahY(blahY<1) = 1;

    lineInds = sub2ind(size(im),round(blahY),round(blahX));     
    lineProf{jj} = im(lineInds);
    radialImg(jj,:) = lineProf{jj};
    lineIndMat(jj,:) = lineInds(:)';
end

varargout{1} = radialImg;
varargout{2} = lineIndMat;
varargout{3} = round(lineThetas*180/pi);

% lineSums(1) = [];
% lineRanges(1) = [];
% lineSums = (lineSums-min(lineSums))/(max(lineSums)-min(lineSums));
% lineRanges = (lineRanges-min(lineRanges))/(max(lineRanges)-min(lineRanges));
% trunkLine = lineSums-lineRanges;
% trunkInd = find(trunkLine==max(trunkLine));
% orientation(frame) = lineThetas(trunkInd(1))*180/pi;
% orientation(frame) = mod(orientation(frame) + 180,360);






end