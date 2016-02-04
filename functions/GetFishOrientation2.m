function orientation = GetFishOrientation2(IM,fishPos,varargin)
%GetFishOrientation - Given an image stack gets the orientaton of the fish
%   in each image of the stack using the fish position input
% orientation = GetFishOrientation(IM,fishPos)
% orientation  = GetFishOrientatin(IM,fishPos,lineLength)

nWorkers = 10;
dTheta = 10;
if nargin < 3
    lineLength = 25;
elseif nargin == 3
    lineLength = varargin{1};
end

x = fishPos(:,1);
y = fishPos(:,2);
lineThetas = (0:dTheta:360)*pi/180;
th = {};
cartCoords = {};
imgFrames = 2:size(IM,3);
% imgFrames = 1000:1010;
orientation = zeros(size(x));
orLine = {};
orLine{1} = 0;

if 0
    figure
blankImg = zeros(size(IM(:,:,2)));
end
thetaInds = 1:numel(lineThetas);
% imagesc(zeros(size(IM(:,:,2)))), axis image, colormap(gray)
disp('Computing fish orientation...')
tic

for frame = imgFrames
%     imagesc(IM(:,:,imgNum)), colormap(gray), axis image
    lineSums = 0;
    lineRanges = 0;
    lineProf = {};
    
    for jj = thetaInds
        th{jj} = repmat(lineThetas(jj),1,lineLength);
        [blahX,blahY] = pol2cart(th{jj},0:lineLength-1);
        blahX = blahX + x(frame);
        yy = y;
%         yy = size(IM,1)-y;
        blahY = blahY  + yy(frame);
%         hold on
%          plot(blahX, blahY,'color',rand(3,1))
      
        blahX(blahX>size(IM,2)) = size(IM,2);
        blahX(blahX<1) = 1;
        blahY(blahY>size(IM,1)) = size(IM,1);
        blahY(blahY<1) = 1;        
        img = -IM(:,:,frame);  
        lineInds = sub2ind(size(img),round(blahY),round(blahX));
        
%         hold on        
% %         img_line = img;
% %         img_line(lineInds) = 0;
% %         imagesc(img_line), colormap(gray),axis image, title(num2str(jj))
% %         pause()

        lineProf{jj} = img(lineInds);
        lineSums = [lineSums; sum(img(lineInds))];
        lineRanges = [lineRanges; max(lineProf{jj})-min(lineProf{jj})];        
        cartCoords{jj} = [blahX(:), blahY(:)];
        
    if 0   
        blankImg(lineInds) = lineProf{jj};
        imagesc(blankImg),axis image, colormap(gray)
    end
    
    end
    lineSums(1) = [];
    lineRanges(1) = [];
    lineSums = (lineSums-min(lineSums))/(max(lineSums)-min(lineSums));
    lineRanges = (lineRanges-min(lineRanges))/(max(lineRanges)-min(lineRanges));
    trunkLine = lineSums-lineRanges;
    trunkInd = find(trunkLine==max(trunkLine));
    orientation(frame) = lineThetas(trunkInd(1))*180/pi;     
    orientation(frame) = mod(orientation(frame) + 180,360);

    if mod(frame,500)== 0
        disp(['Img # ' num2str(frame)])
    end
end
toc


