function orientation = GetFishOrientation(IM,fishPos,varargin)
%GetFishOrientation - Given an image stack gets the orientaton of the fish
%   in each image of the stack using the fish position input
% orientation  = GetFishOrientatin(IM,fishPos,lineLength)
dTh = 10;
if nargin < 3
    lineLength = 25;
elseif nargin == 3
    lineLength = varargin{1};
end
x = fishPos(:,1);
y = fishPos(:,2);
thetas = (0:dTh:360)*pi/180;
th = {};
cartCoords = {};
imgNums = 2:size(IM,3);
orientation = zeros(size(x));
orLine = {};
orLine{1} = 0;
% figure
disp('Computing fish orientation...')
tic
for imgNum = imgNums
%     imagesc(IM(:,:,imgNum)), colormap(gray), axis image
    lineSums = 0;
    lineRanges = 0;
    lineProf = {};
    for jj = 1:numel(thetas)
        th{jj} = repmat(thetas(jj),1,lineLength);
        [blahX,blahY] = pol2cart(th{jj},0:lineLength-1);
        blahX = blahX + x(imgNum);
        yy = y;
%         yy = size(IM,1)-y;
        blahY = blahY  + yy(imgNum);
%         hold on
%          plot(blahX, blahY,'color',rand(3,1))
           
        blahX(blahX>size(IM,2)) = size(IM,2);
        blahX(blahX<1) = 1;
        blahY(blahY>size(IM,1)) = size(IM,1);
        blahY(blahY<1) = 1;        
        img = -IM(:,:,imgNum);  
        lineInds = sub2ind(size(img),round(blahY),round(blahX));
%         hold on
%         img_line = img;
%         img_line(lineInds) = 0;
%         imagesc(img_line), colormap(gray),axis image, title(num2str(jj))
%         pause()
        lineProf{jj} = img(lineInds);
        lineSums = [lineSums; sum(img(lineInds))];
        lineRanges = [lineRanges; max(lineProf{jj})-min(lineProf{jj})];        
        cartCoords{jj} = [blahX(:), blahY(:)];
    end
    lineSums(1) = [];
    lineRanges(1) = [];
    lineSums = zscore(lineSums);
    lineRanges = zscore(lineRanges); 
    inds1 = find(lineRanges < -0.5);
    orSum = max(lineSums(inds1));
    orLine{imgNum} = find(lineSums == orSum);
    if (numel(orLine{imgNum})>1) && (imgNum > 1)
        [~,ind] = min(abs(orLine{imgNum}-orLine{imgNum-1}));
        orLine{imgNum} = orLine{imgNum}(ind);
    elseif (numel(orLine{imgNum})>1) && (imgNum ==1)
        orLine{imgNum} = orLine{imgNum}(1);
    end
     orientation(imgNum) = mod(180+(thetas(orLine{imgNum})*180/pi),360);   
  
%     hold on
%     plot(cartCoords{orLine{imgNum}}(:,1),cartCoords{orLine{imgNum}}(:,2),'color','r','linewidth',2)
%     title(['Img # ' num2str(imgNum),' Orientation: ' num2str(orientation(imgNum))])
%     shg
%     pause()
    if mod(imgNum,1000)== 0
        disp(['Img # ' num2str(imgNum)])
    end
end
toc
