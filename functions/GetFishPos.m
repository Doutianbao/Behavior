function fishPos = GetFishPos(IM,varargin)
%GetFishPos Give an image stack returns the x,y coordinates of the fish in
%   each of the images of the stack
% fishPos = GetFishPos(IM);
% fishpos = GetFishPos(IM,nPixels,method)
% Inputs:
% IM - Image stack
% nPixels - # of bright pixels for centroid determination
% method - Centroid detection method: {'median'}, ['mean']
% 
% Avinash Pujala, HHMI, 2016

nPixels = 30;
method = 'mean';
if nargin == 2
    nPixels = varargin{1};
elseif nargin == 3
    nPixels = varargin{1};
    method = varargin{2};
end


x = zeros(1,size(IM,3));
y = x;
dispChunk = round(size(IM,3)/50)+1;
for jj=1:size(IM,3)
    img= IM(:,:,jj);  
    [r,c] = FishPosInImg(img,nPixels,method);
    x(jj) = c;
    y(jj) = r;
    if mod(jj,dispChunk)==0
        disp(['Img # ' num2str(jj)])
        ShowFishPos(img,[r,c],jj)
    end
end
fishPos = [x; y]';
   

end

 function [r,c] = FishPosInImg(img,nPixels,method)
        [~,maxInds] = sort(img(:),'descend');
        maxInds = maxInds(1:nPixels);
        [r,c] = ind2sub(size(img),maxInds);
        if strcmpi(method,'median')
            r = round(median(r));
            c = round(median(c));
        elseif strcmpi(method,'mean')
%             r = round(mean(r));
%             c = round(mean(c));
            r = round(sum(r(:).*img(maxInds))/sum(img(maxInds)));
            c = round(sum(c(:).*img(maxInds))/sum(img(maxInds)));
        end
        
    end
    function ShowFishPos(img,fishPos,imgNum)
        cla
        imagesc(img), axis image, colormap(gray)
        hold on
        plot(fishPos(2), fishPos(1),'ro')
        title(['Frame # ' num2str(imgNum)])
        drawnow        
        shg
    end