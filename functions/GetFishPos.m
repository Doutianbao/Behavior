function fishPos = GetFishPos(IM,varargin)
%GetFishPos Give an image stack returns the x,y coordinates of the fish in
%   each of the images of the stack
% fishPos = GetFishPos(IM);
% fishpos = GetFishPos(IM,nPixels,method)
% Inputs:
% IM - Image stack
% nPixels - # of bright pixels for centroid determination
% method - Centroid detection method: {'median'}, ['mean']

nPixels = 30;
method = 'median';
if nargin == 2
    nPixels = varargin{1};
elseif nargin == 3
    nPixels = varargin{1};
    method = varargin{2};
end


x = zeros(1,size(IM,3));
y = x;

for jj=1:size(IM,3)
    ii= IM(:,:,jj);  
    [~,maxInds] = sort(ii(:),'descend');
    maxInds = maxInds(1:nPixels);
    [r,c] = ind2sub(size(ii),maxInds);
    if strcmpi(method,'median')
        r = round(median(r));
        c = round(median(c));
    elseif strcmpi(method,'mean')
        r = round(mean(r));
        c = round(mean(c));
    end
    
    x(jj) = r;
    y(jj) = c;
    if mod(jj,100)==0
        disp(['Img # ' num2str(jj)])
    end
end
fishPos = [y; x]';

end

