function fishPos = GetFishPos(IM)
%GetFishPos Give an image stack returns the x,y coordinates of the fish in
%   each of the images of the stack
% fishPos = GetFishPos(IM);

x = zeros(1,size(IM,3));
y = x;

for jj=1:size(IM,3)
    ii= IM(:,:,jj);  
    
    [~,maxInds] = sort(ii(:),'descend');
    maxInds = maxInds(1:30);
    [r,c] = ind2sub(size(ii),maxInds);
    r = round(mean(r));
    c = round(mean(c));
    x(jj) = r;
    y(jj) = c;
    if mod(jj,100)==0
        disp(['Img # ' num2str(jj)])
    end
end
fishPos = [y; x]';

end

