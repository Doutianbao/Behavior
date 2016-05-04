function fishPos = MidlineToFishPos(midlineInds)
% MidlineToFishpos - A useful script to get fish position matrix
%   from midline inds, in case the latter is corrupted or not saved by
%   mistake
% fishPos = MidlineToFishPos(midlineInds);
fishPos = zeros(length(midlineInds),2); 
imgDims =[600,600]; 
for jj = 1:length(midlineInds)
    [r,c]= ind2sub(imgDims,midlineInds{jj}{1});
    fishPos(jj,:) = [c(1),r(1)]; 
end
end
