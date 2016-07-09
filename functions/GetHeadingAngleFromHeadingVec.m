function hAngle = GetHeadingAngleFromHeadingVec(hOrVecInds,imgDims)
%GetHeadingAngleFromHeadingVec - Given the indices of the head direction
%   vector and the image dimensions, returns the heading angle
% hAngle = GetHeadingAngleFromHeadingVec(hOrVecInds,imgDims);
% Inputs:
% hOrVecInds - Indices of head direction vectos as obtained by GetFishPos.
%   hOrVecInds is a cell array, where each cell holds the indices of the
%   heading vector obtained from a single image
% imgDims - Image dimensions (2D).
% Outputs:
% hAngle - Heading angle in image 2D space (0 to 360 degrees)
%
% Avinash Pujala, Koyama lab/HHMI, 2016

imgDims = imgDims(1:2); % Only works in 2D at present
N = length(hOrVecInds);
hAngle = zeros(N,1);
for n  = 1:N
    if size(hOrVecInds{n},2)==1
         [y,x] = ind2sub(imgDims,hOrVecInds{n});
    else
        x = hOrVecInds{n}(:,1);
        y = hOrVecInds{n}(:,2);
    end
   
    hAngle(n) =  angle(x(1)-x(end) + (y(1)-y(end))*1i) * 180/pi;  
end

hAngle(hAngle<0) = 360-(hAngle(hAngle<0) + 360);
dH = diff(hAngle);
dH(dH> 180) = 360-dH(dH>180);
hAngle = cumsum([hAngle(1);dH(:)]);

end

