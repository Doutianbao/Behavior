function orientation = GetFishOrientationFromMidlineInds(midlineInds,imgDims)
%GetFishOrientationFromMidlineInds - Given midlineInds generated by
%   GetMidline.m, and image dimensions, returns a vector of fish
%   orientations, where each element of the vector corresponds to the
%   fish's orientation in a single image frame
% orientation = GetFishOrientationFromMidlineInds(midlineInds,imageDims)
% orientation  = GetFishOrientationFromMidlineInds(midlineInds,imageDims)
% Inputs:
% midlineInds - This is the output variable returned by the function
%   GetMidline. This variable is a cell array of length T, where T is
%   number of image frames (time points). Each cell is itself a cell array
%   of length L, where L is the number of line segments superimposed on
%   the fish. The number of points in each cell, l, of the L cell arrays
%   corresponds to the length of that particular line segment superposed on
%   the fish. The first cell array l1 corresponds to the line used to
%   determine the fish's orientation
% imgDims - Dimensions of the each from which midlineInds were extracted.
%   Img Height x Img Width
% Outputs:
% orientation - Length T vector where each element is the orientation in
%   single image frame of the image stack of size imgDims

if numel(imgDims)>2
    error('Image dimensions must be a 2 element vector')
end

orientation = nan(length(midlineInds{1}),length(midlineInds));
for seg = 1:size(orientation,2)
    for imgNum = 1:size(orientation,1)
        mlInds = midlineInds{seg}{imgNum};
        [x, y] = ImageIndsToXYCoords(mlInds,imgDims);
        x = x-x(1);
        y = y-y(1);
        orientation(imgNum,seg) = cart2Angle(x(end),y(end));
    end
end
orientation(orientation<0)= orientation(orientation<0) + 360;
orientation = mod(orientation+180,360);
orientation = orientation';
end

function [x, y] =  ImageIndsToXYCoords(imgInds,imgDims)
[I,J] = ind2sub(imgDims,imgInds);
x = J;
y = I;
% y = imgDims(1)-I;
end

function angle = cart2Angle(x,y)
[theta,~] = cart2pol(x,y);
angle  = theta*180/pi;
end