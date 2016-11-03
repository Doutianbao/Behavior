function varargout = FishSwim_func(varargin)
%FishSwim_func  A function for tracking fish in an image sequence, cropping
%   images around the tracked fish, getting kinematic timeseries, and
%   storing in a matfile in a subdirectory within the image directory.
% procData = FishSwim_func();
% procData = FishSwim_func(imgDir);
% procData = FishSwim_func(imgDir,'imgExt',imgExt,'imgInds',imgInds,'spatialFilt',spatialFish,
%   'fps',fps,'blockSize',blockSize,'cropWid',cropWid,'saveOrNot',saveOrNot);
% Inputs:
% imgDir - Image directory. If empty, interactively allows choosing of
%   image directory by way of selection of a file within that directory
% imgExt - Image exension, in case a directory has mixed image type.
% imgInds - Indices of the sequence of images to read.
% spatialFilt - A 2, or 1 element variable that specifies how images should
%    be filtered before detecting fish in them. If numel(spatialFilt) ==2,
%    then bandpasses, else smooths with kernel of width specified by
%    spatialFilt.
% 



end

