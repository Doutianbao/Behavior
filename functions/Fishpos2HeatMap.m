function varargout = Fishpos2HeatMap(fishPos,imgDims)
%FishPos2HeatMap Returns a heatmap (matrix) which indicates how long the
%   fish has spent in different parts of the arena
% trajHeatMap = Fishpos2HeatMap(fishPos,imgDims)

heatMap = zeros(imgDims(1:2));
fishPos = floor(fishPos); 
fishPos(fishPos==0)=1;
for fp = 1:size(fishPos,1)
    heatMap(fishPos(fp,1), fishPos(fp,2)) = heatMap(fishPos(fp,1), fishPos(fp,2)) + 1;
end

varargout{1} = heatMap;

end

