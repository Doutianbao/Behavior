function out = TransposeCell(cellData)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

maxSize = 0;
for jj = 1:length(cellData)
    maxSize = max(length(cellData{jj}),maxSize);
end

out = cell(maxSize,1);
for ii = 1: maxSize
    blah = [];
    for jj = 1:length(cellData)
        if iscell(cellData(jj))
            var = cellData{jj};
        else
            var = cellData(jj);
        end
        if ii <= length(var)
             blah = [blah; var(ii)];
        end       
    end
    out{ii} = blah;
end

end

