function pks = GetSwimPksByCycle(X,maxPkDist,pkThr)
%GetSwimPksByCycle Gets swim cyle pks separately for the first 3 cycles
% pks = GetSwimPksByCycle(X,maxPkDist,pkThr);
%Inputs:
% X - m x n timeseries variable(s) where m = # of observations, n = # of
%   variables
% maxPkDist - Maximum distance between peaks before they are treated as
%   belonging to different swim episodes
% pkThr - Threshold for peak detection as defined in "peakdet.m" 
% Outputs:
% pks - Indices of peaks. pks is a cell array of size p x n x 3 x 2, 
%   where p = # of peaks, n = # of timeseries (variables) in input X, 
%   3 corresponds to the first 3 peaks, and 2 beccause the maxima and
%   minima are returned separately
% 
% Avinash Pujala, Koyamalab/HHMI, 2016

if nargin < 3
    error('Minimum 3 inputs required!')
end

sigPks = cell(size(X,2),1);
for sig = 1:size(X,2)
    maxMinPks = cell(2,1);
    [maxtab,mintab]  = peakdet(X(:,sig),pkThr);
    pks = union(maxtab(:,1),mintab(:,1));
    pks_sort = SortPksByCyle(pks,maxPkDist);
%     pks_sort = SortPksByCyle(maxtab(:,1),maxPkDist);
%     maxMinPks{1} = pks_sort;
%     pks_sort = SortPksByCyle(mintab(:,1),maxPkDist);
%     maxMinPks{2} = pks_sort;
%     sigPks{sig} = maxMinPks;
    sigPks{sig} = pks_sort;
end

pks = sigPks;

end

function pks_sort =SortPksByCyle(pks,maxPkDist)
pks_sort = cell(3,1);
dPks = diff(pks);
pks_sort{1} = union(pks(1),pks(find(dPks>maxPkDist)+1));
for cyc = 1:2
    temp = [];
    count = 1;
    for pkNum = 1:numel(pks_sort{cyc})
        pk = pks_sort{cyc}(pkNum);
        blah = pks;
        blah(blah<=pk)=[];
        dPks = blah-pk;
        dPks(dPks > maxPkDist)=[];
        if ~isempty(dPks);
        [~,closeInd] = min(dPks);
        temp(count) = blah(closeInd);        
        count = count + 1;
        end
    end   
    pks_sort{cyc+1}= temp;
end
end
