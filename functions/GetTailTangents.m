function varargout = GetTailTangends(tailCurv,varargin)
%GetTailDeflections Given midlineInds, and tailCurv generated by
%   GetMidlines and SmoothenMidlines respectively, returns a series of
%   tail deflection angles estimated by tangent lines drawn at different
%   points on the tail. This is inspired by analysis described in Huang et
%   al., 2013
% tailAngles = GetTailDeflections(midlineInds,tailCurv,nAngles);
% [tailAngles,headAngles,tailAngles2] = GetTailDeflections(...);
%
% Inputs:
% midlineInds - Indices of a series of line segments of different lengths
%   placed along the fish tail (output of GetMidlines)
% tailCurv - Smoothened version of midlineInds as returned by
%   SmoothenMidlines
% imgDims - Image dimensions. This is used to convert indices to
% subscripts.
% nAngles - The # of tail angles to return for each time point. By default,
%   nAngles = size(tailCurv,1)-1.
% Outputs:
% tailAngles - Matrix of size nAngles-by-T, where T = size(tailCurv,1) =
%   number of time points. Tail angles are w.r.t the heading vector, which
%   for a given time point t, is given by midlineInds{t}{1};
% headAngles - Heading direction
% dS - distances between successive pts on the tail
%
% Avinash Pujala, Koyama lab/HHMI, 2016


nAngles = size(tailCurv,1)-2;
if nargin ==2
    nAngles = min(varargin{1},nAngles);
end
% imgDims = imgDims(1:2);

segLen = max(1,round(size(tailCurv,1)/(nAngles+2)));
tc = tailCurv(1:segLen:end,:,:);
tAngles = zeros(size(tc,1)-2,size(tailCurv,3));

for tt = 1:size(tAngles,2);
    dTC = diff(squeeze([tc(:,1,tt),tc(:,2,tt)]));
    C = dTC(:,1) + dTC(:,2)*1i;
    A = angle(C(1:end-1).*conj(C(2:end)));
    tAngles(:,tt) = cumsum(A);    
end

tA  = tAngles;
ker = gausswin(max(1,round(size(tA,1)/6)));
ker = ker/sum(ker);
disp('Smoothing tail angles...')
for tt = 1:size(tA,2)
    tA(:,tt) = conv2(tA(:,tt),ker(:),'same');
end

tA_flt = zeros(size(tA));
tA2 = tA;
tA2(isnan(tA))=0;
for ang = 1:size(tA,1)
    tA_flt(ang,:) = chebfilt(tA2(ang,:),1/500,50,'low');
end

tAngles = tA_flt;

varargout{1} = tAngles*180/pi;

end

