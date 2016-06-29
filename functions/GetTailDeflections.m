function varargout = GetTailDeflections(midlineInds,tailCurv,imgDims,varargin)
%GetTailDeflections Given midlineInds, and tailCurv generated by
%   GetMidlines and SmoothenMidlines respectively, returns a series of
%   tail deflection angles estimated by tangent lines drawn at different
%   points on the tail. This is inspired by analysis described in Huang et
%   al., 2013
% tailAngles = GetTailDeflections(midlineInds,tailCurv,nAngles);
% [tailAngles,dS] = GetTailDeflections(...);
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
%
% Avinash Pujala, Koyama lab/HHMI, 2016

plotBool = false;

nAngles = size(tailCurv,1)-1;
if nargin > 3
    nAngles = min(varargin{1},nAngles);
end
imgDims = imgDims(1:2);

hAngles = zeros(length(midlineInds),1);
tAngles = zeros(nAngles,length(midlineInds));
tInds = round(linspace(1,size(tailCurv,1),nAngles+1));
tInds(tInds==0) = 1;
tInds(tInds>size(tailCurv,1))=size(tailCurv,1);
for tt = 1:length(midlineInds)
    [y_h, x_h] = ind2sub(imgDims,midlineInds{tt}{1});
    c1 = (x_h(end)-x_h(1)) + (y_h(end)-y_h(1))*1i;
    hAngles(tt) = angle(c1);
    if plotBool
        cla
        hold on
        axis image
        axis([1 imgDims(1) 1 imgDims(2)])
        plot(tailCurv(:,1,tt),tailCurv(:,2,tt),'o-')
        set(gca,'color','k')
    end
    c2 = []; 
    for ang = 1:nAngles
        x_t = tailCurv(tInds(ang):tInds(ang+1),1,tt);
        y_t = tailCurv(tInds(ang):tInds(ang+1),2,tt);
        c2_old = c2;
        c2 =  (x_t(end)-x_t(1) + (y_t(end)-y_t(1))*1i);
        if ang >1
            tAngles(ang,tt) = angle(c2_old*conj(c2));
%               tAngles(ang,tt) = angle(c2);
        
        else
            tAngles(ang,tt) = angle(c1*conj(c2));
%             tAngles(ang,tt) = angle(c2);
        end
        if plotBool
            plot(x_t,y_fit,'g-')
            plot(x_t,y_t,'r-')
        end
    end    
    if plotBool
        title(num2str(tt))
        shg
%         pause()
    end
end
varargout{1} = tAngles;
varargout{2} = hAngles;

end

