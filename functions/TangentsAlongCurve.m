function varargout = TangentsAlongCurve(curv)
%TangentsAlongCurv - Given a curv (x,y coords), returns tangents along each
%   point on the curve w.r.t the first point
% [th1,th2,th3] = TangentsAlongCurve(curve);
% Inputs:
% curve - 2D curve
% Outputs:
% th1 - Angles of tangent lines drawn at various points along the curve (+ve and -ve)
% 
% Avinash Pujala, Koyama lab/HHMI, 2016

if size(curv,2)~=2
    error('Curve must be 2D, check input!')
end

dC  = diff(curv,[],1);
z = dC(:,1) + dC(:,2)*1i;

thetas = zeros(size(dC,1),1);
thetas_abs = thetas;
for pt = 1:size(dC,1)
    thetas(pt) = angle(conj(z(1))*z(pt))*180/pi; 
    thetas_abs(pt) = angle(z(pt))*180/pi;
end
varargout{1} = thetas;
end
