imgDims =[600, 600];
f = 10;
nPhases = 10;
phi = linspace(0,pi,nPhases);
lenDiag = ceil(sqrt(sum(imgDims.^2)))
x = linspace(0,1,lenDiag);
figure
for jj = 1:nPhases
y = square(2*pi*f*x - phi(jj));
% y = sin(2*pi*f*x - phi(jj));
[X,Y] = meshgrid(x,y);    
cla
Y= imrotate(Y,45,'crop');
imagesc(Y),colormap(gray)
axis([0 imgDims(1) 0 imgDims(2)])
drawnow
pause(0.1)
end