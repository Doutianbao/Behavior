
%% Turn stuff
x = tracexy_flt(1,:);
yy = y;
y = size(IM,1)-tracexy_flt(2,:);
dx = diff(x);
dy = diff(y);
S = sqrt(dx.^2 + dy.^2);
S_sort = sort(S,'ascend');
S_perc = S_sort(round(length(S)/2));
noMovInds = find(S<=mean(S));
dx(noMovInds) = []; 
dy(noMovInds) = [];
S(noMovInds) = [];

c = dx + dy*i;
a = zeros(size(dx));
for jj = 1:length(dx)-1
    a(jj+1) = angle(c(jj)*conj(c(jj+1)))*(180/pi);
end
b = a;
b(abs(b)<=10) = [];

leftVec = zeros(size(dx));
leftVec(a<0) = 1;
leftInds = find(leftVec)+1;

rightVec = zeros(size(dx));
rightVec(a>0) = 1;
rightInds = find(rightVec)+1;

straightInds = setdiff(1:length(x),union(leftInds,rightInds));

% Test turn estimation
% figure
% imagesc(ref),axis image, colormap(gray)
% hold on
% startFrame = 1;
% endFrame = 500;
% 
% for f = startFrame:endFrame
%     plot(tracexy_flt(1,f),tracexy_flt(2,f),'ro')
%     %     title(['Frame: ' num2str(f) ', Orientation: ' num2str(orVec(f)) ])
%     if isempty(intersect(f, noMovInds))
%         if intersect(f,leftInds)
%             title(['Frame: ' num2str(f) ', Left' ])
%         elseif intersect(f,rightInds)
%             title(['Frame: ' num2str(f) ', Right' ])
%         else
%             title(['Frame: ' num2str(f) ', Straight' ])
%         end
%         %     title(['Frame: ' num2str(f) ', Turn: ' num2str(orVec(f)) ])
%         shg
%         pause()
%     end
% end
% % break;

 




clear i angle
% thetas = angle(dx + dy*i)*(180/pi);
thetas = atan(dy./dx)*(180/pi);

dTh = diff(thetas);

leftInds = find(dTh > 0);
rightInds = find(dTh < 0 );
straightInds = setdiff(1:length(x),union(leftInds,rightInds));



orVec = zeros(size(S));
or1 = find(dx >0 & dy > 0);
or2 = find(dx <0 & dy > 0);
or3 = find(dx <0 & dy < 0);
or4 = find(dx > 0 & dy <0);
orVec(or1+1) = 1;
orVec(or2+1) = 2;
orVec(or3+1) = 3;
orVec(or4+1) = 4;
for jj = 2:length(orVec)
    if orVec(jj)==0
        orVec(jj) = orVec(jj-1);
    end
end



% theta2 = atan(dy(2:end)./dx(2:end))*(180/pi);

% theta1 = atan(y(1:end-1)./x(1:end-1)) *(180/pi);
% theta2 = atan(y(2:end)./x(2:end))*(180/pi)

orVec = orVec(1:end-1);
rightVec = zeros(size(thetas));
leftVec = zeros(size(thetas));
% dTh = theta2-theta1;
rightVec(((orVec(2:end) ==1) | (orVec(2:end)==4)) & (thetas(2:end)  < 0))=1;
rightVec(((orVec(2:end) ==3) | (orVec(2:end)==4)) & (thetas(2:end)  >0))=1;

leftVec(((orVec(2:end) ==1) | (orVec(2:end)==2)) & (dTh(2:end)  > 0))=1;
leftVec(((orVec(2:end) ==3) | (orVec(2:end)==4)) & (dTh(2:end)  <0))=1;

% 
% rightVec(theta1>5) = 1;
% leftVec(theta1 <-5) = 1;

rightInds = find(rightVec);
leftInds = find(leftVec);


% Test turn estimation
figure
imagesc(ref),axis image, colormap(gray)
hold on
startFrame = 500;
endFrame = 2000;

for f = startFrame:endFrame
    plot(tracexy_flt(1,f),tracexy_flt(2,f),'ro')
%     title(['Frame: ' num2str(f) ', Orientation: ' num2str(orVec(f)) ])
    if f~=1    
    if intersect(f,leftInds)
      title(['Frame: ' num2str(f) ', Left' ])
    elseif intersect(f,rightInds)
       title(['Frame: ' num2str(f) ', Right' ])
    else
        title(['Frame: ' num2str(f) ', Straight' ])
    end
%     title(['Frame: ' num2str(f) ', Turn: ' num2str(orVec(f)) ])
    shg  
     pause()
    end
end 
break

S_right = S(rightInds+1);
S_left = S(leftInds + 1);

[rCount,rVals] = hist(S_right,50);
[lCount,lVals] = hist(S_left,50);

dRightInds = diff(rightInds);
dLeftInds = diff(leftInds);

contRightInds = find(dRightInds==1);
contLeftInds = find(dLeftInds==1);

S_right_cont = S(contRightInds);
S_left_cont = S(contLeftInds);
[rCount,rVals] = hist(S_right_cont,50);
[lCount,lVals] = hist(S_left_cont,50);

figure
plot(lVals,lCount,'linewidth',2)
hold on
plot(rVals,rCount,'r','lineWidth',2)
legend('Left','Right')
xlim([-inf 50])
set(gca,'tickdir','out')
xlabel('Swim distance before turn to other direction')
ylabel('Count')
box off
title('Histogram of swim distances in one direction before turn to other direction')


SS = nan;
blah = 0;
for jj = 1:length(dRightInds)
    if dRightInds(jj)~=0
        blah = blah+S(jj);
    else
        blah = 0;
    end
    SS = [SS; blah];
end
SS_right = SS;

SS = nan;
blah = 0;
for jj = 1:length(dLeftInds)
    if dLeftInds(jj)~=0
        blah = blah+S(jj);
    else
        blah = 0;
    end
    SS = [SS; blah];
end
SS_left = SS;



%% Plot trajectory
figure('Name','Fish trajectory')
imagesc(ref), colormap(gray), axis image, axis off
hold on
plot(tracexy_flt(1,1:end-3),tracexy_flt(2,1:end-3),'r')
title('Fish trajectories over a 20 min period')

break;



%% Animate trajectory
startFrame = 800;
endFrame = 3000;
figure('Name','Trajectory animation')
imagesc(IM(:,:,startFrame)), axis image, colormap(gray), axis off
hold on
count = 0;
for f = startFrame:endFrame    
    if mod(count,50)==0
        cla
        imagesc(IM(:,:,f)), axis image, colormap(gray), axis off, drawnow
    else
        cla
        imagesc(IM(:,:,f)), axis image, colormap(gray), axis off, drawnow
    end
    if isempty(intersect(f+1,noMovInds))
        if  count == 0
            plot(tracexy_flt(1,f), tracexy_flt(2,f),'ro'), drawnow
        else
            plot(tracexy_flt(1,f), tracexy_flt(2,f),'ro-'), drawnow
        end
        count = count + 1;
%         if intersect(f,rightInds)
%             title(['Frame # ' num2str(f) ', Right'])
%         elseif intersect(f,leftInds)
%             title(['Frame # ' num2str(f) ', Left'])
%         else
%             title(['Frame # ' num2str(f) ', Straight'])
%         end
        title(['Orientation quad ' num2str(orVec(f-1))])
%         title(['dx = ' num2str(dx(f-1)), ', dy  = ' num2str(dy(f-1))])
        shg
        pause()
    end
    
end


