%Present OMR stimulus for any number of intervals of varying duration,
%speed, and angle.

close all
clear all
clc

%User Set
ivlSpeeds=repmat(-0.2,1,8); %negative is OMR inducing, units are percent of total screen (.05 is reasonable rate on Minoru's behavior setup)
durs=5;
rotAngles=[0,90,0,-90,0,90,0,-90]; %degrees from y-axis, min -90, max 90
nBars=10;
screenDim=100; %each side will have this many "pixels"
testMode=0; %in test mode, don't move stimulus to other monitor

%For constant variables, create necessary expanded vector of values
nIvls=max([length(ivlSpeeds),length(durs),length(rotAngles)]);
if length(ivlSpeeds)<nIvls
    ivlSpeeds=ivlSpeeds*ones(nIvls,1);
end
if length(durs)<nIvls
    durs=durs*ones(nIvls,1);
end
if length(rotAngles)<nIvls
    rotAngles=rotAngles*ones(nIvls,1);
end

%Initialize time,duration, and figure
startT=tic;
currentT=toc(startT);
duration = 0;
screenH=figure();

%Start drift intervals
for i= 1:length(ivlSpeeds)
    
    %Decompose rotAngle
    rotAngle=rotAngles(i);
    rotAngleMag=abs(rotAngle);
    
    %adjust rotAngle if = 0 or 90 to prevent undefined error
    if rotAngleMag==0
        rotAngleMag=.000000000001;
    elseif rotAngleMag==90
        rotAngleMag=89.999999;
    end
    rotAngleDir=rotAngle/rotAngleMag;
    
    %Prepare presentation screen (note stimulus will go beyond thesebounds)
    set(gcf,'Color','k')
    if rotAngleDir>=0
        presentationAxisDims=[0,screenDim,0,screenDim];
    elseif rotAngleDir<0
        presentationAxisDims=[-screenDim,0,0,screenDim];
    end
    axis equal;
    axis(presentationAxisDims)
    
    %Move Screen to secondary Display
    if testMode~=1
        figure(screenH)
        set(gcf,'Units','normalized','menubar', 'none','Color','k','NumberTitle','off')
        set(gcf,'position',[.9, 0, 1, 1 ]);
    end
    
    %calculate stim dimensions
    rangeY=screenDim+screenDim*tand(rotAngleMag);
    rangeX=screenDim+screenDim/tand(rotAngleMag);
    lengthStimPerpLine=rangeY*cosd(rotAngleMag);
    xPeriod=rangeX/nBars;
    yPeriod=rangeY/nBars;
    perpLinePeriod= lengthStimPerpLine/(nBars);
    
    %prepare line thickness
    currentaxisunits = get(gca,'units');
    set(gca,'units','points');
    p=get(gca,'position');
    yPtsPerUnit=p(4)/screenDim; %for some reason p(3) gives an overestimate of what Id expect for screen dim in X, so only p(4) (screen dim in Y) works for this purpose
    % xPtsPerUnit=p(4)/screenDim;
    % xBarPts=(xPeriod*xPtsPerUnit);
    yPerPts=(yPeriod*yPtsPerUnit);
    perpLinePerPts=perpLinePeriod*yPtsPerUnit;
    lineWid=perpLinePerPts/2;
    set(gca,'units',currentaxisunits);
    clear currentaxisunits;
    
    %Prepare Speed and duration
    speed=ivlSpeeds(i);
    speed=speed*screenDim;%Convert speed to screen units
    duration=duration+durs(i);
    %initialize counter
    counter=0;
    %Move Bars
    while currentT < duration
        stepSize=toc(startT)*speed;
        %decompose stepsize into x and y coordinates
        stepSize=rem(stepSize,perpLinePeriod); %start over after each period
        pctShift=stepSize/perpLinePeriod;
        xShift=pctShift*xPeriod; %cosd(rotAngle)*stepSize;
        yShift=pctShift*yPeriod;%sind(rotAngle)*stepSize;
        %get x and y intercepts of bars
        xCoords=(-xPeriod:xPeriod:rangeX) + xShift;
        yCoords=(-yPeriod:yPeriod:rangeY) + yShift;
        %add in zero coordinates
        xCoords=[xCoords;zeros(size(xCoords))];
        yCoords=[zeros(size(yCoords));yCoords];
        %if rotation is CCW, flip across x-axis
        if rotAngleDir<0
            xCoords=-xCoords;
        end
        %display image
        figure(screenH);
        plot(xCoords,yCoords,'k', 'LineWidth',lineWid); %plot black bars
        axis equal;
        axis (presentationAxisDims);
        set(gca,'Color','w','Visible','on','xtick',[],'ytick',[]) %axis is white, so if turn it back on make sure to plot bars in black
        if testMode==1 %plot fish orientation (sitm is flipped relative to camera)
            hold on
            plot([50,50],[50,40],'y','LineWidth',8);
            plot([45,50],[45,40],'y','LineWidth',5);
            plot([55,50],[45,40],'y','LineWidth',5);
            hold off
        end
        pause(.005);
        currentT=toc(startT);
        if currentT-sum(durs(1:i-1))>counter
            display(strcat('speed:',num2str(speed),' time:',num2str(counter)));
            counter=counter+1;
        end
    end
end
