

% switch  'LoadingNewFilm'  %'LoadingNewFilm' 'RerunAnalysis' 'LoadingCoordinates'
%     case 'LoadingNewFilm'
%

clear all
close all

cd 'Z:\Avinash\Ablations & Behavior'

readMode =  'fromImages'; %'fromMishVid';

poolSize  = 10;
switch readMode
    case 'fromMishVid'
        [IM, outDir] = ReadMishVid();        
    case 'fromImages'
        imgDir = input('Enter image dir path:  ', 's')
        imgExt = input('Enter image extension, e.g. jpg:  ','s')
        imgInds = input('Enter indices of images to read as a vector:  ');
        fName_prefix = input('Enter fish name, e.g., Fish1: ','s');
%         IM = ReadImgSequence_beta(imgDir,imgExt,imgInds);
         IM = ReadImgSequence(imgDir,imgExt,imgInds);
        outDir = fullfile(imgDir,'spont');
end

if exist(outDir)~=7
    mkdir(outDir)
end

%% Processing images
tic
if matlabpool('size')==0
    matlabpool(poolSize)
end
disp('Processing images...')
IM_proc = ProcessImages(IM);

%% Tracking the fish
disp('Tracking fish...')
fishPos = GetFishPos_parallel(IM_proc, 40);
toc
disp('Creating mean reference frame...')
ref = mean(IM,3);
toc

%% Fish Orientation
disp('Getting fish orientation...')
tic
IM_orient = max(IM_proc(:))-IM_proc;
% midlineInds = GetMidline_template_parallel(IM_orient,fishPos,[30]);

break
tic
midlineInds_parallel = GetMidline_parallel(IM,fishPos,[30 24 18]);
toc

%    orientation_corr = CorrectOrientation(orientation, 90);
imgDims = size(IM_proc);
orientation = GetFishOrientationFromMidlineInds(midlineInds,imgDims(1:2));
orientation = orientation';
orientation_backup = orientation;
toc

%% Motion Info
% motionThr = 5;
% [motionFrames, swimStartFrames] = GetMotionFrames(fishPos,motionThr);
% motionInfo = GetMotionInfo(fishPos,orientation,imgDims(1));

%% Save timeseries
% fName = input('Enter fish name (e.g. Fish7): ','s');
fName_suffix = [num2str(round(imgInds(1)/60/30)) '-' num2str(round(imgInds(end)/60/30)) 'mins'];
fName = strcat(fName_prefix,'_',fName_suffix);
ts = datestr(now,30);
save(fullfile(outDir,[fName, '_orientation_' ts '.mat']),'orientation');
save(fullfile(outDir,[fName, '_imgDims_'  ts '.mat']),'imgDims')
save(fullfile(outDir,[fName, '_midlineInds_' ts '.mat']),'midlineInds');
save(fullfile(outDir,[fName, '_ref_' ts '.mat']),'ref');
save(fullfile(outDir,[fName, '_tracexy_' ts '.mat']),'fishPos');
% save(fullfile(outDir,[fName '_motionInfo_' ts '.mat']),'motionInfo');
disp(['Saved orientation, imgDims ,midlineInds, ref, tracexy at ' outDir])
%

%% Saving processed images
saveOrNot = 'y';
%         saveOrNot = input('Save the variables (y/n)?  ','s');
tic
if strcmpi('y',saveOrNot)
    disp('Saving relevant variables...') 
    savefast(fullfile(outDir,[fName, '_IM_proc.mat']),'IM_proc');    
else
    disp('Data not saved!')
end
toc
if matlabpool('size')>0
    matlabpool close
end
break;

%% Turn angles during swims and histogram
x = fishPos(:,1);
y = fishPos(:,2);
dx = diff(x);
dy = diff(y);
dS = sqrt(dx.^2 + dy.^2);
noMovInds  = find(dS < 10);
movInds = setdiff(1:length(dx),noMovInds);
x_mov = x(movInds);
y_mov = y(movInds);
dx_mov = zeros(size(movInds));
dy_mov = zeros(size(movInds));

orientation_mov = fish.orientation(movInds);
dO_mov = diff(orientation_mov);
dO_mov(dO_mov <-180) = dO_mov(dO_mov < -180) + 360;
dO_mov(dO_mov >180) = dO_mov(fish.dO_mov > 180) - 360;

figure('Name','Turn angle histogram')
hist(fish.dO_mov,75)
box off
set(gca,'tickdir','out','xtick',-180:45:180)
xlabel('Turn angle')
ylabel('Count')
title('Turn angle histogram')


%% Swim distances and histogram
dMovInds = diff(movInds);
movInds_jump = movInds;
movInds_jump(dMovInds==1)=[];
x_mov = x(movInds_jump);
y_mov = y(movInds_jump);
dS_mov = sqrt(x_mov.^2 + y_mov.^2);

figure('Name','Swim distance histogram')
hist(dS_mov,50)
box off
xlabel('Swim distances (px)')
title('Swim dist histogram')
ylabel('Count')
set(gca,'tickdir','out')
title('Swim distance histogram')
shg


break

%
%
%
%
%     case 'RerunAnalysis'
%     case 'LoadingCoordinates'
%         clear all
%         close all
%
%         cd X:\YuMu\free_swim
%
%         fr_rate = 1500; %3000 frames per min
%         [FileName,PathName] = uigetfile('*.mishVid*','Select the mishVid');
%         [~,fName,~]=fileparts(FileName);
%         outDir=[PathName,fName,'\swims\'];
%         load([outDir,fName, '_tracexy.mat']);
%         load([outDir,fName, '_ref.mat']);
%         tic
%         load([outDir,fName, '_IM.mat']);
%         toc
%         y_ = fishPos(:,1)';
%         x_ = fishPos(:,2)';
% end
%
%
% %%
% figure;imagesc(ref); axis image
% colormap(gray);
% [diax diay]=ginput(2);
% dia = sqrt((diff(diax)).^2 + (diff(diay)).^2 );
%
% calculatingR1 = dia/2;
% calculatingR2 = dia/2 * .95;
% calculatingR3 = dia/2*sqrt(0.85);
% calculatingR4 = dia/2 * .9;
% calculatingR5 = dia/2 * .85;
% x0 = size(ref,2)/2;
% y0 = size(ref,1)/2;
%
% hold on
% rectangle('Position',[x0-calculatingR1,y0-calculatingR1,2*calculatingR1,2*calculatingR1],'Curvature', [1 1]);
% rectangle('Position',[x0-calculatingR2,y0-calculatingR2,2*calculatingR2,2*calculatingR2],'Curvature', [1 1]);
% rectangle('Position',[x0-calculatingR3,y0-calculatingR3,2*calculatingR3,2*calculatingR3],'Curvature', [1 1]);
% rectangle('Position',[x0-calculatingR4,y0-calculatingR4,2*calculatingR4,2*calculatingR4],'Curvature', [1 1]);
% rectangle('Position',[x0-calculatingR5,y0-calculatingR5,2*calculatingR5,2*calculatingR5],'Curvature', [1 1]);
%
%
%
% %% smoothing
% if 1
%     loc = [y_;x_];
%
%     L = [];
%     s_fac = 1;  %2 for small chamber
%     s_range = 2;
%     ker = exp(-((-s_range):s_range).^2/(2*s_fac^2));
%     ker = ker/sum(ker);
%     L(1,:) = conv(loc(1,:),ker);
%     L(2,:) = conv(loc(2,:),ker);
%     L = L(:, (2*s_fac+1):end);
%     %%%%
%     % NOT GOOD
%     % TALK TO MISHA
%     %%%
%
%     y = L(1,:);
%     x = L(2,:);
%
% else
%     s_fac = 0;
%     s_range = 0;
%     x = x_;
%     y = y_;
% end
%
% savename = ['Turn_',num2str(s_fac), '_',num2str(s_range),'_',FileName(1:10)];
% tracexy_flt = [y; x];
%
% savefast(fullfile([outDir,fName, '_tracexy_flt.mat']),'tracexy_flt');
%
%
%
% %% swimming speed
% if 1
%     dd = sqrt(diff(x).^2+diff(y).^2);
%     dd_ = dd*120/640*samplingRate; %in mm
%     [a_, b_] = hist(dd_,100);
%     speed = bar(b_,a_);
%     title(['swimming speed']);
%     %set(gca,'PaperPositionMode','auto');
%     set(gca,'XLim',[0.5 max(dd)]);
%     xlabel('number','fontsize',10);
%     ylabel('speed (mm/sec)','fontsize',10);
%     speed_name = [outDir,savename,'_speed.tif'];
%     saveas(speed,fullfile(speed_name),'tif');
%
% end
%set(gca,'XLim',[0 100]);

% %%
% thld2 = 100;
% thld = 2;
% if 1
%     id_max = find(dd==max(dd(100:end-100)));
%     id_max = id_max(1);
%     figure;
%     plot(dd);
%     fr_num = length(dd);
%     if id_max>500 & id_max<fr_num-500
%         semilogy(dd);
%         set(gca, 'XLim', [ id_max-500 id_max+500 ]);%, 'YLim', [ 0 max(dd)*1.25 ] );
%     elseif id_max < 500
%         semilogy(dd);
%         set(gca, 'XLim', [ id_max-500 id_max+500 ]);%[ id_max id_max+1000 ]);%, 'YLim', [ 0 100 ] );
%     elseif id_max >fr_num-500
%         semilogy(dd);
%         set(gca, 'XLim', [ id_max-500 id_max+500 ]);%[ id_max-1000 id_max ]);%, 'YLim', [ 0 100 ] );
%     else
%     end
%     set(gca,'YLim',[1 max(dd)*1.25]);
%
%
%     [timex thld]=ginput(2);
%     thld2 = thld(2);
%     thld = thld(1);
% end
%
%
%
% %%
%
% peaks = find( (dd>thld)  .*  (dd<thld2));
% if x(peaks(1)) ~= x(1)
%     peaks = [1 peaks];
% end
% %this version connects the previous and last point between swim bouts
% dx = diff(x(peaks));
% dy = diff(y(peaks));
%
%
% angle = atan2(dy,dx);%atan2(dx,dy);
% dangle = diff(angle);
%
%
% % THIS IS A HACK !!!!
% dangle = [0 dangle];
%
% for i = 1:length(dangle)
%     if dangle(i)>pi
%         dangle(i)=-2*pi+dangle(i);
%     elseif dangle(i)<-pi
%         dangle(i)=2*pi+dangle(i);
%     end
% end
%
%
% r_peak = sqrt((x(peaks(1:end-1))-x0).^2+(y(peaks(1:end-1))-y0).^2);
%
% angle_thr_small = 10/180*pi; %big turns might be ambiguous with the turn direction, can be from either side
% angle_thr_big = 10/180*pi; %big turns might be ambiguous with the turn direction, can be from either side
%
%
% %%
%
% bin_x = -pi:pi/15:pi;
% rangenum = 5;
% hist_dangle = zeros(length(bin_x),rangenum);
% hist_dangle_norm = zeros(length(bin_x),rangenum);
%
% turn = struct;
% turn_fig = figure('position',[200 200 2200,1000]);
% turn_ind = struct;
%
% for i = 1:rangenum
%     rangename = ['region' num2str(i)];
%     turn.(rangename).selecTurn = find(r_peak < eval(['calculatingR', num2str(i)]));
%     turn.(rangename).dangle = dangle(turn.(rangename).selecTurn);
%     turn.(rangename).T_id = find((turn.(rangename).dangle > angle_thr_small)|(turn.(rangename).dangle < -angle_thr_big));
%     turn.(rangename).T_R = find(turn.(rangename).dangle > angle_thr_small & turn.(rangename).dangle < pi - angle_thr_big);
%     turn.(rangename).T_L = find(turn.(rangename).dangle < -angle_thr_small & turn.(rangename).dangle > -pi + angle_thr_big); % because left right is flipped
%     turn.(rangename).T_id = [1 turn.(rangename).T_id];
%
%
%     turn.(rangename).hist_dangle=hist(turn.(rangename).dangle,bin_x);
%     turn.(rangename).hist_dangle_norm = turn.(rangename).hist_dangle/sum(turn.(rangename).hist_dangle);
%
%     T_L_num = length(turn.(rangename).T_L );
%     T_R_num = length(turn.(rangename).T_R);
%     turn.(rangename).IT_num = (T_L_num - T_R_num)/(T_L_num + T_R_num);
%
%     T_L_angle = turn.(rangename).dangle(turn.(rangename).T_L);
%     T_R_angle = turn.(rangename).dangle(turn.(rangename).T_R);
%
%     turn.(rangename).IT_angle = (abs(sum(T_L_angle)) - abs(sum(T_R_angle)))/(abs(sum(T_L_angle)) + abs(sum(T_R_angle)));
%
%     turn.(rangename).T_freq = (T_L_num + T_R_num)/(fr_num*10/fr_rate);
%     subplot(2,rangenum+1,i)
%     plot(bin_x,turn.(rangename).hist_dangle)
%     xlim([-2,2])
%     hold on
%     title([FileName(1:14),'  range',num2str(i)],'interpreter','none')
%
%     subplot(2, rangenum+1,rangenum+1)
%     plot(bin_x,turn.(rangename).hist_dangle)
%     xlim([-2,2])
%     hold on
%     title([FileName(1:14),' 5 rangs'],'interpreter','none')
%
%     subplot(2,rangenum+1,i+rangenum+1)
%     plot(bin_x,turn.(rangename).hist_dangle_norm)
%     xlim([-2,2])
%     hold on
%
%     IT_num_str= num2str(turn.(rangename).IT_num);
%     IT_angle_str= num2str(turn.(rangename).IT_angle);
%     IT_freq_str= num2str(turn.(rangename).T_freq);
%     titlestr = sprintf('IT-num =  %s ; \n   IT_angle =  %s ; \n  T-freq = %s ',  IT_num_str,IT_angle_str,IT_freq_str);
%     title( titlestr, 'interpreter','none')
%
%     subplot(2, rangenum+1,2*rangenum+2)
%     plot(bin_x,turn.(rangename).hist_dangle_norm)
%     xlim([-2,2])
%     hold on
%
%     turn_ind.IT_num(i) = turn.(rangename).IT_num;
%     turn_ind.IT_angle(i) = turn.(rangename).IT_angle;
%     turn_ind.T_freq(i) = turn.(rangename).T_freq;
%     turn_ind.hist_dangle_norm(:,i) = turn.(rangename).hist_dangle_norm;
% end
%
% turn_ind.thrld = [thld;thld2];
% savefast(fullfile([outDir,fName,'_turn.mat']),'turn');
% savefast(fullfile([outDir,fName,'_turn_ind.mat']),'turn_ind');
%
% tif_name = [outDir,fName,'_multiregion.fig'];
% saveas(turn_fig,fullfile(tif_name),'fig');
%
%
% %%
%
% if 1
%
%     %for correlation
%     p = polyfit(dangle_nobig(1:end-1),dangle_nobig(2:end),1);
%     yfit = polyval(p,dangle_nobig(1:end-1));
%     yresid = dangle_nobig(2:end) - yfit;
%     SSresid = sum(yresid.^2);
%     SStotal = (length(dangle_nobig(2:end))-1) * var(dangle_nobig(2:end));
%     rsq = 1 - SSresid/SStotal;
%
%     T_L_num = length(T_L);
%     T_R_num = length(T_R);
%     IT_num = (T_L_num - T_R_num)/(T_L_num + T_R_num);
%     IT_angle = (sum(dangle(T_L)) - abs(sum(dangle(T_R))))/(sum(dangle(T_L)) + abs(sum(dangle(T_R))));
%     IT_L = mean( dangle(T_L));
%     IT_R = mean( abs(dangle(T_R)));
%     T_freq = (T_L_num + T_R_num)/(fr_num*10/fr_rate);
%
%     T_num = length(T_id);
%     TB_L = 2 + find(((dangle(1:(T_num - 2))>angle_thr_small).*(dangle(3:T_num)>angle_thr_small).*(dangle(2:(T_num - 1))>angle_thr_small))>0);
%     TB_R = 2 + find(((dangle(1:(T_num - 2))<-angle_thr_small).*(dangle(3:T_num)<-angle_thr_small).*(dangle(2:(T_num - 1))<-angle_thr_small)>0));
%     TB = 2 + find((((dangle(1:(T_num - 2))>angle_thr_small).*(dangle(3:T_num)>angle_thr_small).*(dangle(2:(T_num - 1))>angle_thr_small))>0)|((dangle(1:(T_num - 2))<-angle_thr_small).*(dangle(3:T_num)<-angle_thr_small).*(dangle(2:(T_num - 1))<-angle_thr_small)>0));
%
%     TB_L_ind = 3*ones(size(TB_L));
%     TB_R_ind = 3*ones(size(TB_R));
%
%     for i = 1:length(TB_L)
%         for b = 1:20
%             if TB_L(i)+b < length(dangle) & dangle(TB_L(i))*dangle(TB_L(i)+b)>0
%                 TB_L_ind(i)= TB_L_ind(i) + 1;
%             else
%                 break
%             end
%         end
%     end
%
%     for i = 1:length(TB_R)
%         for b = 1:20
%             if TB_R(i)+b < length(dangle) & dangle(TB_R(i))*dangle(TB_R(i)+b)>0
%                 TB_R_ind(i)= TB_R_ind(i) + 1;
%             else
%                 break
%             end
%         end
%     end
%
%     TB_L_num = length(TB_L);
%     TB_R_num = length(TB_R);
%     TB_num = length(TB);
%     TB_L_pow = mean(TB_L_ind);
%     TB_R_pow = mean(TB_R_ind);
%     TB_pow = mean([TB_L_ind TB_R_ind]);
%     TB_L_angle = mean(dangle(TB_L));
%     TB_R_angle = mean(abs(dangle(TB_R)));
%     TB_angle = mean(abs(dangle(TB)));
%
%
%     final_h = figure('Position',[200 50 500 900]);
%     subplot(2,1,1);
%     bar(bin_x,hist_dangle_norm)
%     set(gca,'tickdir','out')
%     set(gca, 'XLim', [ -pi pi ])
%     IT_num_str = num2str(IT_num);
%     IT_angle_str = num2str(IT_angle);
%     T_freq_str = num2str(T_freq);
%     rsq_str = num2str(rsq);
%     if length(IT_num_str)<5
%         IT_num_str = [IT_num_str '00000'];
%     end
%
%     if length(IT_angle_str)<5
%         T_angle_str = [T_angle_str '00000'];
%     end
%
%     if length(T_freq_str)<5
%         T_freq_str = [T_freq_str '00000'];
%     end
%     if length(rsq_str)<8
%         rsq_str = [rsq_str '00000'];
%     end
%     title([FileName(1:17) ':' 'IT-num = ' IT_num_str(1:5) ', T-freq = '  T_freq_str(1:5)],'interpreter','none')
%     subplot(2,1,2);
%     plot(dangle_nobig(1:end-1),dangle_nobig(2:end),'o');
%     grid on
%     set(gca, 'XLim', [ -2 2 ], 'YLim', [ -2 2 ] );
%     title(['R2 = '  rsq_str(1:8)],'interpreter','none')
%     set(final_h,'PaperPositionMode','auto');
%     tif_name = [outDir,savename,'.tif'];
%     saveas(final_h,fullfile(tif_name),'tif');
% end
%
% %% this part is for evaluating the turn detection
% if 1
%     rangename = ['region' num2str(5)];
%     turn.(rangename).dangle;
%     turn_ind = 50;
%     start_ind = turn_ind;
%     st_vid = peaks(turn.(rangename).T_id(turn_ind))-5;
%     track = figure('position',[900 100 800 800]);
%     sample_fr = 200;
%     color_dangle = (turn.(rangename).dangle)'/max(abs(turn.(rangename).dangle));
%     if 1 % whether to save the video
%         ff.FrameRate = 1;
%         movie_folder = [outDir,savename,'_',num2str(ff.FrameRate),'.mpeg'];
%         ff=VideoWriter(movie_folder,'MPEG-4');
%
%         open(ff);
%         ha=axes;
%     end
%     for i = st_vid:(st_vid + sample_fr)
%         imagesc(IM(:,:,i),[0 255]);
%         colormap(gray);
%         title(i)
%         drawnow
%         hold on;
%         plot(y(st_vid:i),x(st_vid:i),'k');
%         drawnow
%         hold on;
%         c = 2+abs(turn.(rangename).dangle (turn.(rangename).T_id(start_ind:turn_ind)))*5;
%
%
%         %         c = 2+5*min([abs(dangle(2:turn_ind)); abs(dangle(2:turn_ind)+2*pi); abs(dangle(2:turn_ind)-2*pi)],[],1);
%
%         %           plot(  y(peaks(1:T_id(turn_ind))),x(peaks(1:T_id(turn_ind))),'y.')
%         scatter(y(peaks(turn.(rangename).T_id(start_ind:turn_ind))),x(peaks(turn.(rangename).T_id(start_ind:turn_ind))),c.^2,'r');
%
%         %         line([y(peaks(T_id(1:(turn_ind - 1))))], [x(peaks(T_id(1:(turn_ind - 1))))],'Color',[(1-color_dangle(T_id(turn_ind)))/2 0 (1+color_dangle(T_id(turn_ind)))/2]);
%         line([y(peaks(turn.(rangename).T_id((turn_ind - 1))))], [x(peaks(turn.(rangename).T_id((turn_ind - 1))))],'Color',[(1-color_dangle(turn.(rangename).T_id(turn_ind)))/2 0 (1+color_dangle(turn.(rangename).T_id(turn_ind)))/2]);
%         title(['frame' num2str(i) '; turn ' num2str(turn_ind)])
%         drawnow
%
%
%         hold off
%         if i == peaks(turn.(rangename).T_id(turn_ind)) & (turn_ind + 1)<length(turn.(rangename).T_id)
%             turn_ind = turn_ind + 1;
%         end
%         %hold on;
%         %line([y_aft, y_bef],[x_aft, x_bef]);
%
%         %waitforbuttonpress
%         if 0 % whether to save the video
%             fdata=getframe(gcf);
%             writeVideo(ff,fdata);
%         end
%
%     end
%     %     close(ff)
%     title(['frame' num2str(i) '; turn ' num2str(turn_ind)])
%     set(track,'PaperPositionMode','auto');
%     sample_time = num2str(sample_fr/fr_rate*60);
%     video_name = [outDir  savename  '-track-' sample_time 'secs''.tif']
%     saveas(track,fullfile(video_name),'tif');
% end
%
% %% this part is for evaluating the turn detection
% if 1
%     rangename = ['region' num2str(5)];
%     turn.(rangename).dangle;
%     turn_ind = 50;
%     start_ind = turn_ind;
%     st_vid = peaks(turn.(rangename).T_id(turn_ind))-5;
%     track = figure('position',[900 100 800 800]);
%     sample_fr = 200;
%     color_dangle = (turn.(rangename).dangle)'/max(abs(turn.(rangename).dangle));
%     if 1 % whether to save the video
%         ff.FrameRate = 1;
%         movie_folder = [outDir,savename,'_',num2str(ff.FrameRate),'.mpeg'];
%         ff=VideoWriter(movie_folder,'MPEG-4');
%
%         open(ff);
%         ha=axes;
%     end
%     for i = st_vid:(st_vid + sample_fr)
%         imagesc(IM(:,:,i),[0 255]);
%         colormap(gray);
%         title(i)
%         drawnow
%         hold on;
%         plot(y(st_vid:i),x(st_vid:i),'k');
%         drawnow
%         hold on;
%         c = 2+abs(turn.(rangename).dangle (turn.(rangename).T_id(start_ind:turn_ind)))*5;
%
%
%         %         c = 2+5*min([abs(dangle(2:turn_ind)); abs(dangle(2:turn_ind)+2*pi); abs(dangle(2:turn_ind)-2*pi)],[],1);
%
%         %           plot(  y(peaks(1:T_id(turn_ind))),x(peaks(1:T_id(turn_ind))),'y.')
%         scatter(y(peaks(turn.(rangename).T_id(start_ind:turn_ind))),x(peaks(turn.(rangename).T_id(start_ind:turn_ind))),c.^2,'r');
%
%         %         line([y(peaks(T_id(1:(turn_ind - 1))))], [x(peaks(T_id(1:(turn_ind - 1))))],'Color',[(1-color_dangle(T_id(turn_ind)))/2 0 (1+color_dangle(T_id(turn_ind)))/2]);
%         line([y(peaks(turn.(rangename).T_id((turn_ind - 1))))], [x(peaks(turn.(rangename).T_id((turn_ind - 1))))],'Color',[(1-color_dangle(turn.(rangename).T_id(turn_ind)))/2 0 (1+color_dangle(turn.(rangename).T_id(turn_ind)))/2]);
%         title(['frame' num2str(i) '; turn ' num2str(turn_ind)])
%         drawnow
%
%
%         hold off
%         if i == peaks(turn.(rangename).T_id(turn_ind)) & (turn_ind + 1)<length(turn.(rangename).T_id)
%             turn_ind = turn_ind + 1;
%         end
%         %hold on;
%         %line([y_aft, y_bef],[x_aft, x_bef]);
%
%         %waitforbuttonpress
%         if 0 % whether to save the video
%             fdata=getframe(gcf);
%             writeVideo(ff,fdata);
%         end
%
%     end
%     %     close(ff)
%     title(['frame' num2str(i) '; turn ' num2str(turn_ind)])
%     set(track,'PaperPositionMode','auto');
%     sample_time = num2str(sample_fr/fr_rate*60);
%     video_name = [outDir  savename  '-track-' sample_time 'secs''.tif']
%     saveas(track,fullfile(video_name),'tif');
% end
%
% %% for more detailed evaluation
%
% if 1
%
%     startFr = 10;
%     endFr = startFr +750;
%     anglefac1 =200;
%     anglefac2 = 0.475;
%     rangename = ['region1'];
%     dangle = turn.(rangename).dangle;
%
%     peaks_scatter = peaks(turn.(rangename).selecTurn);
%
%
%     trace_var.ref = ref;
%     trace_var.startFr = startFr;
%     trace_var.endFr = endFr;
%     trace_var.anglefac1 = anglefac1;
%     trace_var.anglefac2 = anglefac2;
%     trace_var.dangle = dangle;
%     trace_var.peaks = peaks_scatter;
%
%     savefast(fullfile([outDir,fName,'_trace_var.mat']),'trace_var');
%
% else
%     [FileName,PathName] = uigetfile('*.mishVid*','Select the mishVid');
%     [~,fName,~]=fileparts(FileName);
%     outDir=[PathName,fName,'\swims\'];
%     load([outDir,fName, '_trace_var.mat']);
%     load([outDir,fName, '_tracexy_flt.mat']);
%     load([outDir,fName, '_turn_ind.mat']);
%     rangename = ['region1'];
%
%     y = tracexy_flt(1,:);
%     x = tracexy_flt(2,:);
%     ref = trace_var.ref;
%     startFr = trace_var.startFr;
%     endFr = trace_var.endFr;
%     anglefac1 = trace_var.anglefac1;
%     anglefac2 = trace_var.anglefac2;
%     dangle = trace_var.dangle;
%     peaks_scatter = trace_var.peaks;
%
%
% end
%
%
% trace = figure('position',[900 100 800 800]);
% imagesc(ref);
% colormap(gray);
% hold on
% %     rectangle('Position',[x0-calculatingR,y0-calculatingR,2*calculatingR,2*calculatingR],'Curvature', [1 1]);
% %     rectangle('Position',[x0-calculatingR2,y0-calculatingR2,2*calculatingR2,2*calculatingR2],'Curvature', [1 1]);
% %     rectangle('Position',[x0-calculatingR3,y0-calculatingR3,2*calculatingR3,2*calculatingR3],'Curvature', [1 1]);
%
%
% %     for i=600:1180   %080607_01
% % for i=40:368   %080607_04
% %          for i=320:560   %101509_01
% %              for i=1:229   %101509_04
% for i = startFr: min(endFr,length(peaks_scatter))
%     %         plot(y(peaks(i:(i+1))),x(peaks(i:(i+1))),'-g');
%     if (dangle(i)==0)
%         %         plot(y(peaks(i)),x(peaks(i)),'k.');
%     elseif(dangle(i)>0)
%         plot(y(peaks_scatter(i)),x(peaks_scatter(i)),'bo','MarkerSize',(abs(dangle(i))*anglefac1).^(anglefac2),'MarkerEdgeColor',[96,187,77]/255,'MarkerFaceColor',[96,187,77]/255);
%     else
%         plot(y(peaks_scatter(i)),x(peaks_scatter(i)),'ro','MarkerSize',(abs(dangle(i))*anglefac1).^(anglefac2),'MarkerEdgeColor',[172,84,152]/255,'MarkerFaceColor',[172,84,152]/255);
%     end
%
%     %         [ i, peaks(i), angle(i),  dangle(i)]
%     %         waitforbuttonpress;
% end
%
% IT_num_str= num2str(turn.(rangename).IT_num);
% IT_angle_str= num2str(turn.(rangename).IT_angle);
% IT_freq_str= num2str(turn.(rangename).T_freq);
%
% titlestr = sprintf('IT-num =  %s ; \n   IT_angle =  %s ; \n  T-freq = %s \n %s',  IT_num_str,IT_angle_str,IT_freq_str,FileName(1:14));
% title(titlestr,'interpreter','none')
%
% %%
% % saveas(trace,fullfile([outDir,savename,'_trace_',num2str(startFr),'to',num2str(endFr),'scale',num2str(anglefac1),'_',num2str(anglefac2),'.tif']),'tif');
% saveas(trace,fullfile([outDir,savename,'_trace_',num2str(startFr),'to',num2str(endFr),'scale',num2str(anglefac1),'_',num2str(anglefac2),'.fig']),'fig');
% exportfig(trace,[outDir,savename,'_trace_',num2str(startFr),'to',num2str(endFr),'scale',num2str(anglefac1),'_',num2str(anglefac2),'.eps'],'Color' ,'rgb');
% % saveas(trace,fullfile([outDir,savename,'_trace_',num2str(startFr),'to',num2str(endFr),'scale',num2str(anglefac1),'_',num2str(anglefac2),'.eps']),'epsc');
%
% return
%
%%
% for ii = 1:(length(dangle)/2)
%     qq = dangle(ii+1:end) .* dangle(1:end-ii);
%     a1(ii) = sum(qq > 0);
%     a2(ii) = sum(qq < 0);
% end
% n_cor = figure('position',[100, 100,  1000,800]);
% title('corelation with n turns agao ')
% hold on;
% subplot(2,1,1);plot(a1./(a1+a2));
% subplot(2,1,2);plot(a1./(a1+a2));
% xlim([0 10]);
% core_name = [outDir,savename,'-n-correlation-' '.tif']
% saveas(n_cor,fullfile(core_name),'tif');
%
% return
% %%
%
%
% figure
% for i = 1:size(IM,3)
%     hold on
%     plot(x(1:i),y(1:i))
%     if sum(peaks==i) > 0
%         plot(x(i),y(i),'ro')
%         ind = find(peaks==i);
%         title(num2str(angle(ind)))
%         waitforbuttonpress
%     else
%         title('')
%     end
% end
%
%
%
% return
%
%
%
% %h = fopen('4sept11_3.mishVid_640_480');
% clear all
%
% [FileName,PathName] = uigetfile('*.mishVid_640_480','Select file');
% filename = [PathName FileName];
%
% h = fopen([filename]);
% X = fread(h,'uint8');
% fclose(h)
%
% a = repmat((1:640)', 1,480)/640;
%
% figure
% q = zeros(1,640*480);
% i = 0;
% ang = zeros(1,100000);
% fr = zeros(1,100000);
% while length(q) > 0
%     i = i + 1;
%     fr(i) = mean(q(:));
%     q = reshape(q, 640, 480);
%     %    imagesc(q); colormap(gray)
%     %    waitforbuttonpress
%     q = mean(q,2);
%     q = bandpassmu(q,640,3, 50);
%     %    plot(q)
%     %    waitforbuttonpress
%     l = find(q==min(q));
%     l = l(1);
%     %    l = mean(mean(a.*q));
%     ang(i) = l;
%     %    plot(l, 0, 'ro')
%     %    xlim([0 640]);
%     %    ylim([-1 1])
%     %    imagesc(q,[0 255])
%     %    colormap(gray)
%     %    drawnow
%     q = fread(h, 640*480);
%     if mod(i,1000) == 0
%         i
%     end
% end
%
% ang = ang(2:i);
% fr = fr(2:i);
%
% fclose(h)
%
% figure; hold on
% %plot(-0.7+fr / max(abs(fr)),'k')
%
%
% st_ = find(diff((fr>0.2))==1);
% en_ = find(diff((fr>0.2))==-1);
%
% fr = smooth(fr,10);
% fr = (fr-min(fr))/(max(fr)-min(fr));
% fr = (fr>0.3);
% ang = (ang-min(ang))/(max(ang)-min(ang));
% ang = (ang-0.08)/2;
%
% mn_ = 0;
% mx_ = 1;
% for i = 1:length(st_)
%     h = patch([st_(i) st_(i) en_(i) en_(i) st_(i)], [mn_ mx_ mx_ mn_ mn_], [0.7 0.7 0.7]);
%     set(h,'edgecolor',[.7 .7 .7])
% end
%
% plot(ang,'k','linewidth',2)
% ylim([.15 .45])
% xlim([1 length(ang)])
% set(gca,'xtick',[1:120*10:length(ang)],'xticklabel',[0:120*10:length(ang)-1]/120)
% xlabel('time (s)')
% set(gca,'ytick',[])
% ylabel('tail movement')
%
% set(gcf,'paperposition',[0 1.6 4 1.6])
%
% print -depsc fig_electroShockDCN.eps
%
%
% %% plot the color coded direction
% turn_ind = 100;%685;
% start_ind = turn_ind;
% st_vid = peaks(T_id(turn_ind))-5;
% sample_fr = 28000;
% color_dangle = dangle'/max(abs(dangle));
%
% track = figure('position',[900 100 800 800]);
% imagesc(im,[0 255]);
% colormap(gray);
% hold on
% for i = st_vid:(st_vid + sample_fr)
%     %     scatter(y(peaks(T_id(start_ind:turn_ind))),x(peaks(T_id(start_ind:turn_ind))),c.^2,'k');
%     line([y(peaks(T_id((turn_ind - 1):turn_ind)))], [x(peaks(T_id((turn_ind - 1):turn_ind)))],'LineWidth',2,'Color',[(1-color_dangle(T_id(turn_ind)))/2 0.2 (1+color_dangle(T_id(turn_ind)))/2]);
%     title(['frame' num2str(i) '; turn ' num2str(turn_ind)],'interpreter','none')
%
%     if i == peaks(T_id(turn_ind)) & (turn_ind + 1)<=start_ind+120
%         turn_ind = turn_ind + 1;
%     elseif (turn_ind + 1)>start_ind+120
%         break
%     end
%
% end
% %%
% switches = dangle; %signed_turn_amp is a vector containing a dtring of -1's and 1's for left and right turns
% switches(dangle>0) = 1;
% switches(dangle<0) = -1;
% switches = find(diff(switches)~=0)+1;
%
% switches = switches(find(abs(dangle(switches))>15/180*pi)); %angthresh is used to throw out turns below a certain threshold. I use 0.
% % turn_amp = signed_turn_amp; %unnecessary, a relic of a previous analysis. But I ws too lazy to change the code in the loop below
%
% totsize = 10
% switches = switches(switches<=length(dangle)-totsize); % I only consider switches that are within 'totsize' from the end of the turn vector (this matters more for me %when I crop out swim events past a certain radius). Totsize is the number of turns into the future you want to plot.
%
% totcum = zeros(length(switches),totsize); %totsize is the number of turns into the future you want to plot
% for wtf = 1:length(switches);
%     if dangle(switches(wtf)) > 0
%         totcum(wtf,:) = cumsum(dangle(switches(wtf):switches(wtf)+totsize-1))-dangle(switches(wtf));
%     else
%         totcum(wtf,:) = -cumsum(dangle(switches(wtf):switches(wtf)+totsize-1))+dangle(switches(wtf));
%     end
% end
%
% figure;plot(mean(totcum))
% %mean(totcum) is your trace



