% 
% %% Plot tail curvatures atop moving fish
% trl = 2;
% pauseDur = 0.1;
% inds = (trl-1)*750 + 1: (trl-1)*750 + 750;
% var1 = IM_proc_crop(:,:,inds);
% var2 = repmat([71 71],numel(inds),1);
% mlInds2 = midlineInds(inds);
% tC = tailCurv(:,:,inds);
% 
% PlayFishTracking2(IM_proc_crop,'fishPos',[71 71], 'midlineInds',midlineInds,'tailAngles',tA,'plotCurv',1,'tailCurv',tailCurv,'frameInds',inds,'pauseDur',pauseDur);


%%
cd 'S:\Avinash\Tracking demo'

readMode =  'fromImages';
% readMode = 'fromMishVid';

poolSize  = 10;
switch readMode
    case 'fromMishVid'
        fName_prefix = input('Enter fish name, e.g., Fish1: ','s');
        [IM, outDir] = ReadMishVid();
        imgInds = 1:size(IM,3);
    case 'fromImages'
        imgDir = input('Enter image dir path:  ', 's')
        imgExt = input('Enter image extension, e.g. jpg:  ','s')
        imgInds = input('Enter indices of images to read as a vector:  ');
%         fName_prefix = input('Enter fish name, e.g., Fish1: ','s');
%         bpMessg =['Filter bandpass for finding fish pos, e.g. [15 25]. ' ...
%             'Press enter [] to skip (filtering is recommended if using collimated ' ...
%             'light during behavior): '];
%         fps = input('Enter frame rate (frames/sec): ');
%         nFramesInTrl  = input('# of frames in each trl (default = 750): ');
%         bp = input(bpMessg);        
%         outDir = fullfile(imgDir,'proc');
        IM = ReadImgSequence(imgDir,imgExt,imgInds);
end


%%
tic
ref = max(IM,[],3);
IM_proc = IM - repmat(ref,[1,1,size(IM,3)]);
toc
trl = {};
trl{1} = IM_proc(:,:,1:750);
trl{2} = IM_proc(:,:,751:1500);
% [fp{1},~] = GetFishPos(trl{1}, 25,'filter',bp,'process','parallel','lineLen',25);
% [fp{2},~] = GetFishPos(trl{2}, 25,'filter',bp,'process','parallel','lineLen',25);


%% Writing sequences
outDir = 'c:\users\pujalaa\documents\outDir'
procDir = fullfile(imgDir,'trl1');
mkdir(procDir);
disp('Writing processed images...')
tic
trl{1} = Standardize(trl{1});
for jj = 1:size(trl{1},3)
    fName = ['Img_' sprintf('%0.4d', jj) '.bmp'];    
    imwrite(trl{1}(:,:,jj),fullfile(procDir,fName),'bmp')
    disp(num2str(jj))
end
toc


procDir = fullfile(imgDir,'trl2');
mkdir(procDir);
disp('Writing processed images...')
tic
trl{2} = Standardize(trl{2});
for jj = 1:size(trl{2},3)
    fName = ['Img_' sprintf('%0.4d', jj) '.bmp'];    
    imwrite(trl{2}(:,:,jj),fullfile(procDir,fName),'bmp')
    disp(num2str(jj))
end
toc