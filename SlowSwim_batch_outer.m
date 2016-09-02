
fishDir = 'S:\Avinash\Ablations and behavior\Intermediate RS\20160901';

%% Getting the names of fish image dirs
disp('Getting the names of fish dirs...')
tic
blah = dir(fishDir);
imgDirs = {};
count = 0;
for jj = 1:length(fldrs)
    nm = blah(jj).name;
    if isempty(strfind(nm,'.'))
        count = count + 1;
        imgDirs{count} = fullfile(fishDir,nm);
    end
end
toc

%% Go into fast dir
for fishNum = 1:length(imgDirs)
    fishName = ['f' num2str(fishNum)];
    disp(['Fish # ' num2str(fishNum)]);
    blah  = dir(imgDirs{fishNum});
    for jj = 1:length(blah)
        nm = blah(jj).name;
        if isempty(strfind(lower(nm),'.')) && ~isempty(strfind(lower(nm),'fast'))
           path = fullfile(imgDirs{fishNum},nm);
           disp(nm)
           blah = dir(path);
           for kk = 1:length(blah)
               nm = blah(kk).name;               
               if isempty(strfind(lower(nm),'.'))
                   imgDir = fullfile(path,nm);
                   SlowSwim_batch
               end
           end
        end
    end
end
