
% Unpacks grpData which contains the data from all ablation groups and returns 
% creates W, which contains all the wavelet transforms in a more easily
% accessible format

grps = fieldnames(grpData);
stimTypes = {'vib','dark'};
W = struct;
for hh = 1:length(grps)
    for ii = 1:length(stimTypes)
        ctrl = grpData.(grps{hh}).ctrl.(stimTypes{ii}).procData;
        abl = grpData.(grps{hh}).abl.(stimTypes{ii}).procData;
        disp([grps{hh} ', ctrl, ' stimTypes{ii}])
        for jj = 1:length(ctrl)
            blah = ctrl{jj}.W;
            blah = blah.curv.coeff;
            for kk = 1:length(blah)
                if jj == 1 && kk == 1
                    W.ctrl.(grps{hh}).(stimTypes{ii}) = cat(3,[],blah{kk});
                else
                    W.ctrl.(grps{hh}).(stimTypes{ii}) = cat(3, W.ctrl.(grps{hh}).(stimTypes{ii}),blah{kk});
                end
            end
        end
        disp([grps{hh} ', abl, ' stimTypes{ii}])
        for jj = 1:length(abl)
            blah = abl{jj}.W;
            blah = blah.curv.coeff;
            for kk = 1:length(blah)
                if jj == 1
                    W.abl.(grps{hh}).(stimTypes{ii}) = cat(3,[],blah{kk});
                else
                    W.abl.(grps{hh}).(stimTypes{ii}) = cat(3, W.abl.(grps{hh}).(stimTypes{ii}),blah{kk});
                end
            end
        end
    end
end
