
data.abl.vib.mean = data.abl.vib_mean;
data.abl.vib.std = data.abl.vib_std;
data.abl.vib.cv = data.abl.vib_cv;
data.abl.vib.corrVec = data.abl.vib_corrVec;


data.abl.dark.mean = data.abl.dark_mean;
data.abl.dark.std = data.abl.dark_std;
data.abl.dark.cv = data.abl.dark_cv;
data.abl.dark.corrVec = data.abl.dark_corrVec;

delNames = {'vib_mean','vib_std','vib_cv', 'vib_corrVec','dark_mean','dark_std','dark_cv', 'dark_corrVec'};

data.abl = rmfield(ata.abl,delNames)