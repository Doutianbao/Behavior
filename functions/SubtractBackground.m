function varargout = SubtractBackground(varargin)
%SubtractBackground Subtracts background from images
% IM_proc = SubtractBackground(IM)
% [IM_proc, refImg] = SubtractBackground(IM,refFrames);
% Inputs:
% IM - Image stack where 3rd time is time
% refFrames - A vector specifying the frames to average for background
%   subtraction. This is useful for considering only frames in which fish
%   moves
% Outputs:
% IM_proc - Stack with background subtracted images
% refImg = Reference image (median of image stack) used for background
%   subtraction
%
% Avinash Pujala, Koyama lab/HHMI, 2016

poolSize = 10;
IM = varargin{1};
disp('Computing ref frame...')
refFrames = 2:size(IM,3); % For some reason the 1st image doesn't record properly [AP - 20160621]
if (nargin == 2)
    refFrames = varargin{2};
end
tic
% im = median(IM(:,:,refFrames),3);
% im = max(IM(:,:,refFrames),[],3);
ref = mean(IM(:,:,refFrames),3);
toc

disp('Subtracting background...')
if size(IM,3)>=500
    IM_proc = ProcInParallel(IM,ref,poolSize);
else
    IM_proc = ProcInSerial(IM,ref);
end

varargout{1} = IM_proc;
varargout{2} = ref;

end
function IM_proc = ProcInParallel(IM,ref, poolSize)
if strcmpi(class(IM),'mappedTensor')
    IM_proc = MappedTensor(size(IM));
else
    IM_proc = zeros(size(IM));
end

imgFrames= 1:size(IM,3);
poolOpened = 0;
if matlabpool('size')==0
    matlabpool(poolSize)
    poolOpened = 1;
end
dispChunk = round(numel(imgFrames)/30);
parfor jj=imgFrames
    %     IM_proc(:,:,jj) = conv2(-squeeze(IM(:,:,jj))+ref,ones(5)/25,'same');
    IM_proc(:,:,jj) = -squeeze(IM(:,:,jj))+ref;
    if mod(jj, dispChunk)==0
        disp(['Img# ' num2str(jj)])
    end
end
if poolOpened
    matlabpool close
end
end

function IM_proc = ProcInSerial(IM,ref)
if strcmpi(class(IM),'mappedTensor')
    IM_proc = MappedTensor(size(IM));
else
    IM_proc = zeros(size(IM));
end

imgFrames= 1:size(IM,3);
dispChunk = round(numel(imgFrames)/30);
for jj=imgFrames
    %     IM_proc(:,:,jj) = conv2(-squeeze(IM(:,:,jj))+ref,ones(5)/25,'same');
    IM_proc(:,:,jj) = -squeeze(IM(:,:,jj))+ref;
    if mod(jj, dispChunk)==0
        disp(['Img# ' num2str(jj)])
    end
end
end