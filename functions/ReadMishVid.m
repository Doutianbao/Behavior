function varargout = ReadMishVid(imgDir,varargin)
% ReadMishVid Function for reading .mishVid files created by Misha's code
% IM = ReadMishVid(imgDir,varargin)
%      varargin = {fps, imgDims, duration}
% [IM,outDir] = ReadMishVid(...);
% Inputs:
% imgDir - Directory where the .mishVid is stored
% fps - Frames per second. frameRate = [] defaults to 30
% imgDims - [H, W], where H = image height, W = image width; default = [600, 600]
% duration - Video duration (in mins)
% Outputs:
% IM - Image stack of dims H x W x T, where T is time points
% outDir - Directory where variables are to be saved

fps = 30;
if nargin == 4
    fps = varargin{1};
    imgDims = varargin{2};
    duration = varargin{3};
elseif nargin == 3
    fps = varargin{1};
    imgDims = varargin{2};
    if isempty(imgDims)
        imgDims = input('Enter image dims as [height, width], e.g. [600,600] : ');
    end
    duration = input('Enter video duration in mins, e.g., 20 : ');
elseif nargin ==2
    fps= varargin{2};
    duration = input('Enter video duration in mins, e.g., 20 : ');
    imgDims = input('Enter image dims as [height, width], e.g. [600,600] : ');
elseif nargin == 0
    imgDir = input('Enter image directory : ', 's');
    duration = input('Enter video duration in mins, e.g., 20 : ');
    imgDims = input('Enter image dims as [height, width], e.g. [600,600] : ');
end

if isempty(fps)
    fps = 30;
end
cd(imgDir)
[fileName,pathName] = uigetfile('*.mishVid*','Select the mishVid');
dirFile = [pathName fileName];
[~,fname,~]=fileparts(fileName);
outDir=[pathName,fname,'\swims\'];
disp(fileName)
% imgDims = input('Enter image dimensions as [height width] e.g. [500 600] :  ');
h_axis = imgDims(1);
v_axis = imgDims(1);
imgCircShift = [0 0];
h = fopen(dirFile);

% Imaging parameters
start_f = 1;
dur_f =  duration*60*30;
samp_int = 1;
stop_f =  start_f + dur_f;
IM = zeros(h_axis,v_axis,ceil(dur_f/samp_int));

% Reading video frames
clear M
for i=start_f:samp_int:stop_f
    if mod(i,500)==0,disp(i),end
    fseek(h,h_axis*v_axis*i,'bof');
    x=fread(h,h_axis*v_axis,'uint8');
    if isempty(x)
        break
    else
        if length(x) == h_axis*v_axis
            xx=reshape(x(1:h_axis*v_axis),h_axis,v_axis);
            xx = circshift(xx,imgCircShift);
        else
            break;
        end
        if mod(i,10) == 0
            imagesc(xx,[0 255]);colormap(flipud(jet));  %%imaging is flipped: turn left is actually turn right!!! Yu 8/7/2014
            axis image
            title(i)
            drawnow
            shg
        end
        IM(:,:,1+floor(i/samp_int)) = xx;
    end
end
fr_num = floor(i/samp_int);
IM=IM(:,:,1:fr_num-1);
fclose(h)

varargout{1} = IM;
varargout{2} = outDir;
end

