function IM_proc = ProcessImages(IM)
%ProcessImages Smooth and background subtract images
% IM_proc = ProcessImages(IM)
disp('Computing smoothed mean frame...')
tic
imKer = ones(5)/5^2;
im=conv2(squeeze(mean(IM,3)),imKer,'same');
toc

nWorkers= 10;
disp('Processing images...')
tic
IM_proc = zeros(size(IM));
imgFrames= 1:size(IM,3);
% matlabpool(nWorkers)
for jj=imgFrames
    ii=conv2(-squeeze(IM(:,:,jj))+im,ones(5)/25,'same');
    IM_proc(:,:,jj) = ii;
    if mod(jj,500)==0
        disp(['Img # ' num2str(jj)])
    end
end
% matlabpool close
toc

end

