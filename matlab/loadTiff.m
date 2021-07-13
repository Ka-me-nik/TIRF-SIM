% loadTiff - convenience function for loading tiff movies in one go
%
% Filename is specified in the 'fn' parameter.
%
% Supports multiple channels.
function [I,info] = loadTiff(fn)

info = imfinfo(fn);
I = repmat(imread(fn,'Info',info),1,1,1,length(info));
for i = 2:length(info)
    I(:,:,:,i) = imread(fn,'Index',i,'Info',info);
end
I = squeeze(I);
