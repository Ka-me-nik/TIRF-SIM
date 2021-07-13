% saveTiff - helper function for saving tiff movies
%
% Parameters:
%    I - image data to save (3D or 4D)
%    fn - filename of the output tiff movie
%    optional tag, value pairs - define values for custom tiff tags (1 value for each frame)
%
% Supports grayscale or RGB (3-channel) data.
function saveTiff(I,fn,varargin)

ts.ImageLength = size(I,1);
ts.ImageWidth = size(I,2);
rgb = size(I,3)==3;
if size(I,4)==1 && ~rgb
    I = permute(I,[1,2,4,3]);
end
if rgb
    ts.Photometric = Tiff.Photometric.RGB;
    ts.SamplesPerPixel = 3;
else
    ts.Photometric = Tiff.Photometric.MinIsBlack;
    ts.SamplesPerPixel = 1;
end
switch class(I)
    case 'uint8'
        ts.SampleFormat = Tiff.SampleFormat.UInt;
        ts.BitsPerSample = 8;
    case 'uint16'
        ts.SampleFormat = Tiff.SampleFormat.UInt;
        ts.BitsPerSample = 16;
    case 'single'
        ts.SampleFormat = Tiff.SampleFormat.IEEEFP;
        ts.BitsPerSample = 32;
    otherwise
        error(['Class ''' class(I) '''currently not supported.']);
end
ts.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

N = size(I,4);
assert(mod(nargin,2)==0,'Incorrect number of inputs.');
tag = {varargin{1:2:end}};
% assert(all(cellfun(@ischar,tag)),'Incorrect tag parameter.');
val = {varargin{2:2:end}};
assert(all(cellfun(@(x) length(x)==N,val)),'Incorrect tag value parameter size.');

T = Tiff(fn,'w');
for i = 1:N
    if i>1
        writeDirectory(T);
    end
    setTag(T,ts);
    for t = 1:length(tag)
        setTag(T,tag{t},val{t}(i));
    end
    write(T,I(:,:,:,i));
end
