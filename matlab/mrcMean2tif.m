% mrcMean2tif - prepare input data for CME Analysis package
%
% Parameters:
%    idx - which (of the 9) of the input SIM images to use [0 = mean of all 9]
%    mask - filter specific mrc movies ['**/TIRF488*L.mrc' = only 488 movies]
%
% Converts specified movies in mrc format to tiff while either using only
% one SIM image or average all 9 of them (based on 'idx' parameter).
function mrcMean2tif(idx,mask)

if nargin<1
    idx = 0;
end
if nargin<2
    mask = ['**' filesep 'TIRF488*L.mrc'];
end

rd = uigetdir(pwd, 'Select the ''condition'' folder:');
fprintf('Processing "%s":\n',rd);
d = dir([rd filesep mask]);
n = length(rd);
for i = 1:length(d)
    t = 0;
    dt = dir([d(i).folder filesep '*.txt']);
    if ~isempty(dt)
        f = fopen([dt(1).folder filesep dt(1).name],'r');
        for j = 1:12, fgetl(f); end
        t = fscanf(f,'%f %*s %*s %*i %*f %*f %*f %*f %*i %*f %i',[2 Inf]).';
        fclose(f);
        t(t(:,2)<1,:) = [];
        t(t(:,2)>=max(t(:,2)),:) = [];
        j = diff([0;t(:,2)])>0;
        t = t(j,1);
        t = mean(diff(t));
    end
    fi = [d(i).folder filesep d(i).name(1:end-6) '.mrc'];
    fprintf('(%2i/%2i) ...%s',i,length(d),fi(n+1:end));
    if t>0
        fprintf(' - average frame time: %g s', t);
    end
    fprintf('\n');
    if idx<1
        fo = [d(i).folder filesep d(i).name(1:end-6) '-mean.tif'];
    else
        fo = [d(i).folder filesep d(i).name(1:end-6) '-' num2str(idx) '.tif'];
    end
    if exist(fo,'file')==2
        fprintf('   ... already exists.\n');
        continue
    end
    a = ReadMRC(fi);
    fprintf('   ... %i images',size(a,3)/9);
    if mod(size(a,3),9)~=0
        fprintf(' = inconsistent size.\n');
        continue
    end
    if idx<1
        b = rot90(squeeze(mean(reshape(a,size(a,1),size(a,2),9,[]),3,'native')));
    else
        b = rot90(a(:,:,idx:9:end));
    end
    for j = 1:size(b,3)
        imwrite(b(:,:,j),fo,'WriteMode','append');
    end
    if t>0
%         movefile();
    end
    fprintf(' - done.\n');
end
