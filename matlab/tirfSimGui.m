% tirfSimGui - GUI for visualisation of TIRF-SIM data and tracks from CME analysis
%
% Parameters:
%    folder - cell data folder ['' = will ask for the folder]
%    boost - intensity scaling for TIRF-SIM data (hack for visual improvement) [4]
%
% Shows TIRF-SIM data as well as results of "classical CME analysis" and
% enables manual selection of "interesting" (in any sense) tracks for
% further processing.
function tirfSimGui(folder,boost)

marker = 'x';  % displayed marker for CME track positions
rCCP = 10;     % cut-out radius (size = 2*r+1)
simIdx = 0;    % 0 = mean, 1-9 = SIM index
if nargin<2
    boost = 4; % intensity scale
end

if nargin<1 || isempty(folder)
    folder = uigetdir(pwd, 'Select the ''cell'' folder:');
    if folder==0
        return
    end
end
if ~strcmp(folder(end), filesep)
    folder = [folder filesep];
end
% detect framerate
fr = regexp(folder, '_(\d+)?(\.)?\d+s', 'match');
if ~isempty(fr)
    framerate = str2double(fr{1}(2:end-1));
else
    fr = regexp(folder, '_\d+ms', 'match');
    if ~isempty(fr)
        framerate = str2double(fr{1}(2:end-2))/1000;
    else
        framerate = 2; % default: 2s
    end
end

fn = dir([folder 'TIRF*_L.mrc']);
if numel(fn)<1
    error('Cannot find SIM reconstruction.');
end
wl = cellfun(@(s) str2double(s(5:7)),{fn.name});
if wl(1)~=488
    error('First (master) wavelength is not 488.');
end
nCh = numel(wl);
N = 0;
imgsz = [];
fprintf('Loading image data ... ');
I = cell(2,nCh);
for c = 1:nCh
    a = ReadMRC([folder fn(c).name(1:end-6) '.mrc']);
    if c==1
        N = size(a,3)/9;
        imgsz = size(a(:,:,1));
    else
        if any(size(a)~=[imgsz,9*N])
            error('Inconsistent input data size.\n');
        end
    end
    a = single(reshape(a,size(a,1),size(a,2),9,[]));
    if simIdx < 1
        I{1,c} = rot90(squeeze(mean(a,3)));
    else
        I{1,c} = rot90(squeeze(a(:,:,simIdx,:)));
    end
end
clear a
fnr = dir([folder 'TIRF*_L*reg.tif']);
if numel(fnr)>0
    wlr = cellfun(@(s) str2double(s(5:7)),{fnr.name});
    i = arrayfun(@(w) find(wl==w),wlr);
    fn(i) = fnr;
end
fn = {fn.name};
BC = zeros(N,nCh);
for c = 1:nCh
    if strcmpi(fn{c}(end-2:end),'mrc')
        I{2,c} = rot90(ReadMRC([folder fn{c}]));
    else
        I{2,c} = loadTiff([folder fn{c}]);
    end
    if any(size(I{2,c})~=[2*imgsz,N])
        erorr('Inconsistent reconstruction data sizes.');
    end
    b = double(squeeze(mean(mean(I{2,c},1),2)));
    x = (1:N).';
    ft = fit(x,b,'exp2');
    b = ft.a*exp(ft.b*x)+ft.c*exp(ft.d*x);
    BC(:,c) = mean(b)./b; % preserve mean value
%     BC(:,c) = 1./b; % normalize mean to 1
end
movies = fn;
m = zeros(N,nCh,2);
for i = 1:nCh
    m(:,i,1) = squeeze(min(min(I{1,i},[],1),[],2));
    m(:,i,2) = squeeze(max(max(I{1,i},[],1),[],2));
end
% mI2 = [min(I{2,1}(:)),min(I{2,2}(:))];
fprintf('done.\n');
idx = 1;

% fn = [data.source 'Detection' filesep 'detection_v2.mat'];
% if exist(fn,'file')~=2
%     error('Detection data not found.');
% end
% D = cell(1,N);
% det = load(fn);
% for i = 1:N
%     D{i} = [2*det.frameInfo(i).x;2*det.frameInfo(i).y;det.frameInfo(i).A].';
% end

msk = imread([folder 'Detection' filesep 'cellmask.tif']);

%load tracks
fprintf('Loading tracks ... ');
fn = [folder 'Tracking' filesep 'ProcessedTracks.mat'];
if exist(fn,'file')~=2
    error('Tracking data not found.');
end
trk = load(fn);
trk = trk.tracks;
trk = trk([trk.nSeg]==1);
[~,i] = sort([trk.lifetime_s],'descend');
trk = trk(i);
nt = numel(trk);
% apply cell mask
x = NaN(1,nt);
y = NaN(1,nt);
for t = 1:nt
    trk(t).x = trk(t).x(1,:);
    trk(t).y = trk(t).y(1,:);
    x(t) = round(nanmean(trk(t).x));
    y(t) = round(nanmean(trk(t).y));
end
trk = trk(msk(sub2ind(size(msk), y, x))==1);
trk = trk([trk.lifetime_s] >= 5*framerate);
nt = numel(trk);
T = cell(1,N);
for i = 1:N
    tt = [];
    for t = 1:nt
        ii = find(trk(t).f==i);
        if sum(ii)>0
            tt = [tt [2*trk(t).x(ii);2*trk(t).y(ii);t*ones(size(ii))]];
        end
    end
    T{i} = tt.';
%     T{i} = [2*det.frameInfo(i).x;2*det.frameInfo(i).y;det.frameInfo(i).A].';
end
fprintf('done.\n');
idxT = 0;

fn = [folder 'Tracking' filesep 'exportedTracks.mat'];
if exist(fn,'file')==2
    expTrk = load(fn);
    expTrk = expTrk.tracks;
else
    z = cellfun(@(x)zeros(size(x)),{trk.f},'Uniform',0);
    expTrk = struct('f',{trk.f},'x',cellfun(@(x)2*x,{trk.x},'Uniform',0),'y',cellfun(@(x)2*x,{trk.y},'Uniform',0),'cx',z,'cy',z,'status',num2cell(zeros(size(z))));
    et = dir([folder 'Tracking' filesep '*Track_*_' num2str(wl(1)) '.tif']);
    en = cellfun(@(x)str2double(regexpi(x,'(?<=_)\d+(?=_\d{3})','match','once')),{et.name});
    [expTrk(en).status] = deal(1);
end

S = cell(1,2);
chST = cell(1,2);
chCh = cell(nCh,2);
img = cell(1,2);
ax = zeros(1,2);
f = figure('Name',sprintf('TIRF-SIM'),'Units','normalized','Position',[.04 .04 .9 .9],'Menubar','none','Toolbar','figure','WindowKeyPressFcn',@keyPress,'CloseRequestFcn',@closeFig);
if ~verLessThan('matlab','9.4')
    set(f,'WindowState','maximized');
end
for i = 1:2
    chST{i} = uicontrol(f,'Style','checkbox','String',['Show track points (' num2str(i) ')'],'Units','normalized','Position',[0.02+(2-i)*0.49,0.02,0.2,0.02],'Callback',@chSDChange);
    for j = 1:nCh
        chCh{j,i} = uicontrol(f,'Style','checkbox','String',[num2str(wl(j)) ' (' num2str(i) ')'],'Units','normalized','Position',[0.1+j*0.1+(2-i)*0.49,0.02,0.2,0.02],'Callback',@chSDChange,'Value',1);
    end
    S{i} = subplot('Position',[0.02+(2-i)*0.49, 0.04, 0.47, 0.96]);
    id = imresize(I{i,1}(:,:,idx),3-i);
    id = id - min(id(:));
    id = id/max(id(:));
    img{i} = image(repmat(id,[1,1,3]));
    set(img{i},'ButtonDownFcn',@imgClick);
    axis equal off
    ax(i) = img{i}.Parent;
    zoom reset
    hold on
end
linkaxes(ax,'xy');
txt = uicontrol(f,'Style','text','String',['Frame: 1 / ' num2str(N)],'Units','normalized','Position',[0 0 .05 .02]);
sld = uicontrol(f,'Style','slider','Min',1,'Max',N,'SliderStep',[1 1]./(N-1),'Value',idx,'Units','normalized','Position',[.05 0 .95 .02],'Callback',@sldChange);
txtT = uicontrol(f,'Style','text','String','Track:','Units','normalized','Position',[0.4 .985 0.2 .015]);
uicontrol(f,'Style','pushbutton','String','Export','Units','normalized','Position',[0.6 .98 0.1 .02],'Callback',@export);
uicontrol(sld);
L = cell(1,4);
sldChange(sld);
% waitfor(f); % wait till the form is closed

function sldChange(obj,~)
    idx = round(get(obj,'Value'));
    set(txt,'String',sprintf('Frame: %i / %i',idx,N));
    idata = zeros([2*imgsz,3,2]);
    for k = 1:nCh
        if chCh{k,1}.Value>0
            tmp = imresize(I{1,k}(:,:,idx),2);
            idata(:,:,4-k,1) = (tmp-m(idx,k,1))/(m(idx,k,2)-m(idx,k,1));
        end
        if chCh{k,2}.Value>0
            tmp = I{2,k}(:,:,idx);
            if boost==0
                tmp = tmp-min(tmp(:));
                tmp = tmp/max(tmp(:));
            else
                tmp = (boost*tmp-m(idx,k,1))/(m(idx,k,2)-m(idx,k,1));
            end
            idata(:,:,4-k,2) = tmp;
        end
    end
    if nCh==1
        idata(:,:,1:2,:) = repmat(idata(:,:,3,:),[1,1,2,1]);
    end
    img{1}.CData = idata(:,:,:,1);
    img{2}.CData = idata(:,:,:,2);
    for k = 1:4
        if ~isempty(L{k})
            delete(L{k});
            L{k} = [];
        end
    end
    if ~isempty(T{idx})
        nat = T{idx}(T{idx}(:,3)~=idxT,:);
        at = T{idx}(T{idx}(:,3)==idxT,:);
        for k = 1:2
            subplot(S{k});
            if chST{k}.Value>0
                L{k} = scatter(nat(:,1),nat(:,2),100,[1 1 1],marker,'ButtonDownFcn',@imgClick);
                L{k}.MarkerEdgeAlpha = .5;
            end
            L{k+2} = scatter(at(:,1),at(:,2),100,[1 1 1],'s');
            L{k+2}.MarkerEdgeAlpha = .5;
        end
        if idxT>0
            set(txtT,'String',sprintf('Track: %i (%g / %g s)',idxT,(idx-trk(idxT).start+.5)*framerate,trk(idxT).lifetime_s));
        else
            set(txtT,'String','Track:');
        end
    end
end
function imgClick(obj,~)
    p = get(obj.Parent,'CurrentPoint');
    p = p(1,1:2);
    d = sqrt(sum((T{idx}(:,1:2)-p).^2,2));
    [d,di] = min(d);
    if (d<20)
        idxT = T{idx}(di,3);
    else
        idxT = 0;
    end
    chSDChange;
end
function chSDChange(~,~)
    sldChange(sld);
    uicontrol(sld);
end
function keyPress(~,e)
    switch e.Key
        case '1'
            chST{1}.Value = 1-chST{1}.Value;
            chSDChange;
        case '2'
            chST{2}.Value = 1-chST{2}.Value;
            chSDChange;
    end
end
function closeFig(~,~)
    ss.tracks = expTrk;
    for ci = 1:nCh
        ss.channels(ci).movie = movies{ci};
        ss.channels(ci).wavelength = wl(ci);
        ss.channels(ci).bleaching = BC(:,ci);
    end
    save([folder 'Tracking' filesep 'exportedTracks.mat'],'-struct','ss');
    delete(gcf);
end
function export(~,~)
    if idxT>0
        XY = -ones(trk(idxT).end-trk(idxT).start+1,2);
        for it = trk(idxT).start:trk(idxT).end
            itf = find(T{it}(:,3)==idxT);
            if ~isempty(itf)
                XY(it-trk(idxT).start+1,:) = round(T{it}(itf(1),1:2)+.5);
            end
        end
        E = cell(1,nCh);
        for ich = 1:nCh
            E{ich} = zeros([rCCP*2+1,rCCP*2+1,size(XY,1)],'single');
            for it = 1:size(XY,1)
                if XY(it,1)>=0
                    tmpI = padarray(I{2,ich}(:,:,it+trk(idxT).start-1),rCCP*[1,1]);
                    tmpI = tmpI * BC(it+trk(idxT).start-1,ich);
                    E{ich}(:,:,it) = tmpI(XY(it,2)+(0:2*rCCP),XY(it,1)+(0:2*rCCP));
                end
            end
            saveTiff(E{ich},[folder 'Tracking' filesep 'Track_' num2str(idxT) '_' num2str(wl(ich)) '.tif']);
        end
        expTrk(idxT).status = 1;
    end
end
end
