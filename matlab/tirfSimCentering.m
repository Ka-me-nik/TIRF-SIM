% tirfSimCentering - GUI for final polishing of tracks exported by tirfSimGui
%
% Optional parameter:
%    Folder - cell data folder ['' = will ask for the folder]
%
% Named parameters:
%    Radius - cut-out radius (size = 2R+1) [10]
%    Marker - displayed marker for CME track positions ['x']
%    Color - cut-out overlay color [[.5 .5 1]]
function tirfSimCentering(varargin)

ip = inputParser;
ip.addOptional('Folder', '', @ischar);
ip.addParameter('Radius', 10, @(x)validateattributes(x,{'numeric'},{'positive','scalar','integer'}));
ip.addParameter('Marker', 'x', @ischar);
ip.addParameter('Color', [.5 .5 1], @(x) ischar(x) || (isnumeric(x) && length(x)==3));
ip.parse(varargin{:});
folder = ip.Results.Folder; % cell folder
rCCP = ip.Results.Radius;   % cut-out radius (size = 2*r+1)
marker = ip.Results.Marker; % displayed marker for CME track positions
cMark = ip.Results.Color;   % cut-out overlay color

rMark = 1000/65/2; % 1/4 um in TIRF-SIM resolution
keyBegin = 'b';
keyEnd = 'e';
% tags
tags = struct('name',{'split','ring'},'tag',{1,2},'key',{'x','o'},'marker',{'x','o'});

if isempty(folder)
    folder = uigetdir(pwd, 'Select the ''cell'' folder:');
    if folder==0
        return
    end
end
if ~strcmp(folder(end), filesep)
    folder = [folder filesep];
end

fn = [folder 'Tracking' filesep 'exportedTracks.mat'];
assert(exist(fn,'file')==2,'Exported tracks data not found. Run ''tirfSimGui'' first.');
expTrk = load(fn);
chnl = expTrk.channels;
nCh = length(chnl);
ID = find([expTrk.tracks.status]>0);
assert(~isempty(ID),'No exported tracks. Use ''tirfSimGui'' to export some.');
if ~isfield(expTrk.tracks,'tag')
    for i = 1:length(expTrk.tracks)
        expTrk.tracks(i).tag = zeros(size(expTrk.tracks(i).f));
    end
end
trk = expTrk.tracks(ID);
changed = false(size(ID));
% track correction (remove gaps/splitting/merging)
err = find(cellfun(@(x)any(isnan(x)),{trk.f}));
if ~isempty(err)
    for i = err
        f = min(trk(i).f):max(trk(i).f);
        [uf,ii] = unique(trk(i).f);
        ii = ii(~isnan(uf)); % ignore the NaN
        [~,uf] = setdiff(f,uf);
        while ~isempty(uf)
            ii = ii([1:uf(1)-1,uf(1)-1:end]);
            uf(1) = [];
        end
        trk(i).f = f;
        trk(i).x = trk(i).x(ii);
        trk(i).y = trk(i).y(ii);
        trk(i).cx = trk(i).cx(ii);
        trk(i).cy = trk(i).cy(ii);
        trk(i).tag = trk(i).tag(ii);
        changed(i) = true;
    end
    fprintf('Corrupted tracks repaired: %i\n', length(err));
end
err = find(arrayfun(@(t)min(t.x+t.cx)<0,trk)|arrayfun(@(t)min(t.y+t.cy)<0,trk));
if ~isempty(err)
    for i = err
        ii = (trk(i).x+trk(i).cx)<0 | (trk(i).y+trk(i).cy)<0;
        trk(i).f(ii) = [];
        trk(i).x(ii) = [];
        trk(i).y(ii) = [];
        trk(i).cx(ii) = [];
        trk(i).cy(ii) = [];
        trk(i).tag(ii) = [];
        changed(i) = true;
    end
    fprintf('Invalid starts/ends repaired: %i\n', length(err));
end
% end track correction
idx = 1;
frm = cellfun(@(x)x(1),{trk.f});

fprintf('Loading data ...');
I = cell(1,nCh);
mi = zeros(nCh,2);
for i = 1:nCh
    if strcmpi(chnl(i).movie(end-2:end),'mrc')
        I{i} = rot90(ReadMRC([folder chnl(i).movie])).*reshape(chnl(i).bleaching,1,1,[]);
    else
        I{i} = loadTiff([folder chnl(i).movie]).*reshape(chnl(i).bleaching,1,1,[]);
    end
    mi(i,:) = [min(I{i}(:)),max(I{i}(:))];
    I{i} = padarray(I{i},rCCP*[1,1]);
end
N = size(I{1},3);
fprintf(' done.\n');

% prepare masks (circles) for radial function computation
C = zeros(2*rCCP+1,2*rCCP+1,2*rCCP);
[x,y] = meshgrid(-rCCP-.495:.01:rCCP+.495);
for r = 1:2*rCCP
    tc = (x.^2+y.^2<(r/2)^2);
    C(:,:,r) = reshape(mean(reshape(reshape(mean(reshape(tc,100,[])),2*rCCP+1,[]).',100,[])),size(C(:,:,1)));
end
for r = 2*rCCP:-1:2
    C(:,:,r) = C(:,:,r)-C(:,:,r-1);
end
C = C./sum(sum(C));

fprintf('Calculating intensities ...');
A = cell(1,length(trk));
tc = zeros([2*rCCP+1,2*rCCP+1,nCh]);
for i = 1:length(trk)
    for j = 1:length(trk(i).f)
        for c = 1:nCh
            tc(:,:,c) = I{c}(round(trk(i).y(j)+trk(i).cy(j))+(0:2*rCCP),round(trk(i).x(j)+trk(i).cx(j))+(0:2*rCCP),trk(i).f(j));
            A{i}(j,c) = max(sum(sum(tc(:,:,c).*C)));
        end
    end
end
fprintf(' done.\n');

f = figure('Name',sprintf('TIRF-SIM cut-out centering'),'Units','normalized','Position',[.04 .04 .9 .9],'Menubar','none','Toolbar','none','WindowKeyPressFcn',@keyPress,'CloseRequestFcn',@closeFig);
if ~verLessThan('matlab','9.4')
    f.WindowState = 'maximized';
end
subplot('Position',[0.01,0.03,0.63,0.96]);
img = image(I{1}(rCCP+1:end-rCCP,rCCP+1:end-rCCP,frm(idx)),'ButtonDownFcn',@imgClick);
colormap(gray(256));
set(gca,'CLim',mi(1,:));
axis equal off
zoom reset
fi = frm(idx)-trk(idx).f(1)+1;
rec = rectangle('Position',[trk(idx).x(fi)+trk(idx).cx(fi)-rCCP-0.5,trk(idx).y(fi)+trk(idx).cy(fi)-rCCP-0.5,2*rCCP+1,2*rCCP+1],'EdgeColor','g');
subplot('Position',[0.65,0.03,0.23,0.4]);
cut{1} = image(zeros([2*rCCP+1,2*rCCP+1,3]));
axis equal off
zoom reset
hold on
tfi = 0:pi/50:2*pi;
plot(rMark*cos(tfi)+rCCP+1,rMark*sin(tfi)+rCCP+1,':','Color',cMark);
plot(rMark/2*cos(tfi)+rCCP+1,rMark/2*sin(tfi)+rCCP+1,':','Color',cMark);
MD = plot(rCCP+1,rCCP+1,marker,'Color',cMark);
% plot(rCCP+1,rCCP+1,'o','MarkerSize',100,'Color',[.5 .5 1]);
line([.5,rCCP+1;2*rCCP+1.5,rCCP+1],[rCCP+1,.5;rCCP+1,2*rCCP+1.5],'Color',cMark,'LineStyle',':');
subplot('Position',[0.89,0.03,0.1,0.2]);
cut{2} = image(zeros([2*rCCP+1,2*rCCP+1,3]),'ButtonDownFcn',@cutButDown);
axis equal off
zoom reset
subplot('Position',[0.89,0.24,0.1,0.2]);
cut{3} = image(zeros([2*rCCP+1,2*rCCP+1,3]),'ButtonDownFcn',@cutButDown);
ci = [1,2,3];
axis equal off
zoom reset
ax(1) = subplot('Position',[0.65,0.455,0.23,0.2]);
hold on
RI = cell(1,2);
RI{1} = plot((1:2*rCCP)/2,zeros(1,2*rCCP),'g');
RI{2} = plot((1:2*rCCP)/2,zeros(1,2*rCCP),'r');
ax(2) = subplot('Position',[0.65,0.68,0.23,0.2],'ButtonDownFcn',@graphButDown);
hold on
TA = cell(1,2);
TA{1} = plot(1:N,zeros(1,N),'g');
TA{2} = plot(1:N,zeros(1,N),'r');
TA{3} = line([1 1],[0,1],'Color',[.5 .5 .5]);
for i = 1:length(tags)
    TA{i+3} = scatter([],[],50,'b',tags(i).marker);
end
lst = uicontrol(f,'Style','listbox','String',arrayfun(@(x)sprintf('%3i   (%3i - %3i)',ID(x),trk(x).f(1),trk(x).f(end)),1:length(ID),'Uniform',0),'Units','normalized','Position',[0.89 0.45 0.1 0.54],'Value',idx,'Callback',@lstChange);
txt = uicontrol(f,'Style','text','String',['Frame: ' num2str(frm(idx)) ' / ' num2str(N)],'Units','normalized','Position',[0 0 .05 .02]);
tagTxt = [];
for i = 1:length(tags)
    tagTxt = sprintf('%s, %s - %s',tagTxt,tags(i).key,tags(i).name);
end
uicontrol(f,'Style','text','String',sprintf('Up / Down - select track\nLeft / Right - selet frame\nHome / End - select first / last frame\nShift + Arrows - shift cut-out\n%s / %s - set begin / end of track (+-1 frame)\nTags: %s',keyBegin,keyEnd,tagTxt(2:end)),'Units','normalized','Position',[0.64 0.89 0.16 0.1],'Max',2);
% uicontrol(f,'Style','pushbutton','String','Reset track to CME detection','Units','normalized','Position',[0.78 0.96 .1 .03],'Max',2,'Callback',@resetTrack);
bSave = uicontrol(f,'Style','pushbutton','String','Save changes','Units','normalized','Position',[0.8 0.96 .08 .03],'Enable','off','Callback',@saveTrks);
uicontrol(f,'Style','pushbutton','String','Resave all cut-outs','Units','normalized','Position',[0.8 0.92 .08 .03],'Callback',@resaveAll);
% chCh = uicontrol(f,'Style','checkbox','String','Show all channels','Value',0,'Units','normalized','Position',[0.8 0.885 .08 .03],'Callback',@showChannels);
sld = uicontrol(f,'Style','slider','Min',1,'Max',N,'SliderStep',[1 1]./(N-1),'Value',frm(idx),'Units','normalized','Position',[.05 0 .95 .02],'Callback',@sldChange);
sldChange(sld);
if any(changed)
    bSave.Enable = 'on';
end

function graphButDown(~,e)
    sld.Value = max(min(round(e.IntersectionPoint(1)),sld.Max),sld.Min);
    sldChange(sld);
end
function cutButDown(s,~)
    if s==cut{2}
        ci = 3-ci;
        ci(ci==0) = 3;
    else
        ci = 4-ci;
    end
    showChannels();
end
function imgClick(obj,~)
    p = obj.Parent.CurrentPoint(1,1:2);
    eidx = find(arrayfun(@(x)any(x.f==frm(idx)),trk));
    eidx2 = arrayfun(@(x)find(x.f==frm(idx),1),trk(eidx));
    crd = [arrayfun(@(t,i)t.x(i)+t.cx(i),trk(eidx),eidx2);arrayfun(@(t,i)t.y(i)+t.cy(i),trk(eidx),eidx2)].';
    d = sqrt(sum((crd-p).^2,2));
    [d,di] = min(d);
%     if (d<20)
        lst.Value = eidx(di);
        frm(eidx(di)) = frm(idx);
        lstChange(lst);
%     end
end
function sldChange(obj,~)
    frm(idx) = round(get(obj,'Value'));
    set(txt,'String',sprintf('Frame: %i / %i',frm(idx),N));
    img.CData = I{1}(rCCP+1:end-rCCP,rCCP+1:end-rCCP,frm(idx));
    posChange;
    uicontrol(txt);
end
function posChange()
    fi = frm(idx)-trk(idx).f(1)+1;
    if fi<1 || fi>length(trk(idx).f)
        rec.LineStyle = ':';
        MD.MarkerIndices = [];
    else
        rec.LineStyle = '-';
        MD.MarkerIndices = 1;
    end
    fi = max(min(fi,length(trk(idx).f)),1);
    rec.Position = [trk(idx).x(fi)+trk(idx).cx(fi)-rCCP-0.5,trk(idx).y(fi)+trk(idx).cy(fi)-rCCP-0.5,rec.Position(3:4)];
    showChannels();
end
function lstChange(obj,~)
    idx = get(obj,'Value');
    sld.Value = frm(idx);
    sldChange(sld);
end
function keyPress(~,e)
    fi = frm(idx)-trk(idx).f(1)+1;
    if length(e.Modifier)==1 && strcmpi(e.Modifier{1},'shift')
        if fi>=1 && fi<=length(trk(idx).f)
            switch e.Key
                case 'leftarrow'
                    trk(idx).cx(fi) = trk(idx).cx(fi) - 1;
                    changed(idx) = true;
                    posChange;
                case 'rightarrow'
                    trk(idx).cx(fi) = trk(idx).cx(fi) + 1;
                    changed(idx) = true;
                    posChange;
                case 'uparrow'
                    trk(idx).cy(fi) = trk(idx).cy(fi) - 1;
                    changed(idx) = true;
                    posChange;
                case 'downarrow'
                    trk(idx).cy(fi) = trk(idx).cy(fi) + 1;
                    changed(idx) = true;
                    posChange;
            end
        end
    else
        switch e.Key
            case 'leftarrow'
                sld.Value = max(sld.Value-1,sld.Min);
                sldChange(sld);
            case 'rightarrow'
                sld.Value = min(sld.Value+1,sld.Max);
                sldChange(sld);
            case 'home'
                sld.Value = trk(idx).f(1);
                sldChange(sld);
            case 'end'
                sld.Value = trk(idx).f(end);
                sldChange(sld);
            case 'uparrow'
                lst.Value = max(idx-1,1);
                lstChange(lst);
            case 'downarrow'
                lst.Value = min(idx+1,length(ID));
                lstChange(lst);
            case keyBegin
                if frm(idx)==trk(idx).f(1)-1
                    trk(idx).f = [frm(idx),trk(idx).f];
                    trk(idx).x = [-1000,trk(idx).x];
                    trk(idx).cx = [trk(idx).x(2)+trk(idx).cx(1)+1000,trk(idx).cx];
                    trk(idx).y = [-1000,trk(idx).y];
                    trk(idx).cy = [trk(idx).y(2)+trk(idx).cy(1)+1000,trk(idx).cy];
                    trk(idx).tag = [0,trk(idx).tag];
                    A{idx} = [0,0;A{idx}];
                    lst.String{idx} = sprintf('%3i   (%3i - %3i)',ID(idx),trk(idx).f(1),trk(idx).f(end));
                    posChange();
                    changed(idx) = true;
                elseif frm(idx)==trk(idx).f(1)+1
                    trk(idx).f(1) = [];
                    trk(idx).x(1) = [];
                    trk(idx).cx(1) = [];
                    trk(idx).y(1) = [];
                    trk(idx).cy(1) = [];
                    trk(idx).tag(1) = [];
                    A{idx}(1,:) = [];
                    lst.String{idx} = sprintf('%3i   (%3i - %3i)',ID(idx),trk(idx).f(1),trk(idx).f(end));
                    posChange();
                    changed(idx) = true;
                end
            case keyEnd
                if frm(idx)==trk(idx).f(end)+1
                    trk(idx).f = [trk(idx).f,frm(idx)];
                    trk(idx).x = [trk(idx).x,-1000];
                    trk(idx).cx = [trk(idx).cx,trk(idx).x(end-1)+trk(idx).cx(end)+1000];
                    trk(idx).y = [trk(idx).y,-1000];
                    trk(idx).cy = [trk(idx).cy,trk(idx).y(end-1)+trk(idx).cy(end)+1000];
                    trk(idx).tag = [trk(idx).tag,0];
                    A{idx} = [A{idx};0,0];
                    lst.String{idx} = sprintf('%3i   (%3i - %3i)',ID(idx),trk(idx).f(1),trk(idx).f(end));
                    posChange();
                    changed(idx) = true;
                elseif frm(idx)==trk(idx).f(end)-1
                    trk(idx).f(end) = [];
                    trk(idx).x(end) = [];
                    trk(idx).cx(end) = [];
                    trk(idx).y(end) = [];
                    trk(idx).cy(end) = [];
                    trk(idx).tag(end) = [];
                    A{idx}(end,:) = [];
                    lst.String{idx} = sprintf('%3i   (%3i - %3i)',ID(idx),trk(idx).f(1),trk(idx).f(end));
                    posChange();
                    changed(idx) = true;
                end
        end
        if length(e.Key)==1 && fi>=1 && fi<=length(trk(idx).f)
            keyIdx = find([tags.key]==e.Key,1);
            if ~isempty(keyIdx)
                trk(idx).tag(fi) = bitxor(trk(idx).tag(fi),tags(keyIdx).tag);
                showChannels();
                changed(idx) = true;
            end
        end
    end
    if any(changed)
        bSave.Enable = 'on';
    end
end
% function resetTrack(~,~)
%     trk(idx).cx(:) = 0;
%     trk(idx).cy(:) = 0;
%     changed(idx) = true;
%     posChange;
% end
function resaveAll(~,~)
    changed(:) = true;
end
function showChannels(~,~)
    fi = frm(idx)-trk(idx).f(1)+1;
    cutI = zeros([2*rCCP+1,2*rCCP+1,nCh]);
    for iii = 1:nCh
        cutI(:,:,iii) = I{iii}(round(rec.Position(2)+rCCP+0.5)+(0:2*rCCP),round(rec.Position(1)+rCCP+0.5)+(0:2*rCCP),frm(idx));
        if fi>=1 && fi<=length(trk(idx).f)
            A{idx}(fi,iii) = max(sum(sum(cutI(:,:,iii).*C)));
        end
    end
    RI{1}.YData = squeeze(sum(sum(cutI(:,:,1).*C)));
    RI{2}.YData = squeeze(sum(sum(cutI(:,:,2).*C)));
    TA{3}.YData = [min(min(A{idx}(:,1:2))),max(max(A{idx}(:,1:2)))];
    TA{1}.XData = trk(idx).f;
    TA{1}.YData = A{idx}(:,1);
    TA{2}.XData = trk(idx).f;
    TA{2}.YData = A{idx}(:,2);
    TA{3}.XData = frm(idx)*[1,1];
    for ti = 1:length(tags)
        TA{ti+3}.XData = find(bitand(trk(idx).tag,tags(ti).tag))+trk(idx).f(1)-1;
        TA{ti+3}.YData = TA{1}.YData(TA{ti+3}.XData-trk(idx).f(1)+1);
    end
    ax(2).YLimMode = 'auto';
    ax(2).YLim(1) = 0;
%     ax(2).YLim = [0,max(A{idx}(:))];
    TA{3}.YData = ax(2).YLim;
    ax(1).YLim = ax(2).YLim;
    if fi>=1 && fi<=length(trk(idx).f)
        MD.XData = trk(idx).x(fi)-round(trk(idx).x(fi))-trk(idx).cx(fi)+rCCP+1;
        MD.YData = trk(idx).y(fi)-round(trk(idx).y(fi))-trk(idx).cy(fi)+rCCP+1;
    end
%     if chCh.Value>0
% %         ti = (cutI-reshape(mi(:,1),1,1,[]))./(reshape(mi(:,2)-mi(:,1),1,1,[]));
        ti = (cutI-ax(2).YLim(1))./diff(ax(2).YLim);
        cut{ci(1)}.CData = cat(3,ti(:,:,2),ti(:,:,1),zeros(2*rCCP+1));
%     else
        cut{ci(2)}.CData = cutI(:,:,1);
        cut{ci(3)}.CData = cutI(:,:,2);
%     end
end
function saveTrks(~,~)
    expTrk.tracks(ID) = trk;
    save([folder 'Tracking' filesep 'exportedTracks.mat'],'-struct','expTrk');
    for it = find(changed)
        for ch = 1:length(chnl)
            E = zeros([2*rCCP+1,2*rCCP+1,length(trk(it).f)],'single');
            for fr = 1:length(trk(it).f)
                E(:,:,fr) = I{ch}(round(trk(it).y(fr)+trk(it).cy(fr))+(0:2*rCCP),round(trk(it).x(fr)+trk(it).cx(fr))+(0:2*rCCP),trk(it).f(fr));
            end
            saveTiff(E,[folder 'Tracking' filesep 'Track_' num2str(ID(it)) '_' num2str(chnl(ch).wavelength) '.tif']);
        end
    end
    changed(:) = false;
    bSave.Enable = 'off';
end
function closeFig(~,~)
    if any(changed)
        res = questdlg('Save changes?','Finish centering');
        if strcmpi(res,'Yes')
            saveTrks();
        elseif ~strcmpi(res,'No')
            return;
        end
    end
    delete(gcf);
end

end
