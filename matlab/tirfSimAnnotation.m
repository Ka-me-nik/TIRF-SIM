% tirfSimAnnotation - GUI for TIRF-SIM tracks annotation
%
% Optional parameter:
%    Folder - cell data folder ['' = will ask for the folder]
%    Annotator - name of the annotator ['' = will ask for your name]
%
% Named parameters:
%    Radius - cut-out radius (size = 2R+1) [10]
%    Marker - displayed marker for track positions ['x']
%    ColorCutout - cut-out overlay color [[.5 .5 1]]
%    ColorTracks - tracks marker color for individual status states [[1 .2 .2;1 1 .2;.2 1 .2]]
function tirfSimAnnotation(varargin)

ip = inputParser;
ip.addOptional('Folder', '', @ischar);
ip.addOptional('Annotator', '', @ischar);
ip.addParameter('Radius', 10, @(x)validateattributes(x,{'numeric'},{'positive','scalar','integer'}));
ip.addParameter('MarkerCutout', 'x', @ischar);
ip.addParameter('ColorCutout', [.5 .5 1], @(x) ischar(x) || (isnumeric(x) && length(x)==3));
ip.addParameter('MarkerTracks', ['x';'x';'x'], @(x)(ischar(x) && size(x,2)==1));
ip.addParameter('MarkerSizeTracks', [6; 6; 6], @(x)(isnumeric(x) && size(x,2)==1));
ip.addParameter('ColorTracks', [1 .5 .9;1 1 0;.5 1 .5], @(x)(isnumeric(x) && size(x,2)==3));
ip.addParameter('ColorGrid', [.5 .5 1], @(x) ischar(x) || (isnumeric(x) && length(x)==3));
ip.parse(varargin{:});
folder = ip.Results.Folder;        % cell folder
annotator = ip.Results.Annotator;  % annotator name
rCCP = ip.Results.Radius;          % cut-out radius (size = 2*r+1)
marker = ip.Results.MarkerCutout;  % displayed marker in the cutout
cMark = ip.Results.ColorCutout;    % cut-out overlay color
mTracks = ip.Results.MarkerTracks; % tracks marker (for individual status states)
sTracks = ip.Results.MarkerSizeTracks; % tracks marker size (for individual status states)
cTracks = ip.Results.ColorTracks;  % tracks marker color (for individual status states)
cGrid = ip.Results.ColorGrid;      % grid color
if size(mTracks,1)==1
    mTracks = repmat(mTracks,3,1);
end
if size(sTracks,1)==1
    sTracks = repmat(sTracks,3,1);
end
if size(cTracks,1)==1
    cTracks = repmat(cTracks,3,1);
end

rMark = 1000/65/2; % 1/4 um in TIRF-SIM resolution
keyBegin = 'b';
keyEnd = 'e';
% tags
tags = struct('name',{'ring','cluster'},'tag',{1,2},'key',{'q','c'},'marker',{'o','x'});
initCheck = 'ag4'; % (a)ll,(s)nakes, (d)istances, (g)rid, (4)88, (5)60

caption = 'TIRF-SIM annotation';
if isempty(folder)
    folder = uigetdir(pwd, 'Select the ''cell'' folder:');
    if folder==0
        error('Cell folder not specified.');
    end
end
if ~strcmp(folder(end), filesep)
    folder = [folder filesep];
end

fann = [];
if exist('annotator.txt','file')==2
    f = fopen('annotator.txt','r');
    fann = fgetl(f);
    fclose(f);
end
if isempty(annotator)
    if isempty(fann)
        annotator = inputdlg('Enter annotator name:',caption);
        if isempty(annotator)||isempty(annotator{1})
            error('Annotator name not specified.');
        end
        annotator = annotator{1};
    else
        annotator = fann;
    end
end
if ~strcmp(annotator,fann)
    f = fopen('annotator.txt','w');
    fprintf(f,'%s\n',annotator);
    fclose(f);
end

fn = [folder 'tracks-' annotator '.mat'];
if exist(fn,'file')~=2
    d = dir([folder 'tracks-*.mat']);
    assert(~isempty(d),'No initial tracks found.');
    ini = listdlg('PromptString','Choose initialization:','ListString',cellfun(@(t)t(8:end-4),{d.name},'unif',0),'SelectionMode','single');
    assert(~isempty(ini),'Initialization must be specified, cannot continue.');
    copyfile([d(ini).folder filesep d(ini).name],fn);
end
expTrk = load(fn);
chnl = expTrk.channels;
nCh = length(chnl);
trk = expTrk.tracks;
changed = false;
idx = 1;
frm = [trk.start];

fprintf('Loading data ...');
I = cell(1,nCh);
mi = zeros(nCh,2);
cRange = zeros(nCh,2);
pSize = 0;
for i = 1:nCh
    if strcmpi(chnl(i).movie(end-2:end),'mrc')
        [m,s] = ReadMRC([folder chnl(i).movie]);
        pSize = s.rez*1000;
        I{i} = rot90(m).*reshape(chnl(i).bleaching,1,1,[]);
    else
        [m,s] = loadTiff([folder chnl(i).movie]);
        I{i} = m.*reshape(chnl(i).bleaching,1,1,[]);
        if pSize==0
            pSize = 1000/s.XResolution;
        end
    end
    mi(i,:) = [min(I{i}(:)),max(I{i}(:))];
    cRange(i,:) = [0,mean(max(max(I{i})))];
    I{i} = padarray(I{i},rCCP*[1,1]);
end
if pSize==0
    pSize = 30.65;
end
cc = cRange(1,2)/cRange(end,2);
N = size(I{1},3);
fprintf(' done.\n');

% prepare masks (circles) for radial function computation
C = zeros(2*rCCP+1,2*rCCP+1,2*rCCP);
[x,y] = meshgrid(-rCCP-.495:.01:rCCP+.495);
for r = .25:.5:rCCP+.25
    tc = (x.^2+y.^2<r^2);
    C(:,:,r*2+.5) = reshape(mean(reshape(reshape(mean(reshape(tc,100,[])),2*rCCP+1,[]).',100,[])),size(C(:,:,1)));
end
for r = 2*rCCP+1:-1:2
    C(:,:,r) = C(:,:,r)-C(:,:,r-1);
end
C = C./sum(sum(C));

fprintf('Calculating intensities ...');
A = cell(1,length(trk));
tc = zeros([2*rCCP+1,2*rCCP+1,nCh]);
for i = 1:length(trk)
    trk(i).x = round(trk(i).x);
    trk(i).y = round(trk(i).y);
    for j = 1:length(trk(i).x)
        for c = 1:nCh
            tc(:,:,c) = I{c}(trk(i).y(j)+(0:2*rCCP),trk(i).x(j)+(0:2*rCCP),trk(i).start+j-1);
            A{i}(j,c) = max(sum(sum(tc(:,:,c).*C)));
        end
    end
end
fprintf(' done.\n');

% fprintf('Calculating distances ...');
% D = NaN(length(trk),length(trk),N);
% for i = 1:length(trk)
%     p = [trk(i).x;trk(i).y];
%     p = [p(:,1)*ones(1,trk(i).start-1),p,p(:,end)*ones(1,N-trk(i).start+1-size(p,2))];
%     for j = 1:length(trk)
%         D(i,j,trk(j).start:trk(j).start+numel(trk(j).x)-1) = sqrt((trk(j).x-p(1,trk(j).start:trk(j).start+numel(trk(j).x)-1)).^2+(trk(j).y-p(2,trk(j).start:trk(j).start+numel(trk(j).x)-1)).^2);
%     end
% end
% fprintf(' done.\n');

f = figure('Name','TIRF-SIM CCP annotation','Units','normalized','Position',[.04 .04 .9 .9],'Menubar','none','Toolbar','none','WindowKeyPressFcn',@keyPress,'CloseRequestFcn',@closeFig);
if ~verLessThan('matlab','9.4')
    f.WindowState = 'maximized';
end
axi = subplot('Position',[0.01,0.03,0.63,0.96]);
img = image(axi,I{1}(rCCP+1:end-rCCP,rCCP+1:end-rCCP,frm(idx)),'ButtonDownFcn',@imgClick);
colormap(gray(256));
set(axi,'CLim',mi(1,:));
axis(axi,'equal','off');
hold(axi,'on');
zoom reset
axi.Toolbar.Visible = 'on';
fi = frm(idx)-trk(idx).start+1;
s = size(img.CData).'.*[1,2,3]/4+.5;
grid = plot([s(2,:),.5,.5,.5;s(2,:),(size(img.CData,2)+.5)*[1,1,1]],[.5,.5,.5,s(1,:);(size(img.CData,1)+.5)*[1,1,1],s(1,:)],'Color',cGrid,'LineStyle','--');
rec = rectangle(axi,'Position',[trk(idx).x(fi)-rCCP-0.5,trk(idx).y(fi)-rCCP-0.5,2*rCCP+1,2*rCCP+1],'EdgeColor','g');
for i = 1:3
    all{i} = plot(axi,1,1,'LineStyle','none','Marker',mTracks(i),'MarkerEdgeColor',cTracks(i,:),'MarkerSize',sTracks(i),'MarkerIndices',find([trk.status]==i-1),'ButtonDownFcn',@imgClick,'Visible','off');
end
% axd = subplot('Position',[0.01,0.04,0.63,0.19],'ButtonDownFcn',@graphButDown,'Visible','off','XLim',[1,N],'YLim',[0,20]);
axd = subplot('Position',[0.01,0.02,0.63,0.01],'ButtonDownFcn',@graphButDown,'Visible','off','XLim',[1,N],'YLim',[0,20]);
hold(axd,'on');
axd.Position = [0.01,0.04,0.63,0.19];
% plot([1,1],axd.YLim,'Color',[.5 .5 .5]);
% for i = length(trk):-1:1
%     plot(trk(i).start:trk(i).start+numel(trk(i).x)-1,zeros(1,numel(trk(i).x)),'ButtonDownFcn',{@selectTrack,i});
% end
subplot('Position',[0.65,0.03,0.23,0.4]);
cut{1} = image(zeros([2*rCCP+1,2*rCCP+1,3]),'ButtonDownFcn',@cutClick);
axis equal off
zoom reset
hold on
tfi = 0:pi/50:2*pi;
plot(rMark*cos(tfi)+rCCP+1,rMark*sin(tfi)+rCCP+1,':','Color',cMark,'ButtonDownFcn',@cutClick);
plot(rMark/2*cos(tfi)+rCCP+1,rMark/2*sin(tfi)+rCCP+1,':','Color',cMark,'ButtonDownFcn',@cutClick);
MD = plot(rCCP+1,rCCP+1,marker,'Color',cMark,'ButtonDownFcn',@cutClick);
% plot(rCCP+1,rCCP+1,'o','MarkerSize',100,'Color',[.5 .5 1]);
line([.5,rCCP+1;2*rCCP+1.5,rCCP+1],[rCCP+1,.5;rCCP+1,2*rCCP+1.5],'Color',cMark,'LineStyle',':','ButtonDownFcn',@cutClick);
subplot('Position',[0.89,0.03,0.1,0.2]);
cut{2} = image(zeros([2*rCCP+1,2*rCCP+1,3]));
axis equal off
zoom reset
chs{1} = uicontrol(f,'Style','checkbox','Value',any(initCheck=='4')||(nCh==1),'String',chnl(1).wavelength,'ForegroundColor','g','Units','normalized','Position',[0.89 0.02 0.1 .02],'HorizontalAlignment','left','Callback',@channelSelect);
if nCh==1
    chs{1}.Enable = 'inactive';
end
subplot('Position',[0.89,0.24,0.1,0.2]);
cut{3} = image(zeros([2*rCCP+1,2*rCCP+1,3]));
axis equal off
zoom reset
if nCh==1
    ch2w = 560;
else
    ch2w = chnl(2).wavelength;
end
chs{2} = uicontrol(f,'Style','checkbox','Value',any(initCheck=='5'),'String',ch2w,'ForegroundColor','r','Units','normalized','Position',[0.89 0.23 0.1 .02],'HorizontalAlignment','left','Callback',@channelSelect);
ax(1) = subplot('Position',[0.65,0.455,0.23,0.2],'XLim',[-1,1]*(rCCP+.5)*pSize);
hold(ax(1),'on');
RI = cell(1,2);
RI{1} = plot(ax(1),pSize*(-2*rCCP:2*rCCP)/2,zeros(1,4*rCCP+1),'g');
RI{2} = plot(ax(1),pSize*(-2*rCCP:2*rCCP)/2,zeros(1,4*rCCP+1),'r');
% RI{2} = plot(ax(1),pSize*(1:2*rCCP)/2,zeros(1,2*rCCP),'r');
ax(2) = subplot('Position',[0.65,0.68,0.23,0.2],'ButtonDownFcn',@graphButDown);
hold(ax(2),'on');
TA = cell(1,3+length(tags));
TA{1} = plot(ax(2),1:N,zeros(1,N),'g');
TA{2} = plot(ax(2),1:N,zeros(1,N),'r');
TA{3} = line(ax(2),[1 1],[0,1],'Color',[.5 .5 .5]);
for i = 1:length(tags)
    TA{i+3} = scatter(ax(2),[],[],50,'b',tags(i).marker);
end
lst = uicontrol(f,'Style','listbox','String',arrayfun(@(x)listTxt(x),1:length(trk),'Uniform',0),'Units','normalized','Position',[0.89 0.45 0.1 0.54],'Value',idx,'Callback',@lstChange);
txt = uicontrol(f,'Style','text','String',['Frame: ' num2str(frm(idx)) ' / ' num2str(N)],'Units','normalized','Position',[0 0 .05 .02]);
tagTxt = [];
for i = 1:length(tags)
    tagTxt = sprintf('%s, %s - %s',tagTxt,tags(i).key,tags(i).name);
end
uicontrol(f,'Style','text','FontSize',7,'String',sprintf('Up, Down / Left, Right - select track / frame\nHome / End - select first / last frame\nCtrl+click / Delete / x / m - create / delete / split / merge track\nShift + Arrows - shift cut-out\n%s / %s - set begin / end of track\nspace / n - change status / next unsolved track\nTags: %s\na / s / d / g / 4 / 5 - all / snakes / distances / grid / 488 / 560\nClick - select nearest track / shift cutout',keyBegin,keyEnd,tagTxt(2:end)),'Units','normalized','Position',[0.64 0.89 0.16 0.1],'Max',2);
bSave = uicontrol(f,'Style','pushbutton','String','Save changes','Units','normalized','Position',[0.8 0.96 .08 .03],'Enable','off','Callback',@saveTrks);
chAll = uicontrol(f,'Style','checkbox','String','Show all tracks','Value',any(initCheck=='a'),'Units','normalized','Position',[0.8 0.945 .08 .015],'Callback',@posChange);
chSnake = uicontrol(f,'Style','checkbox','String','Show snakes','Value',any(initCheck=='s'),'Units','normalized','Position',[0.8 0.93 .08 .015],'Callback',@posChange);
chDist = uicontrol(f,'Style','checkbox','String','Show tracks distance','Value',any(initCheck=='d'),'Units','normalized','Position',[0.8 0.915 .08 .015],'Callback',@posChange);
chGrid = uicontrol(f,'Style','checkbox','String','Show grid','Value',any(initCheck=='g'),'Units','normalized','Position',[0.8 0.90 .08 .015],'Callback',@gridChange);
gridChange;
% uicontrol(f,'Style','pushbutton','String','Resave all cut-outs','Units','normalized','Position',[0.8 0.92 .08 .03],'Callback',@resaveAll);
% chCh = uicontrol(f,'Style','checkbox','String','Show all channels','Value',0,'Units','normalized','Position',[0.8 0.885 .08 .03],'Callback',@showChannels);
sld = uicontrol(f,'Style','slider','Min',1,'Max',N,'SliderStep',[1 1]./(N-1),'Value',frm(idx),'Units','normalized','Position',[.05 0 .95 .02],'Callback',@sldChange);
sldChange(sld);
if any(changed)
    bSave.Enable = 'on';
end

function txt = listTxt(index)
%     col = '000000';
    sta = '-';
    switch trk(index).status
        case 0
%             col = 'BF0000';
            sta = '  ';
        case 2
%             col = '00BF00';
            sta = 'F';
    end
    if any(mod(floor(trk(index).tag/2),2))
        sta = ['c',sta];
    else
        sta = ['  ',sta];
    end
    if any(mod(trk(index).tag,2))
        sta = ['r',sta];
    else
        sta = ['  ',sta];
    end
%     txt = sprintf('<html><body style="background-color:#%s">%3i   (%3i - %3i)</body></html>',col,index,trk(index).start,trk(index).start+numel(trk(index).x)-1);
%     txt = sprintf('<html><font color=#%s><b>%3i   (%3i - %3i)</b></font></html>',col,index,trk(index).start,trk(index).start+numel(trk(index).x)-1);
    txt = sprintf('%s  %3i   (%3i - %3i)',sta,index,trk(index).start,trk(index).start+numel(trk(index).x)-1);
end
function gridChange(~,~)
    if chGrid.Value>0
        val = 'on';
    else
        val = 'off';
    end
    for ig = 1:length(grid)
        grid(ig).Visible = val;
    end
end
function graphButDown(~,e)
    sld.Value = max(min(round(e.IntersectionPoint(1)),sld.Max),sld.Min);
    sldChange(sld);
end
function channelSelect(s,~)
    if chs{1}.Value==0 && chs{2}.Value==0
        if s==chs{1}
            chs{2}.Value = 1;
        else
            chs{1}.Value = 1;
        end
    end
    if nCh==1
        chs{1}.Value = 1;
    end
    showChannels();
end
function cutClick(obj,e)
    p = round(obj.Parent.CurrentPoint(1,1:2))-rCCP-1;
    if e.Button==1
        fi = frm(idx)-trk(idx).start+1;
        trk(idx).x(fi) = trk(idx).x(fi) + p(1);
        trk(idx).y(fi) = trk(idx).y(fi) + p(2);
        trk(idx).status = 1;
        changed = true;
        posChange;
        calcCirc(idx,frm(idx));
        lst.String{idx} = listTxt(idx);
    end
end
function imgClick(obj,e)
    p = round(obj.Parent.CurrentPoint(1,1:2));
    if e.Button==1
        if ismember('control',get(f,'currentModifier'))            
            idx = length(trk)+1;
            frm(idx) = round(sld.Value);
            trk(idx).start = frm(idx);
            trk(idx).x = p(1);
            trk(idx).y = p(2);
            trk(idx).status = 1;
            trk(idx).tag = 0;
            calcCirc(idx,frm(idx));
            lst.String{idx} = listTxt(idx);
            changed = true;
            lst.Value = idx;
            lstChange(lst);
        else
            eidx = find(arrayfun(@(x)any(x.start<=frm(idx)&frm(idx)<x.start+numel(x.x)),trk));
            eidx2 = arrayfun(@(x)frm(idx)-x.start+1,trk(eidx));
            crd = [arrayfun(@(t,i)t.x(i),trk(eidx),eidx2);arrayfun(@(t,i)t.y(i),trk(eidx),eidx2)].';
            d = sqrt(sum((crd-p).^2,2));
            [~,di] = min(d);
%             if (d<20)
                lst.Value = eidx(di);
                frm(eidx(di)) = frm(idx);
                lstChange(lst);
%             end
        end
    end
end
function sldChange(obj,~)
    frm(idx) = round(get(obj,'Value'));
    set(txt,'String',sprintf('Frame: %i / %i',frm(idx),N));
    img.CData = I{1}(rCCP+1:end-rCCP,rCCP+1:end-rCCP,frm(idx));
    posChange;
    uicontrol(txt);
end
function selectTrack(~,~,id)
    lst.Value = id;
    lstChange(lst);
end
function posChange(~,~)
    fi = frm(idx)-trk(idx).start+1;
    if fi<1 || fi>length(trk(idx).x)
        rec.LineStyle = ':';
        MD.MarkerIndices = [];
    else
        rec.LineStyle = '-';
        MD.MarkerIndices = 1;
    end
    fi = max(min(fi,length(trk(idx).x)),1);
    rec.Position = [trk(idx).x(fi)-rCCP-0.5,trk(idx).y(fi)-rCCP-0.5,rec.Position(3:4)];
    if chAll.Value>0
        ti = arrayfun(@(t)t.start<=frm(idx)&&frm(idx)<t.start+numel(t.x),trk);
        for ia = 1:3
            all{ia}.XData = arrayfun(@(t)t.x(frm(idx)-t.start+1),trk(ti));
            all{ia}.YData = arrayfun(@(t)t.y(frm(idx)-t.start+1),trk(ti));
            all{ia}.MarkerIndices = find([trk(ti).status]==ia-1);
            all{ia}.Visible = 'on';
        end
    else
        for ia = 1:3
            all{ia}.Visible = 'off';
        end
    end
    delete(axi.Children(1:end-5-6));
    nco = length(axi.ColorOrder);
    if chSnake.Value>0
        for k=1:length(trk)
            if chSnake.Value>0 && trk(k).start<=frm(idx) && trk(k).start+numel(trk(k).x)>frm(idx)
                plot(axi,trk(k).x,trk(k).y,'ButtonDownFcn',{@selectTrack,k},'Color',axi.ColorOrder(mod(k,nco)+1,:));
            end
        end
    end
    cla(axd);
    if chDist.Value>0
        p = [trk(idx).x;trk(idx).y];
        p = [p(:,1)*ones(1,trk(idx).start-1),p,p(:,end)*ones(1,N-trk(idx).start+1-size(p,2))];
        d = arrayfun(@(t)[t.start:t.start+numel(t.x)-1;sqrt((t.x-p(1,t.start:t.start+numel(t.x)-1)).^2+(t.y-p(2,t.start:t.start+numel(t.x)-1)).^2)],trk,'unif',0);
        for k=1:length(d)
%             axd.Children(k).YData = d{k}(2,:);
            plot(axd,d{k}(1,:),d{k}(2,:),'ButtonDownFcn',{@selectTrack,k},'Color',axi.ColorOrder(mod(k,nco)+1,:));
%             plot(axd,1:N,D(idx,k,:),'ButtonDownFcn',{@selectTrack,k});
        end
%         axd.Children(length(trk)+1).XData = [1,1]*frm(idx);
        plot(axd,frm(idx)*[1,1],axd.YLim,'Color',[.5 .5 .5]);
        axi.Position = [0.01,0.23,0.63,0.76];
        axd.Visible = 'on';
    else
        axd.Visible = 'off';
        axi.Position = [0.01,0.03,0.63,0.96];
    end
    showChannels();
end
function lstChange(obj,~)
    idx = get(obj,'Value');
    sld.Value = frm(idx);
    sldChange(sld);
end
function keyPress(~,e)
    fi = frm(idx)-trk(idx).start+1;
    if length(e.Modifier)==1 && strcmpi(e.Modifier{1},'shift') && contains(e.Key,'arrow')
        if fi==0
            trk(idx).start = trk(idx).start-1;
            trk(idx).x = [trk(idx).x(1),trk(idx).x];
            trk(idx).y = [trk(idx).y(1),trk(idx).y];
            trk(idx).tag = [0,trk(idx).tag];
            A{idx} = [0,0;A{idx}];
            fi = 1;
        elseif fi==length(trk(idx).x)+1
            trk(idx).x = [trk(idx).x,trk(idx).x(end)];
            trk(idx).y = [trk(idx).y,trk(idx).y(end)];
            trk(idx).tag = [trk(idx).tag,0];
            A{idx} = [A{idx};0,0];
        end
        if fi>=1 && fi<=length(trk(idx).x)
            switch e.Key
                case 'leftarrow'
                    trk(idx).x(fi) = trk(idx).x(fi) - 1;
                case 'rightarrow'
                    trk(idx).x(fi) = trk(idx).x(fi) + 1;
                case 'uparrow'
                    trk(idx).y(fi) = trk(idx).y(fi) - 1;
                case 'downarrow'
                    trk(idx).y(fi) = trk(idx).y(fi) + 1;
            end
            trk(idx).status = 1;
            changed = true;
            posChange;
            calcCirc(idx,frm(idx));
            lst.String{idx} = listTxt(idx);
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
                sld.Value = trk(idx).start;
                sldChange(sld);
            case 'end'
                sld.Value = trk(idx).start+numel(trk(idx).x)-1;
                sldChange(sld);
            case 'uparrow'
                lst.Value = max(idx-1,1);
                lstChange(lst);
            case 'downarrow'
                lst.Value = min(idx+1,length(trk));
                lstChange(lst);
            case 'n'
                ti = find([trk.status]<2);
                if any(ti>idx)
                    lst.Value = ti(find(ti>idx,1));
                    lstChange(lst);
                elseif ~isempty(ti)
                    lst.Value = ti(1);
                    lstChange(lst);
                else
                    questdlg('Wow, you''ve done it!','Finished','OK FANTASTIC','OK FANTASTIC');
                end
            case 'space'
                if trk(idx).status<2
                    trk(idx).status = 2;
                else
                    trk(idx).status = 1;
                end
                changed = true;
                lst.String{idx} = listTxt(idx);
                posChange;
            case 'a'
                chAll.Value = 1-chAll.Value;
                posChange;
            case 'd'
                chDist.Value = 1-chDist.Value;
                posChange;
            case 's'
                chSnake.Value = 1-chSnake.Value;
                posChange;
            case 'g'
                chGrid.Value = 1-chGrid.Value;
                gridChange;
            case '4'
                chs{1}.Value = 1-chs{1}.Value;
                channelSelect(chs{1});
            case '5'
                chs{2}.Value = 1-chs{2}.Value;
                channelSelect(chs{2});
            case keyBegin
                if frm(idx)<trk(idx).start
                    d = trk(idx).start-frm(idx);
                    if d>1 && ~strcmp(questdlg(sprintf('Stretch start of track by %i frames?',d),'Track start','OK','Cancel','OK'),'OK')
                        return;
                    end
                    trk(idx).start = trk(idx).start-d;
                    trk(idx).x = [trk(idx).x(1)*ones(1,d),trk(idx).x];
                    trk(idx).y = [trk(idx).y(1)*ones(1,d),trk(idx).y];
                    trk(idx).tag = [zeros(1,d),trk(idx).tag];
                    A{idx} = [zeros(d,2);A{idx}];
                    for ii=frm(idx):frm(idx)+d-1
                        calcCirc(idx,ii);
                    end
                    trk(idx).status = 1;
                    changed = true;
                    posChange();
                    lst.String{idx} = listTxt(idx);
                elseif frm(idx)>trk(idx).start && frm(idx)<=trk(idx).start+numel(trk(idx).x)-1
                    d = frm(idx)-trk(idx).start;
                    if d>1 && ~strcmp(questdlg(sprintf('Shrink start of track by %i frames?',d),'Track start','OK','Cancel','OK'),'OK')
                        return;
                    end
                    trk(idx).start = trk(idx).start+d;
                    trk(idx).x(1:d) = [];
                    trk(idx).y(1:d) = [];
                    trk(idx).tag(1:d) = [];
                    A{idx}(1:d,:) = [];
                    trk(idx).status = 1;
                    changed = true;
                    posChange();
                    lst.String{idx} = listTxt(idx);
                end
            case keyEnd
                if frm(idx)>=trk(idx).start+numel(trk(idx).x)
                    d = frm(idx)-trk(idx).start-numel(trk(idx).x)+1;
                    if d>1 && ~strcmp(questdlg(sprintf('Stretch end of track by %i frames?',d),'Track end','OK','Cancel','OK'),'OK')
                        return;
                    end
                    trk(idx).x = [trk(idx).x,trk(idx).x(end)*ones(1,d)];
                    trk(idx).y = [trk(idx).y,trk(idx).y(end)*ones(1,d)];
                    trk(idx).tag = [trk(idx).tag,zeros(1,d)];
                    A{idx} = [A{idx};zeros(d,2)];
                    for ii=frm(idx)-d+1:frm(idx)
                        calcCirc(idx,ii);
                    end
                    trk(idx).status = 1;
                    changed = true;
                    posChange();
                    lst.String{idx} = listTxt(idx);
                elseif frm(idx)<trk(idx).start+numel(trk(idx).x)-1 && frm(idx)>=trk(idx).start
                    d = trk(idx).start+numel(trk(idx).x)-1-frm(idx);
                    if d>1 && ~strcmp(questdlg(sprintf('Shrink end of track by %i frames?',d),'Track end','OK','Cancel','OK'),'OK')
                        return;
                    end
                    trk(idx).x(end-d+1:end) = [];
                    trk(idx).y(end-d+1:end) = [];
                    trk(idx).tag(end-d+1:end) = [];
                    A{idx}(end-d+1:end,:) = [];
                    trk(idx).status = 1;
                    changed = true;
                    posChange();
                    lst.String{idx} = listTxt(idx);
                end
            case 'x'
                if fi<2 || fi>numel(trk(idx).x)-2 || ~strcmp(questdlg('Really split this track?',sprintf('Track %i',idx),'Yes','No','Yes'),'Yes')
                    return;
                end
                idx2 = length(trk)+1;
                frm(idx2) = round(sld.Value)+1;
                trk(idx2).start = frm(idx2);
                trk(idx2).x = trk(idx).x(fi+1:end);
                trk(idx2).y = trk(idx).y(fi+1:end);
                trk(idx2).status = 1;
                trk(idx2).tag = trk(idx).tag(fi+1:end);
                A{idx2} = A{idx}(fi+1:end,:);
                lst.String{idx2} = listTxt(idx2);
                trk(idx).x(fi:end) = [];
                trk(idx).y(fi:end) = [];
                trk(idx).status = 1;
                trk(idx).tag(fi:end) = [];
                A{idx}(fi:end,:) = [];
                lst.String{idx} = listTxt(idx);
                changed = true;
                lst.Value = idx2;
                lstChange(lst);
            case 'm'
                crd = round([trk(idx).x(fi),trk(idx).y(fi)]);
                ti = find(arrayfun(@(t)t.start<=frm(idx)&&frm(idx)<t.start+numel(t.x),trk));
                idx2 = ti(arrayfun(@(t)round(t.x(frm(idx)-t.start+1))==crd(1)&&round(t.y(frm(idx)-t.start+1))==crd(2),trk(ti)));
                idx3 = find(idx2==idx);
                if ~isempty(idx3) && length(idx2)==2
                    idx2(idx3) = [];
                    idx3 = [idx;idx2];
                    if trk(idx2).start+length(trk(idx2).x)-trk(idx).start < trk(idx).start+length(trk(idx).x)-trk(idx2).start
                        idx3 = [idx2;idx];
                    end
                    if strcmp(questdlg('Really merge these tracks?',sprintf('Track %i + %i',idx3(1),idx3(2)),'Yes','No','Yes'),'Yes')
                        frm(idx3(1)) = frm(idx);
                        trk(idx3(1)).x = [trk(idx3(1)).x(1:frm(idx)-trk(idx3(1)).start+1),trk(idx3(2)).x(frm(idx)-trk(idx3(2)).start+2:end)];
                        trk(idx3(1)).y = [trk(idx3(1)).y(1:frm(idx)-trk(idx3(1)).start+1),trk(idx3(2)).y(frm(idx)-trk(idx3(2)).start+2:end)];
                        trk(idx3(1)).tag = [trk(idx3(1)).tag(1:frm(idx)-trk(idx3(1)).start+1),trk(idx3(2)).tag(frm(idx)-trk(idx3(2)).start+2:end)];
                        trk(idx3(1)).status = 1;
                        A{idx3(1)} = [A{idx3(1)}(1:frm(idx)-trk(idx3(1)).start+1,:);A{idx3(2)}(frm(idx)-trk(idx3(2)).start+2:end,:)];
                        trk(idx3(2)) = [];
                        frm(idx3(2)) = [];
                        A(idx3(2)) = [];
                        lst.String = arrayfun(@(x)listTxt(x),1:length(trk),'Uniform',0);
                        changed = true;
                        lst.Value = idx3(1);
                        if idx3(2)<idx3(1)
                            lst.Value = idx3(1)-1;
                        end
                        lstChange(lst);
                    end
                end
            case 'delete'
                if length(trk)<2 || ~strcmp(questdlg('Really delete this track?',sprintf('Track %i',idx),'Yes','No','Yes'),'Yes')
                    return;
                end
                trk(idx) = [];
                frm(idx) = [];
                A(idx) = [];
                lst.String = arrayfun(@(x)listTxt(x),1:length(trk),'Uniform',0);
                changed = true;
                lst.Value = min(idx,length(trk));
                lstChange(lst);
        end
        if length(e.Key)==1 && fi>=1 && fi<=length(trk(idx).x)
            keyIdx = find([tags.key]==e.Key,1);
            if ~isempty(keyIdx)
                trk(idx).tag(fi) = bitxor(trk(idx).tag(fi),tags(keyIdx).tag);
                showChannels();
                trk(idx).status = 1;
                changed = true;
                lst.String{idx} = listTxt(idx);
                posChange;
            end
        end
    end
    if changed
        bSave.Enable = 'on';
    end
end
% function resetTrack(~,~)
%     trk(idx).cx(:) = 0;
%     trk(idx).cy(:) = 0;
%     changed(idx) = true;
%     posChange;
% end
% function resaveAll(~,~)
%     changed(:) = true;
% end
function showChannels(~,~)
    fi = frm(idx)-trk(idx).start+1;
    cutI = zeros([2*rCCP+1,2*rCCP+1,nCh]);
    for iii = 1:nCh
        cutI(:,:,iii) = I{iii}(round(rec.Position(2)+rCCP+0.5)+(0:2*rCCP),round(rec.Position(1)+rCCP+0.5)+(0:2*rCCP),frm(idx));
%         if fi>=1 && fi<=length(trk(idx).x)
%             A{idx}(fi,iii) = max(sum(sum(cutI(:,:,iii).*C)));
%         end
        rad = squeeze(sum(sum(cutI(:,:,iii).*C)));
        RI{iii}.YData = [rad(end:-1:2);rad];
        TA{iii}.XData = trk(idx).start:trk(idx).start+numel(trk(idx).x)-1;
        TA{iii}.YData = A{idx}(:,iii);
    end
    TA{3}.YData = [min(min(A{idx}(:,1:nCh))),max(max(A{idx}(:,1:nCh)))];
    TA{3}.XData = frm(idx)*[1,1];
    for ti = 1:length(tags)
        TA{ti+3}.XData = find(bitand(trk(idx).tag,tags(ti).tag))+trk(idx).start-1;
        TA{ti+3}.YData = TA{1}.YData(TA{ti+3}.XData-trk(idx).start+1);
    end
    ax(2).YLimMode = 'auto';
    ax(2).YLim(1) = 0;
%     ax(2).YLim = [0,max(A{idx}(:))];
    TA{3}.YData = ax(2).YLim;
    ax(1).YLim = ax(2).YLim;
    if fi>=1 && fi<=length(trk(idx).x)
        MD.XData = trk(idx).x(fi)-round(trk(idx).x(fi))+rCCP+1;
        MD.YData = trk(idx).y(fi)-round(trk(idx).y(fi))+rCCP+1;
    end
    
    cut{2}.CData = cutI(:,:,1);
    if nCh>1
        cutI(:,:,2) = cutI(:,:,2)*cc;
        cut{3}.CData = cutI(:,:,2);
    end
    ti = (cutI-ax(2).YLim(1))./diff(ax(2).YLim);
    if nCh==1
        ti(:,:,2) = 0;
    end
    if chs{1}.Value==1 && chs{2}.Value==1
        cut{1}.CData = cat(3,ti(:,:,2),ti(:,:,1),zeros(2*rCCP+1));
    elseif chs{1}.Value==1
        cut{1}.CData = cutI(:,:,1);
    else
        cut{1}.CData = cutI(:,:,2);
    end
end
function saveTrks(~,~)
    expTrk.tracks = trk;
    save(fn,'-struct','expTrk');
    changed = false;
    bSave.Enable = 'off';
end
function closeFig(~,~)
    if changed
        res = questdlg('Save changes?','Finish centering');
        if strcmpi(res,'Yes')
            saveTrks();
        elseif ~strcmpi(res,'No')
            return;
        end
    end
    delete(f);
end
function calcCirc(iTrk,iFrm)
    if iFrm>=trk(iTrk).start && iFrm<trk(iTrk).start+numel(trk(iTrk).x)
        fi = iFrm-trk(iTrk).start+1;
        for ci = 1:nCh
            A{iTrk}(fi,ci) = max(sum(sum(I{ci}(round(trk(iTrk).y(fi))+(0:2*rCCP),round(trk(iTrk).x(fi))+(0:2*rCCP),iFrm).*C)));
        end
    end
end
function calcDist(iTrk,iFrm)
end

end
