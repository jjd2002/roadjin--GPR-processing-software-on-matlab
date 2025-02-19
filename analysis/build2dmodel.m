function build2dmodel()
%
% BUILD2DMODEL : GUI utility to construct a 2D model for processing
%                with GET2DVELOVITY/SPLITSTEPMIG/MIGPSPI and
%                SPLITSTEP2DMODEL  
%
%     Usage    : build2dmodel;
%
%     REQUIRES : setpolygon.m, getpoint.m, snappolytoaxes.m, centroid.m,
%                getprops.m, matgprwindows.m 
%
% ==> The model comprices an ensemble of objects with polygonal or circular
%     cross-sections. There is no limit to the number of vertices. The
%     order in which objects are introduced in the structural the model
%     determines whether they will appear in the velocity model, or not.
%     The program implements a FRONT to BACK object hierarchy:
%     ==> Foreground objects will appear whole in the model. 
%     ==> Background objects that are totally obscured will not appear in
%         the velocity model.
%     ==> Background objects that are partially obscured will partially
%         appear in the velocity model.
%     ==> Therefore, the objects appearing LATER in the list may obscure
%         those appearing earlier in the list! 
%
% ==> The object coordinates and properties are introduced in three ways: 
%     1. Graphically, by poining and clicking on the model figure, 
%     2. Through dialog boxes, and,
%     3. From ASCII files on disk. These are very simple in structure,
%        comprising only a two column list with the horizontal (x) and
%        vertical (z) coordinates of the object.
% ==> BUILD2DMODEL provides for the graphical editing of object shape and 
%     properties (insert, move or remove vertices and review EM properties, 
%     with unlimited levels of undo). 
%
% ==> The program comprises a parent GUI function with menus, the items of 
%     which call function handle functions. Information between the main 
%     and handle functions is communicated with global variables, but is
%     otherwise kept strictly within the scope of BUILD2DMODEL.  These
%     variables are:
%     nbody    : The number of bodies in the model
%     inpmod   : Input mode (1 for typing, 2 for graphically entering 
%                information)
%     bodies   : (NBODY x 4)  cell array holding body data. 
%                bodies{j,1}  is a vector holding the X-vertices of the 
%                             j'th body 
%                bodies{j,2}  is a vector holding the Z-vertices of the 
%                             j'th body 
%                bodies{j,3}  is a 4-vector holding the EM properties of 
%                             the j'th body 
%                bodies{j,4}  is a 5-vector holding the ASCII code of the 
%                             j'th body's tag
%     handles  : (NBODY x 4)  cell array holding the handles of the gpaphic 
%                             objects
%                handles{j,1} is a vector holding the handles of the vertex 
%                             of the j'th body
%                handles{j,2} is a vector holding the patch of the j'th
%                             body
%                handles{j,3} is a vector holding the graphic data for the 
%                             "ID tag" of the j'th body 
%                handles{j,4} is a vector holding the gpahic objects for the 
%                             "EM properties" of the j'th body
%     backup   : (M x NBODY x 4) cell array holding backup body data for undo 
%                 actions. M is arbitrary and self updating
%
% ==> The free format ASCII model file is structured as follows:
%         antenna                float, antenna frequency.
%         Xmin Xmax Zmin Zmax    float, horizontal and vertical extent of
%                                the model.
%         nbody                  int, number of bodies in the model.
%         tag                    string, Tag of 1st body.
%         nv, K, Q, mu, rho      int + 4 x float. Respectively are:      
%                                a) The # vertices in 1st body, 
%                                b) its relative dielectric constant,
%                                c) its quality factor, 
%                                d) its relative magnetic permeability,
%                                e) its resistivity 
%         x(1)       z(1)        float, float: coordinates of the first
%                                vertex 
%         x(2)       z(2)        float, float: coordinates of the second
%                                vertex 
%         ...        ...
%         x(nv)      z(nv)       float, float: coordinates of the last 
%                                vertex 
%         tag                    string, Tag of the 2nd body.
%         nv, K, Q, mu, rho      int + 4 x float: Respectively # of
%                                vertices, relative dielectric constant,
%                                quality factor, relative permeability and
%                                resistivity of the 2nd body
%         x(1)       z(1)        float, float: coordinates of the first
%                                vertex of the 2nd body 
%         .... and so on 
%
%  Author : Andreas Tzanis
%           Department of Geophysics, 
%           University of Athens
%           atzanis@geol.uoa.gr
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%

global ENVAR
global bodies handles nbody inpmod backup model_was_saved
model_was_saved = 0;
%%%%%   initialize 
inpmod = 2;
nbody = 0;
bodies  = cell(1,4);                     % initialize cell arrays
backup  = {};
handles = cell(1,4);
%%%%%    Create model axes
modf   = figure('Name','Model Builder', ...
                'tag','modelfigure', ... 
                'Position', [150 150 800 500], ...
                'Numbertitle','off', ...
                'menubar','none',  ...
                'DeleteFcn',['matgprwindows(''updateremove'', ' ...
                             'findobj(''tag'',''modelfigure''));']);
figuretools;
%%% Set up UI MENUS
%%% -----------------------------------------------------------------------
%%% I/O OPTIONS
iomenu = uimenu('Label',' Data');
uimenu(iomenu,'Label','New Model',   'Callback', @newmodel);
uimenu(iomenu,'Label','Import Model','Callback', @readmodel);
uimenu(iomenu,'Label','Save Model',  'Callback', @savemodel);
uimenu(iomenu,'Label','Quit',        'Callback', @quitprg);
%%% -----------------------------------------------------------------------
%%% PROCESSING OPTIONS
crmenu = uimenu('Label',' Actions');
% input mode
inpm = uimenu(crmenu,'Label','Input Mode');
  uimenu(inpm,'Label','Use Fingers','Callback', {@inputmode,'fingers'});
  uimenu(inpm,'Label','Use Pointer','Callback', {@inputmode,'usepointer'});
  uimenu(inpm,'Label','From File  ','Callback', {@inputmode,'readfile'});
% Insert structural elements
uimenu(crmenu,'Label','Add Polygon', 'Callback', @addpolygon);
uimenu(crmenu,'Label','Add Circle',  'Callback', @addcircle);
uimenu(crmenu,'Label','Undo Shape',  'Callback', @delbody1);
% Modify bodies and object properties
uimenu(crmenu,'Label','Review EM Properties', 'Separator', 'on', ...
                                     'Callback', @reviewproperties1);
%RELOCATE VERTEX     
uimenu(crmenu,'Label','Move Vertex', 'Callback', @relocatevertex);
%ADD VERTEX     
uimenu(crmenu,'Label','Insert Vertex','Callback',@addvertex);
%REMOVE VERTEX     
uimenu(crmenu,'Label','Remove Vertex','Callback',@removevertex);
% UNDO     
uimenu(crmenu,'Label','Undo Changes','Callback', @undochanges);
% Presentation utilities 
uimenu(crmenu,'Label','Peek with cursor','Separator','on', ...
              'Callback', 'getpoint('' meters'', '' meters'');');
uimenu(crmenu,'Label','Recolor','Separator','on', ...
                                     'Callback', @recolor);
mrkm = uimenu(crmenu,'Label','Markers','Separator','on');       
       uimenu(mrkm,'Label','Hide',   'Callback', @hidemarkers); 
       uimenu(mrkm,'Label','Show',   'Callback', @showmarkers);
       uimenu(mrkm,'Label','Solid',  'Callback', @solidmarkers);
       uimenu(mrkm,'Label','Hollow', 'Callback', @hollowmarkers);
tagm = uimenu(crmenu,'Label','Tags');       
       uimenu(tagm,'Label','Hide',   'Callback', @hidetags); 
       uimenu(tagm,'Label','Show',   'Callback', @showtags);
propm = uimenu(crmenu,'Label','EM Properties');       
       uimenu(propm,'Label','Show',  'Callback', @showemprops);
       uimenu(propm,'Label','Hide',  'Callback', @hideemprops);
% -------------------------------------------------------------------------
% Setup the MATGPR fast window switching service
matgprwindows('setup');
matgprwindows('updateadd');
% -------------------------------------------------------------------------
% ANTENNA FREQUENCY  
uicontrol('BackGroundColor','y','Style','edit','tag', 'antennaf', ...
    'Units','normalized','Position',[.70 .93 .2 .05],  ...
    'String', 'Set antenna frequency in MHz');
return
%%% END FUNCTION build2dmodel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newmodel(vhandle,eventdata) 
%%%%%   Start New model
global bodies handles nbody inpmod backup model_was_saved
inpmod = 2;                           % default input mode is "POINTER"
bodies  = cell(1,4);                  % initialize cell arrays
backup  = {};
handles = cell(1,4);
mark1   = [];
modax = findobj('tag','modelaxes');   % Ensure model axes is current object
if ~isempty(modax),
    cla
end
%%%%%   Define the horizontal and vertical extent of the model (background)
z0 = 0;    zn = 10;   
x0 = 0;    xn = 20;                   % default starting values
clb = cell(3,1);                                          
clb(1) = cellstr(['LEFT  limit']);
clb(2) = cellstr(['RIGHT limit']);
clb(3) = cellstr(['Vertical (depth) extent']);
cldef(1) = cellstr(num2str(x0));
cldef(2) = cellstr(num2str(xn));
cldef(3) = cellstr(num2str(zn));
answer = inputdlg(clb,'Horizontal and Vertical extent',1,cldef);
if isempty(answer),                    % Canceled
    return; 
end;
lb = char(answer);                                                            
for i=1:3
    comma = findstr(lb(i,:),',');
    if ~isempty(comma), 
        for j=1:length(comma); 
            lb(i,comma) = '.'; 
        end;
        clear comma
    end
end
x0 = str2num(lb(1,:));
xn = str2num(lb(2,:));
zn = str2num(lb(3,:));
%%%%%   Draw axes
plot([x0 xn xn x0], [z0 z0 zn zn ], 'k', 'tag', 'dummy'); 
axis([x0,xn,z0,zn]); axis ij;
set(gca,'pos',[0.1 0.12 0.8 0.8],'tag','modelaxes');
set(gca,'nextplot','add');
ylabel('Depth (m)'); 
xlabel('Distance (m)')
%-------------------------------------------------------------------------
% Set up the default background uniform structure - This is the FIRST object 
nbody = 1;
nmark = 4;
x = [x0 xn xn x0]';
z = [z0 z0 zn zn ]';
for i = 1:nmark;
%%% Plot the outline of the polygon            
    mark1(i) = plot(x(i),z(i),'sb','markerfacecolor','w'); %#ok<AGROW>
    set(mark1(i),'MarkerSize',4,'LineWidth',1.0)
end
%%% Draw the background object
patch3 = patch(x,z,[rand(1,1), rand(1,1) rand(1,1)],'edgealpha',0);
%%% Antenna central frequency for new model
clb = cellstr('Please set the Antenna Frequency. Default = 250 MHz');
answer = inputdlg(clb,'New Model: Antenna Frequency',1, {'250'} );     
if isempty(answer),                        % Abandon creation of new model!  
    delete(mark1(1:nmark));   delete(patch3);
    bodies(nbody,:) = [];     backup = {};
    nbody = nbody - 1;
end;
antenna = str2num(answer{1});
set(findobj('tag','antennaf'),'String',answer{1})   
%%% Background EM properties 
bodyprops = getprops(antenna);
if isempty(bodyprops),          % Last minute decision to abandon new model
    delete(mark1(1:nmark));   delete(patch3);
    bodies(nbody,:) = [];     backup = {};
    nbody = nbody - 1;
    return
end
%%% Fill cell-arrays with background object
bodies{nbody,1}  = x(1:nmark);
bodies{nbody,2}  = z(1:nmark);
bodies{nbody,3}  = bodyprops;
handles{nbody,1} = mark1(1:nmark);
handles{nbody,2} = patch3;
deftag = cellstr('Background');
bodies{nbody,4} = double('Background');
[tx, tz] = centroid(x(1:nmark),z(1:nmark));
txt = text(tx,tz,'Background','horizontalal','center', ... 
     'fontweight','bold','fontsize',9,'erase','xor');
handles{nbody,3} = txt;
% draw and hide object properties 
str = ['r=' num2str(bodyprops(4)) ', K=' num2str(bodyprops(1)) ...
    ', m=' num2str(bodyprops(3))];
txt = text(tx,tz,str,'horizontalal','center', ...
     'fontname','symbol','fontweight','bold','fontsize',10,...
     'erase','xor','visible','off');
handles{nbody,4} = txt;
%--------------------------------------------------------------------------
model_was_saved = 0;
return
%%% END HANDLE FUNCTION newmodel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function readmodel(vhandle,eventdata) 
%%%%%   Import 2d model from disk file
global ENVAR
global bodies handles nbody inpmod backup model_was_saved
bodies  = cell(1,4);                  % initialize cell arrays
handles = cell(1,4);
modax = findobj('tag','modelaxes');   % Ensure model axes is current object
if ~isempty(modax),
    cla
end
[mfname, mpname]= uigetfile('*.dat','Give model file name', ...
    ENVAR.currentworkdir);
if mfname == 0, 
    return; 
end
fid = fopen([mpname mfname],'r');
antenna = fscanf(fid,' %f \n',1);        % get antenna frequency
set(findobj('tag','antennaf'),'String',num2str(antenna))   
dummy = fscanf(fid,' %f %f %f %f \n',4); % read axis bounds
xl = [dummy(1) dummy(2)];
zl = [dummy(3) dummy(4)];
%%%%%   Create model axes 
plot([xl(1) xl(2) xl(2) xl(1)], [zl(1) zl(1) zl(2) zl(2) ], ...
    'k', 'tag', 'dummy'); 
axis([xl zl]); axis ij;
set(gca,'pos',[0.1 0.12 0.8 0.8],'tag','modelaxes');
set(gca,'nextplot','add');
ylabel('Depth (m)'); 
xlabel('Distance (m)')
%%%%%   Import model
nbody = fscanf(fid,' %d ',1);                         % read # bodies
for ib=1:nbody,
    x = [];
    z = [];
    mark1 = [];
    %dummy = fscanf(fid,' %s \n',1);                   % read tag
    dummy = fgetl(fid);                                % read tag
    bodies{ib,4} = double(dummy);
    dummy = fscanf(fid,' %d %f %f %f %f \n',5);
    nmark = dummy(1);                                 % # of vertices 
    bodies{ib,3} = dummy(2:5)';                       % object properties
    for iv = 1:nmark,
        dummy = fscanf(fid,' %f  %f \n',2);
        x(iv)  = dummy(1);                            % read vertices
        z(iv)  = dummy(2);
        mark1(iv) = plot(x(iv),z(iv),'sb','era','back');   %#ok<AGROW>
        set(mark1(iv),'MarkerSize',[4],'LineWidth',[1.0])
    end
    bodies{ib,1} = x';
    bodies{ib,2} = z';
    patch3 = patch(x,z,[rand(1,1), rand(1,1) rand(1,1)],'edgealpha',0);
% Store new handles         
    handles{ib,1} = mark1(1:nmark);
    handles{ib,2} = patch3;
% draw the tag
    [tx, tz] = centroid(x,z);
    txt = text(tx,tz,char(bodies{ib,4}),'horizontalal','center', ... 
        'fontweight','bold','fontsize',9,'erase','xor');
    handles{ib,3} = txt;
% draw and hide object properties 
    dummy = bodies{ib,3};
    str = ['r=' num2str(dummy(4)) ', K=' num2str(dummy(1)) ...
        ', m=' num2str(dummy(3))];
    txt = text(tx,tz,str,'horizontalal','center', ...
        'fontname','symbol','fontweight','bold','fontsize',10,...
        'erase','xor','visible','off');
    handles{ib,4} = txt;
end
fclose(fid);
backup  = {};                                % initialize backup cell array
model_was_saved = 1;
return
%%% END HANDLE FUNCTION readmodel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savemodel(vhandle,eventdata) 
%%%%%   Save current model to ascii file
global ENVAR
global bodies handles nbody inpmod backup model_was_saved
if nbody == 0,
    helpdlg('THE MODEL IS STILL BLANK! NOTHING TO DO!','HELP');
    uiwait
    return
end
modax = findobj('tag','modelaxes');   % Ensure model axes is current object
xl = get(modax,'xlim');               % and get axes bounds
zl = get(modax,'ylim');
antenna=str2num(get(findobj('tag','antennaf'),'String'));  % get frequency
[mfname, mpname]= uiputfile('*.dat','Give model file name',...
    ENVAR.currentworkdir);
if mfname == 0, 
    return 
end;
idat = findstr(mfname,'.dat');
if isempty(idat),
    mfname = [mfname '.dat'];
end
fid = fopen([mpname mfname],'w');                 % and open O/P file
fprintf(fid,' %f \n',antenna);                    % write antenna frequency
fprintf(fid,' %f %f %f %f \n',[xl zl]);           % write axis bounds
fprintf(fid,' %d \n',nbody);                      % write # bodies
for ib=1:nbody,
    x = bodies{ib,1};
    z = bodies{ib,2};
    nmark = length(x);
    fprintf(fid,' %s \n',char(bodies{ib,4}));     % write tag
    fprintf(fid,' %d %f %f %f %f \n',[nmark bodies{ib,3}]); % properties
    for iv = 1:nmark,
        fprintf(fid,' %f  %f \n',[x(iv) z(iv)]);  % write vertices
    end
end
fclose(fid);
model_was_saved = 1;
return
%%% END HANDLE FUNCTION savemodel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quitprg(vhandle,eventdata) 
%%%%%   Terminates program
global bodies handles nbody inpmod backup model_was_saved
if nbody == 0 | isempty(findobj('tag','modelaxes')), 
    modf = findobj('tag','modelfigure'); 
    delete(modf); 
    clear global bodies handles inpmod nbody backup model_was_saved; 
    return; 
end; 
if ~model_was_saved,
    qbut = questdlg(['The model has not been saved! ' ...
        'Do you wish to continue?'], 'WARNING','YES','NO','NO');
    if isequal(qbut,'NO'), 
        return
    end
end
modf = findobj('tag','modelfigure'); 
delete(modf); 
clear global modf bodies handles inpmod nbody backup model_was_saved
return
%%% END HANDLE FUNCTION quitprg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inputmode(vhandle,eventdata,mode) 
%%%%%   Select input mode
global inpmod 
modax = findobj('tag','modelaxes');
if strcmp(mode,'fingers'),
    inpmod = 1;
    set(get(modax,'title'), 'string', 'Input Mode: Use Fingers')
elseif strcmp(mode,'usepointer'),
    inpmod = 2;
    set(get(modax,'title'), 'string', 'Input Mode: Use Pointer')
elseif strcmp(mode,'readfile'),
    inpmod = 3;
    set(get(modax,'title'), 'string', 'Input Mode: From File')
    erh = warndlg(['Note: This input mode applies only to ' ...
                    'POLYGONAL objects. To enter circular objects,' ...
                    'select one of the other modes.'],'NOTE');
    uiwait(erh)
else
    erh = errordlg('Wrong choise! Default set to use pointer','ERROR');
    uiwait(erh)
    inpmod = 2;
end
return
%%% END HANDLE FUNCTION inputmode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function addpolygon(vhandle,eventdata) 
%%%%%   ADD A BODY OF POLYGONAL CROSS SECTION
global ENVAR
global bodies handles nbody inpmod backup model_was_saved
x       = [];                            
z       = [];                            
mark1   = [];     
lin2    =[];
modax = findobj('tag','modelaxes');   % Ensure model axes is current object
if isempty(modax),
    helpdlg('THE MODEL IS BLANK! NOTHING TO DO!','HELP');
    uiwait
    return
end
xl = get(modax,'xlim');
zl = get(modax,'ylim');
dx = diff(xl);
dz = diff(zl);
antenna=str2num(get(findobj('tag','antennaf'),'String'));  % get frequency
if isempty(antenna),
     errordlg('Please, first set the antenna frequency','ERROR');
     uiwait
     return
end
if antenna > 1e6,                         % Ensure frequency in MH
    antenna = antenna/1e6; 
end; 
nbody = nbody + 1;                        % the new object can now be set
%%%%%   Type polygon vertices
if inpmod == 1,
    answer = inputdlg('Give Number of Vertices',['Body ' num2str(nbody)],1);
    if isempty(answer), nbody = nbody-1; return; end;      % Cancelled
    nmark  = str2num(answer{1}); 
    x = zeros(nmark,1);
    z = zeros(nmark,1);
    clb = cell(1,1); 
    for i=1:nmark,
        clb = cellstr(['Vertex ' num2str(i) ... 
            ': Give X and Z co-ordinates                 .']);
        answer = inputdlg(clb,['Body ' num2str(nbody) ...
            ': Vertex co-ordinates'],1 );     
        if isempty(answer),                                % Cancelled
            nbody = nbody-1; 
            return; 
        end;
        lb = char(answer);                                          
        comma = findstr(lb(1,:),',');
        if ~isempty(comma), 
            for j=1:length(comma); 
                lb(1,comma) = '.'; 
            end;
            clear comma
        end
        dummy = str2num(lb(1,:));
        x(i) = dummy(1); 
        z(i) = dummy(2);
    end
%%%%%   Define polygonal object using Pointer
elseif inpmod == 2,
%%%% Display help message %%%
    msgtxt = cell(1);
    msgtxt(1) = cellstr('To set vertices click the LEFT mouse button. ');
    msgtxt(2) = cellstr('To correct mistakes click the MIDDLE button. ');
    msgtxt(3) = cellstr('To finish, click the RIGHT mouse button! ');
    msg = msgbox(msgtxt);
    pos=get(msg,'position');
    set(msg,'pos',[20 20 pos(3) pos(4)]);
    pause(0.5)
%%% make sure that right figure is on focus
    figure(findobj('tag','modelfigure'));  
    pause(0.5)
    [x, z] = setpolygon('meters', 'meters');
    nmark = size(x,1);
    if exist('msg') && ishandle(msg),              % kill help message
        delete(msg);
    end
%%%%%   Load polygonal object from file
elseif inpmod == 3,
    [pfname, ppname]= uigetfile('*.dat;*.txt',...
        'Please give polygonal object''s file name', ENVAR.currentworkdir);
    if pfname == 0,                           % Canceled
        return;  
    end; 
    poly = importdata([ppname pfname]);       % get x,z vertices of polygon
    x = poly(:,1);
    z = poly(:,2);
    nmark = size(x,1);
    clear global poly
end                                           % inpmod
%%% Process the object
%%% If x(i) or z(i) are very close to, or beyond axes limits, snap them to 
%%% axes limits by linear interpolation / extrapolation. Do not bother if
%%% inpmod = 2, as this has been taken care of by "setpolygon". 
%%%%%%%%%%%%
if isempty(x),
    nbody = nbody-1;
    return
end
if inpmod ~= 2,
    [x, z ] = snappolytoaxes(x,z,xl,zl);
end
%%% If [x(i), z(i)] are very close to an existing vertex of an adjacent
%%% polygon, snap them together.
for i = 1:nmark
    if nbody > 1,
        for ib = 1:nbody-1,
            testx = bodies{ib,1};
            testz = bodies{ib,2};
            d = abs((testx-x(i))/dx) + abs((testz-z(i))/dz);
            %d = sqrt( (testx - x(i)).^2 + (testz - z(i)).^2);
            [dmin,jmin] = min(d);
            if dmin < 0.01,
                x(i) = testx(jmin);
                z(i) = testz(jmin);
            end
        end
    end
%%% Plot the outline of the polygon            
    mark1(i) = plot(x(i),z(i),'sb','markerfacecolor','w');
    set(mark1(i),'MarkerSize',4,'LineWidth',1.0)
end
%%% Draw the new polygonal body
patch3 = patch(x,z,[rand(1,1), rand(1,1) rand(1,1)],'edgealpha',0);
%%%%%   Fill cell-arrays with vertices of polygonal bodies
bodies{nbody,1}  = x(1:nmark);
bodies{nbody,2}  = z(1:nmark);
%%%%%   Take object properties
if isempty(bodies{nbody,3}),
      bodyprops = getprops(antenna,bodies{nbody,3});
      if isempty(bodyprops),    % Last minute decision to abandon this one
          delete(mark1(1:nmark));
          delete(patch3);
          bodies(nbody,:) = [];
          backup = {};
          nbody = nbody - 1;
          return
     end
end
bodies{nbody,3}  = bodyprops;
handles{nbody,1} = mark1(1:nmark);
handles{nbody,2} = patch3;
set(modax,'xlim',xl,'ylim',zl)                          % force axes limits 
% Give a name to the object
deftag = cellstr(['Poly' num2str(nbody)]);
answer = inputdlg('Please give a name (tag) to this object','REQUEST',...
    1,deftag );  
if isempty(answer),
    bodies{nbody,4} = double(deftag{1});
else
    bodies{nbody,4} = double(answer{1}); 
end
% Project the tag
[tx, tz] = centroid(x(1:nmark),z(1:nmark));
txt = text(tx,tz,answer{1},'horizontalal','center', ... 
     'fontweight','bold','fontsize',9,'erase','xor');
handles{nbody,3} = txt;
% draw and hide object properties 
str = ['r=' num2str(bodyprops(4)) ', K=' num2str(bodyprops(1)) ...
    ', m=' num2str(bodyprops(3))];
txt = text(tx,tz,str,'horizontalal','center', ...
     'fontname','symbol','fontweight','bold','fontsize',10,...
     'erase','xor','visible','off');
handles{nbody,4} = txt;
backup = {};               % Re-initialize backup - model structure changed
model_was_saved = 0;
return
%%% END HANDLE FUNCTION addpolygon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function addcircle(vhandle,eventdata) 
%%%%%   ADD A BODY OF CIRCULAR CROSS SECTION
global bodies handles nbody inpmod backup model_was_saved
x       = [];
z       = [];
mark1   = [];     
lin2    =[];
modax = findobj('tag','modelaxes');   % Ensure model axes is current object
if isempty(modax),
    helpdlg('THE MODEL IS BLANK! NOTHING TO DO!','HELP');
    uiwait
    return
end
xl = get(modax,'xlim');
zl = get(modax,'ylim');
dx = diff(xl);
dz = diff(zl);
antenna=str2num(get(findobj('tag','antennaf'),'String'));  % get frequency
if isempty(antenna),
     errordlg('Please, first set the antenna frequency','ERROR');
     uiwait
     return
end
if antenna > 1e6,                             % Ensure frequency in MHz
    antenna = antenna/1e6; 
end;
nbody = nbody + 1;                            % the new object can now be set
%%%%%   Type centre and radius data 
if inpmod == 1,
    clb = cell(3,1);                                      
    clb(1) = cellstr(['Give X-coordinate of the CENTRE']);
    clb(2) = cellstr(['Give Z-coordinate of the CENTRE']);
    clb(3) = cellstr(['Give the RADIUS of the object']);
    for i=1:3
        cldef(i) = cellstr(num2str(0));
    end
    answer = inputdlg(clb,['Body ' num2str(nbody) ' (Circle)'],1,cldef );
    if isempty(answer), 
        nbody = nbody-1; 
        return, 
    end;
    lb = char(answer);                                                   
    for i=1:3
        comma = findstr(lb(i,:),',');
        if ~isempty(comma), 
            for j=1:length(comma); 
                lb(i,comma) = '.'; 
            end;
            clear comma
        end
    end
    xc     = str2num(lb(1,:));
    zc     = str2num(lb(2,:));
    radius = str2num(lb(3,:));
%%%%%   Define circle using mouse and getpoint
elseif inpmod == 2,
    msgtxt = cell(1);
    msgtxt(1) = cellstr('To set the circular object click at the CENTRE ');
    msgtxt(2) = cellstr('and at distance RADIUS from the CENTRE');
    msg = msgbox(msgtxt);
    pos=get(msg,'position');
    set(msg,'pos',[20 20 pos(3) pos(4)]);
    pause(0.5)
%%% make sure that right figure is on focus
    figure(findobj('tag','modelfigure'));  
    pause(0.5)
%%%
    [xc, zc] = getpoint(' meters', ' meters');
    markc = plot(xc,zc,'kh','markerfacecolor','w','LineWidth',1.0);
    [xr, zr] = getpoint(' meters', ' meters');
    radius = sqrt((xr-xc)^2 + (zr-zc)^2);
%%%    
    if exist('msg') & ishandle(msg),              % kill help message
        delete(msg);
    end
end                                               % inpmod
%%%%%   Compute circle as polygon            
add   = 2*pi/24.0;    
alf   =[ 0 : add : 2*pi]';
nmark = length(alf);
x     = xc + radius*sin(alf);
z     = zc + radius*cos(alf);
%%% If x(i) or z(i) are beyond axes limits, snap them to axes by linear
%%% interpolation or extrapolation.  
[x, z ] = snappolytoaxes(x,z,xl,zl);
%%% Draw polygon
if exist('markc') && ishandle(markc),
    delete(markc);
end
for i=1:nmark
    mark1(i) = plot(x(i),z(i),'sb','era','back');  
    set(mark1(i),'MarkerSize',[3],'LineWidth',[1.0])
end;
patch3 = patch(x,z,[rand(1,1), rand(1,1) rand(1,1)],'edgealpha',0);
%%%%%   Fill cell-arrays with object data
bodies{nbody,1}  = x(1:nmark)';
bodies{nbody,2}  = z(1:nmark)';
%   Take object properties
if isempty(bodies{nbody,3}),                
    bodyprops = getprops(antenna,bodies{nbody,3});
    if isempty(bodyprops),   % Last minute decision to abandon object
        delete(mark1(1:nmark));
        delete(patch3);
        bodies(nbody,:) = [];
        backup = {};
        nbody = nbody - 1;
        return
    end
end
bodies{nbody,3}  = bodyprops;
handles{nbody,1} = mark1(1:nmark);
handles{nbody,2} = patch3;
set(modax,'xlim',xl,'ylim',zl)                          % force axes limits 
% Give a name to the object
deftag = cellstr(['Circ ' num2str(nbody)]);
answer = inputdlg('Please give a name (tag) to this object','REQUEST',1,deftag );  
if isempty(answer),
    bodies{nbody,4} = double(deftag{1});
else
    bodies{nbody,4} = double(answer{1}); 
end
% and project it
[tx, tz] = centroid(x,z); 
txt = text(tx,tz,answer{1},'horizontalal','center', ... 
      'fontweight','bold','fontsize',9,'erase','xor');
handles{nbody,3} = txt;
% draw and hide object properties 
    str = ['r=' num2str(bodyprops(4)) ', K=' num2str(bodyprops(1)) ...
            ', m=' num2str(bodyprops(3))];
    txt = text(tx,tz,str,'horizontalal','center', ...
        'fontname','symbol','fontweight','bold','fontsize',10,...
        'erase','xor','visible','off');
    handles{nbody,4} = txt;
backup = {};               % Re-initialize backup - model structure changed
model_was_saved = 0;
return
% END HANDLE FUNCTION addcircle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delbody1(vhandle,eventdata) 
%%%%%   UNDO (DELETE) A BODY
global bodies handles nbody inpmod backup model_was_saved
if nbody == 0 || isempty(findobj('tag','modelaxes')),
    helpdlg('THE MODEL IS STILL BLANK! NOTHING TO DO!','HELP');
    uiwait
    return
end
str = cell(0);
for ib=1:nbody,
    str = [str cellstr(['UNDO ' char(bodies{ib,4})])];
end
str = [ str  cellstr('Cancel')];
uicontrol('Style', 'popup', 'backgroundcolor','y', 'tag', 'deletebody', ...
          'String', str, 'units','normalized',...
          'Position', [0.1 0.93 0.2 0.05], 'Callback', @delbody2); 
return 
% END HANDLE FUNCTION delbody1

function delbody2(vhandle,eventdata) 
%%%%%   UNDO (DELETE) A BODY
global bodies handles nbody inpmod backup model_was_saved
delbody = findobj('tag','deletebody');  
ib=get(delbody,'value'); 
if ib == nbody + 1; 
    delete(delbody);  
    delete(findobj('tag','deltxt')); 
    return; 
elseif ib <= nbody; 
      bodies(ib,:) = []; 
      delete(handles{ib,1}); delete(handles{ib,2});  
      delete(handles{ib,3}); delete(handles{ib,4}); 
      handles(ib,:) = []; 
      nbody = nbody -1; 
      delete(delbody); 
      delete(findobj('tag','deltxt')); 
end; 
backup = {};              % Re-initialize backup - model structure changed
model_was_saved = 0;
return 
% END HANDLE FUNCTION delbody2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function reviewproperties1(vhandle,eventdata) 
%%%%%   CHANGE BODY PROPERTIES
global bodies handles nbody inpmod backup model_was_saved
if nbody == 0 || isempty(findobj('tag','modelaxes')),
    helpdlg('THE MODEL IS STILL BLANK! NOTHING TO DO!','HELP');
    uiwait
    return
end
%%%%%   Store current model state in backup arrays for future UNDO's
if exist('backup'),                             
    nbk = length(backup); 
    backup{nbk+1} = bodies;
end
str = cell(0);
for ib=1:nbody,
    str = [str cellstr(['Change ' char(bodies{ib,4})])];
end
str = [ str  cellstr('Cancel')];
uicontrol('Style', 'popup', 'backgroundcolor','y', 'tag', 'chbps', ...
          'String', str, 'units','normalized',...
          'Position', [0.1 0.93 0.2 0.05], 'Callback', @reviewproperties2); 
return
% END HANDLE FUNCTION reviewproperties1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function reviewproperties2(vhandle,eventdata) 
%%%%%   CHANGE BODY PROPERTIES
global bodies handles nbody inpmod backup model_was_saved
chbp = findobj('tag','chbps');  
ib=get(chbp,'value'); 
if ib == nbody + 1; 
    delete(chbp);  
    delete(findobj('tag','deltxt')); 
    return; 
elseif ib <= nbody; 
    antenna=str2num(get(findobj('tag','antennaf'),'String')); 
    bodyprops = getprops(antenna,bodies{ib,3}); 
    if isempty(bodyprops),          
        delete(chbp); 
        delete(findobj('tag','deltxt')); 
        return; 
    end; 
    bodies{ib,3} = bodyprops; 
    dummy = bodies{ib,3}; 
    [tx,tz] = centroid(bodies{ib,1},bodies{ib,2}); 
    str = ['r=' num2str(dummy(4)) ', K=' num2str(dummy(1)) ...
        ', m=' num2str(dummy(3))]; 
    txt = text(tx,tz,str,'horizontalal','center','fontname','symbol'); 
    set(txt,'fontweight','bold','fontsize',10,'erase','xor',...
        'visible','off'); 
    delete(handles{ib,4}); 
    handles{ib,4} = txt; 
    delete(chbp); 
    delete(findobj('tag','deltxt')); 
    model_was_saved = 0;
end; 
return 
% END HANDLE FUNCTION reviewproperties2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function relocatevertex(vhandle,eventdata) 
%%%%%   RELOCATE A VERTEX
global bodies handles nbody inpmod backup model_was_saved
modax = findobj('tag','modelaxes');  % Ensure model axes is current object
if nbody == 0 | isempty(modax),
    helpdlg('THE MODEL IS STILL BLANK! NOTHING TO DO!','HELP');
    uiwait
    return
end
xl = get(modax,'xlim');
zl = get(modax,'ylim');
dx = diff(xl);
dz = diff(zl);
% Click on a point near the vertex
[xp, zp] = getpoint(' meters', ' meters');
% Determine the # of vertices near clicked point
if nbody > 0,
    nvfound = 0;
    ipsv  = [];
    jvsv  = [];
    for ib = 1:nbody,
        testx = bodies{ib,1};
        testz = bodies{ib,2};
        d = abs((testx-xp)/dx) + abs((testz-zp)/dz);
        %d = sqrt( (testx - xp).^2 + (testz - zp).^2);
        [dmin,jmin] = min(d);
% Flag vertices found
        if dmin < 0.01,
            nvfound = nvfound + 1;
            ipsv(nvfound)=ib;
            jvsv(nvfound)=jmin;
            xpv = testx(jmin);
            zpv = testz(jmin);
        end
    end
end
if nvfound == 0, 
    errordlg('No vertices found near selected point','ERROR');    
    uiwait
    return
end
%disp(['Found ' num2str(nvfound) ' vertices at [' num2str(xpv), ', '...
%num2str(zpv) ' ]']);     
%%%%%   Store current model state in backup arrays for future UNDO's
if exist('backup'),                             
    nbk = length(backup); 
    backup{nbk+1} = bodies;
end
if inpmod == 1 || inpmod == 3,
    clb = cell(2,1);                                      
    clb(1) = cellstr(['Give NEW X-coordinate of the vertex']);
    clb(2) = cellstr(['Give NEW Z-coordinate of the vertex']);
    cldef = {num2str(xpv); num2str(zpv)};
    answer = inputdlg(clb,'Give NEW location',1,cldef );
    if isempty(answer), 
        return, 
    end;
    lb = checkcomma(answer);                                                   
    xv = str2num(lb(1,:));
    zv = str2num(lb(2,:));
    
elseif inpmod == 2,
  
% Click on new vertex position      
    pause(0.5)
    [xv, zv] = getpoint(' meters', ' meters');
end

% Copy this position 
for k=1:nvfound,
    tmp = bodies{ipsv(k),1};
    tmp(jvsv(k)) = xv;
    bodies{ipsv(k),1} = tmp;
    tmp = bodies{ipsv(k),2};
    tmp(jvsv(k)) = zv;
    bodies{ipsv(k),2} = tmp;
end
% Redraw the model
for ib=1:nbody,
% Clear ib'th object and delete its handles but first save patch colours
    delete(handles{ib,1});
    clr = get(handles{ib,2},'facecolor');
    delete(handles{ib,2});
    delete(handles{ib,3});
    delete(handles{ib,4});
% Redraw the ib'th object
    x = bodies{ib,1};
    z = bodies{ib,2};
    nmark = length(x);
    for i=1:nmark,
        mark1(i) = plot(x(i),z(i),'sb','era','back');  
        set(mark1(i),'MarkerSize',[4],'LineWidth',[1.0])
    end;
    patch3 = patch(x,z,clr,'edgealpha',0);
% Store new handles         
    handles{ib,1} = mark1(1:nmark);
    handles{ib,2} = patch3;
% Restore axes limits 
    set(modax,'xlim',xl,'ylim',zl)      
% Draw the tag
    [tx, tz] = centroid(x,z);
    txt = text(tx,tz,char(bodies{ib,4}),'horizontalal','center', ... 
        'fontweight','bold','fontsize',9,'erase','xor');
    handles{ib,3} = txt;
% draw and hide object properties 
    dummy = bodies{ib,3};
    str = ['r=' num2str(dummy(4)) ', K=' num2str(dummy(1)) ', m=' ...
            num2str(dummy(3))];
    txt = text(tx,tz,str,'horizontalal','center', ...
        'fontname','symbol','fontweight','bold','fontsize',10,'erase',...
        'xor','visible','off');
    handles{ib,4} = txt;
end
model_was_saved = 0;
return
% END HANDLE FUNCTION relocatevertex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function addvertex(vhandle,eventdata) 
%%%%%   ADD A VERTEX
global bodies handles nbody inpmod backup model_was_saved
modax = findobj('tag','modelaxes');   % Endure model axes is current object
if nbody == 0 || isempty(modax),
    helpdlg('THE MODEL IS STILL BLANK! NOTHING TO DO!','HELP');
    uiwait
    return
end
xl = get(modax,'xlim');
zl = get(modax,'ylim');
dx = diff(xl);
dz = diff(zl);
if inpmod == 1 || inpmod == 3,
    qbut = questdlg( ...
        'Click on the two vertices ADJACENT to the new point. ', ...
        'INFO','OK','OK');
elseif inpmod == 2,
    qbut = questdlg( ...
        ['Click on the two vertices ADJACENT to the new point. ' ...
         'Then click on the location of the NEW vertex.'],'INFO','OK','OK');
end
pause(0.5)
% Click on the first vertex and determine the # of bodies sharing vertex
[xp, zp] = getpoint(' meters', ' meters');
nvfound1 = 0;
ipsv1  = [];
for ib = 1:nbody,
    x = bodies{ib,1};
    z = bodies{ib,2};
    d = abs((x-xp)/dx) + abs((z-zp)/dz);
    [dmin,jmin] = min(d);
    if dmin < 0.01,
        nvfound1 = nvfound1 + 1;
        ipsv1(nvfound1)=ib;
        xp1 = x(jmin);
        zp1 = z(jmin);
    end
end
if nvfound1 == 0, 
    errordlg('No vertices found near the first selected point','ERROR');
    return
end
% Click on the second vertex and determine the # of bodies sharing vertex
[xp, zp] = getpoint(' meters', ' meters');
nvfound2 = 0;
ipsv2  = [];
for ib = 1:nbody,
    x = bodies{ib,1};
    z = bodies{ib,2};
    d = abs((x-xp)/dx) + abs((z-zp)/dz);
    [dmin,jmin] = min(d);
    if dmin < 0.01,
        nvfound2 = nvfound2 + 1;
        ipsv2(nvfound2)=ib;
        xp2 = x(jmin);
        zp2 = z(jmin);
    end
end
if nvfound2 == 0,
    errordlg('No vertices found near second selected point','ERROR');
    return
end
% The side { (xp1,zp1), (xp2,zp2) } will be on the intersection of bodies
% ipsv1 and ipsv2 
ipsv = intersect(ipsv1,ipsv2);
if isempty(ipsv),
    errordlg('Selected vertices do not form the side of a polygon','ERROR');
    return
end
nvfound = length(ipsv);

if inpmod == 1 || inpmod == 3,
    clb = cell(2,1);                                      
    clb(1) = cellstr(['Give X-coordinate of the NEW vertex']);
    clb(2) = cellstr(['Give Z-coordinate of the NEW vertex']);
    answer = inputdlg(clb,'Give NEW vertex location',1 );
    if isempty(answer), 
        return, 
    end;
    lb  = checkcomma(answer);                                                   
    xp3 = str2num(lb(1,:));
    zp3 = str2num(lb(2,:));
    
elseif inpmod == 2,
  
    % Click on position of new vertex 
    %[xp3,zp3]=ginput(1);
    [xp3, zp3] = getpoint(' meters', ' meters');
end

% If [xp3 zp3] are very close to an existing vertex of a polygon, snap them
% together 
for ib = 1:nbody-1,
    x = bodies{ib,1};
    z = bodies{ib,2};
    d = abs((x-xp3)/dx) + abs((z-zp3)/dz);
    %d = sqrt( (x - xp3).^2 + (z - zp3).^2);
    [dmin,jmin] = min(d);
    if dmin < 0.01,
        xp3 = x(jmin);
        zp3 = z(jmin);
    end
end
%%%%%   Store current model state in backup arrays for future UNDO's
if exist('backup'), 
    nbk = length(backup); 
    backup{nbk+1} = bodies;
end
% Now insert the new vertex to polygons ipsv(k);   
for k = 1:nvfound
    x = bodies{ipsv(k),1};
    z = bodies{ipsv(k),2};
    nn = length(x);
    ii = find(x==xp1 & z==zp1);
    jj = find(x==xp2 & z==zp2);
% Insert (xp3, zp3) in ipsv(k) object
    if ii == 1 && jj == nn,
        x = [x; xp3];
        z = [z; zp3];
    elseif ii == nn & jj == 1,
        x = [xp3; x];
        z = [zp3; z];
    elseif abs(ii-jj)~=1,
        errordlg(['DID NOT CLICK ON ADJACENT VERTICES OF OBJECT ' ...
            num2str(ipsv(k))],'ERROR'); 
% Since action is aborted, restore backup array to its previous state
        nbk = length(backup); 
        backup = backup(1:nbk-1);   
        return
    elseif ii < jj,
        x = [x(1:ii); xp3; x(jj:nn)];
        z = [z(1:ii); zp3; z(jj:nn)];
    elseif ii > jj
        x = [x(1:jj); xp3; x(ii:nn)];
        z = [z(1:jj); zp3; z(ii:nn)];
    end
    bodies{ipsv(k),1} = x;
    bodies{ipsv(k),2} = z;
end
% Redraw the model
for ib=1:nbody,
% Clear ib'th object and delete its handles but first save patch colours
    delete(handles{ib,1});
    clr = get(handles{ib,2},'facecolor');
    delete(handles{ib,2});
    delete(handles{ib,3});
    delete(handles{ib,4});
% Redraw the ib'th object
    x = bodies{ib,1};
    z = bodies{ib,2};
    nmark = length(x);
    for i=1:nmark,
        mark1(i) = plot(x(i),z(i),'sb','era','back');  
        set(mark1(i),'MarkerSize',[4],'LineWidth',[1.0])
    end;
    patch3 = patch(x,z,clr,'edgealpha',0);
% Store new handles 
    handles{ib,1} = mark1(1:nmark);
    handles{ib,2} = patch3;
% Set axes limits 
    set(modax,'xlim',xl,'ylim',zl)  
% Draw tag   
    [tx, tz] = centroid(x,z);
    txt = text(tx,tz,char(bodies{ib,4}),'horizontalal','center', ... 
         'fontweight','bold','fontsize',9,'erase','xor');
    handles{ib,3} = txt;
% draw and hide object properties 
    dummy = bodies{ib,3};
    str = ['r=' num2str(dummy(4)) ', K=' num2str(dummy(1)) ...
        ', m=' num2str(dummy(3))];
    txt = text(tx,tz,str,'horizontalal','center', ...
        'fontname','symbol','fontweight','bold','fontsize',10,...
        'erase','xor','visible','off');
    handles{ib,4} = txt;
end
model_was_saved = 0;
return
% END HANDLE FUNCTION addvertex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function removevertex(vhandle,eventdata) 
%%%%%   REMOVE A VERTEX
global bodies handles nbody inpmod backup model_was_saved
modax = findobj('tag','modelaxes');   % Ensure model axes is current object
if nbody == 0 || isempty(modax),
    helpdlg('THE MODEL IS STILL BLANK! NOTHING TO DO!','HELP');
    uiwait
    return
end
xl = get(modax,'xlim');
zl = get(modax,'ylim');
dx = diff(xl);
dz = diff(zl);
% Click on a point near the vertex
pause(0.5)
[xp, zp] = getpoint(' meters', ' meters');
% Determine the # of vertices near clicked point
if nbody > 0,
    nvfound = 0;
    ipsv  = [];
    jvsv  = [];
    for ib = 1:nbody,
        testx = bodies{ib,1};
        testz = bodies{ib,2};
        d = abs((testx-xp)/dx) + abs((testz-zp)/dz);
        [dmin,jmin] = min(d);
        if dmin < 0.01,
            nvfound = nvfound + 1;
            ipsv(nvfound)=ib;
            jvsv(nvfound)=jmin;
        end
    end
end
if nvfound == 0, 
    errordlg('No vertices found near selected point','ERROR');
    return
end
if nvfound > 2, 
    errordlg('Cannot remove a junction of multiple bodies','ERROR');
    return
end
%%%%%   Store current model state in backup arrays for future UNDO's
if exist('backup'), 
    nbk = length(backup); 
    backup{nbk+1} = bodies;
end
% Clear selected vertex
for k=1:nvfound,
    testx = bodies{ipsv(k),1};
    testz = bodies{ipsv(k),2};
    nn = length(testx);
    if jvsv(k) == 1,
        testx = testx(2:nn);
        testz = testz(2:nn);
    elseif jvsv(k) == nn,
        testx = testx(1:nn-1);
        testz = testz(1:nn-1);
    else
        testx = [testx(1:jvsv(k)-1); testx(jvsv(k)+1:nn)];
        testz = [testz(1:jvsv(k)-1); testz(jvsv(k)+1:nn)];
    end
    bodies{ipsv(k),1} = testx;
    bodies{ipsv(k),2} = testz;
end
% Redraw the model
for ib=1:nbody,
% Clear ib'th object and delete its handles but first save patch colours
    delete(handles{ib,1});
    clr = get(handles{ib,2},'facecolor');
    delete(handles{ib,2});  
    delete(handles{ib,3});
    delete(handles{ib,4});
% Redraw the ib'th object
    x = bodies{ib,1};
    z = bodies{ib,2};
    nmark = length(x);
    for i=1:nmark,
        mark1(i) = plot(x(i),z(i),'sb','era','back');  
        set(mark1(i),'MarkerSize',[4],'LineWidth',[1.0])
    end;
    patch3 = patch(x,z,clr,'edgealpha',0);
% Store new handles 
    handles{ib,1} = mark1(1:nmark);
    handles{ib,2} = patch3;
% Set axes limits 
    set(modax,'xlim',xl,'ylim',zl)  
% Draw tag
    [tx, tz] = centroid(x,z);
    txt = text(tx,tz,char(bodies{ib,4}),'horizontalal','center', ... 
        'fontweight','bold','fontsize',9,'erase','xor');
    handles{ib,3} = txt;
% draw and hide object properties 
    dummy = bodies{ib,3};
    str = ['r=' num2str(dummy(4)) ', K=' num2str(dummy(1)) ...
            ', m=' num2str(dummy(3))];
    txt = text(tx,tz,str,'horizontalal','center', ...
        'fontname','symbol','fontweight','bold','fontsize',10,...
        'erase','xor','visible','off');
    handles{ib,4} = txt;
end
model_was_saved = 0;
return
% END HANDLE FUNCTION removevertex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function undochanges(vhandle,eventdata) 
%%%%%   UNDO VERTEX OPERATIONS
global bodies handles nbody inpmod backup model_was_saved
if isempty(backup) || isempty(backup{1}),
    helpdlg('NOTHING DONE TO UNDO!','HELP');
    uiwait
    return
end
modax = findobj('tag','modelaxes');  % Ensure model axes is current object
xl = get(modax,'xlim');
zl = get(modax,'ylim');
clr = [];                            % Save patch colours
nbb = size(handles,1); 
for ib=1:nbb,
    clr(ib,1:3) = get(handles{ib,2},'facecolor');
end
delete(handles{:,1});                % delete handles but do not clear the
delete(handles{:,2});                % variable name because it is global
delete(handles{:,3}); 
delete(handles{:,4}); 
%clear handles
nbk = length(backup);                % make previous model state current
bodies = backup{nbk};
nbody = size(bodies);    
nbody = nbody(1);
backup = backup(1:nbk-1);            % clear previous backup model state
for ib=1:nbody,                      % Redraw current model 
    x = bodies{ib,1};
    z = bodies{ib,2};
    nmark = length(x);
    for i=1:nmark,
        mark1(i) = plot(x(i),z(i),'sb','era','back');  
        set(mark1(i),'MarkerSize',[4],'LineWidth',[1.0])
    end;
    patch3 = patch(x,z,clr(ib,:),'edgealpha',0);
    handles{ib,1} = mark1(1:nmark);   
    handles{ib,2} = patch3;
    set(modax,'xlim',xl,'ylim',zl)      
    [tx, tz] = centroid(x,z);
    txt = text(tx,tz,char(bodies{ib,4}),'horizontalal','center', ... 
       'fontweight','bold','fontsize',9,'erase','xor');
    handles{ib,3} = txt;
% draw and hide object properties 
    dummy = bodies{ib,3};
    str = ['r=' num2str(dummy(4)) ', K=' num2str(dummy(1)) ... 
            ', m=' num2str(dummy(3))];
    txt = text(tx,tz,str,'horizontalal','center', ...
        'fontname','symbol','fontweight','bold','fontsize',10,...
        'erase','xor','visible','off');
    handles{ib,4} = txt;
end
model_was_saved = 0;
return
% END HANDLE FUNCTION undochanges

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function recolor(vhandle,eventdata) 
%%%%%   Recolor model
global bodies handles nbody inpmod backup model_was_saved
for ib=1:nbody; 
    set(handles{ib,2},'facecolor',[rand(1,1), rand(1,1) rand(1,1)]); 
end;
return
% END HANDLE FUNCTION recolor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hidemarkers(vhandle,eventdata) 
%%%%%   hide vertex markers
global bodies handles nbody inpmod backup model_was_saved
for ib=1:nbody; 
    set(handles{ib,1},'visible','off');
end; 
for ib=1:nbody; 
    set(handles{ib,2},'visible','off'); 
    set(handles{ib,2},'visible','on');
end;
return
% END HANDLE FUNCTION hidemarkers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showmarkers(vhandle,eventdata) 
%%%%%   show vertex markers
global bodies handles nbody inpmod backup model_was_saved
for ib=1:nbody; 
    set(handles{ib,1},'visible','on'); 
end;
return
% END HANDLE FUNCTION showmarkers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function solidmarkers(vhandle,eventdata) 
%%%%%   fill faces of vertex markers
global bodies handles nbody inpmod backup model_was_saved
for ib=1:nbody; 
    set(handles{ib,1},'markerfacecolor','w'); 
end;
return
% END HANDLE FUNCTION solidmarkers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hollowmarkers(vhandle,eventdata) 
%%%%%   empty faces of vertex markers
global bodies handles nbody inpmod backup model_was_saved
for ib=1:nbody; 
    set(handles{ib,1},'markerfacecolor','none'); 
end; 
for ib=1:nbody; 
    set(handles{ib,2},'visible','off'); 
    set(handles{ib,2},'visible','on'); 
end;
return
% END HANDLE FUNCTION hollowmarkers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hidetags(vhandle,eventdata) 
%%%%%   hide object tags
global bodies handles nbody inpmod backup model_was_saved
for ib=1:nbody; 
    set(handles{ib,3},'visible','off'); 
end;
return
% END HANDLE FUNCTION hidetags

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showtags(vhandle,eventdata) 
%%%%%   show object tags
global bodies handles nbody inpmod backup model_was_saved
for ib=1:nbody; 
    set(handles{ib,3},'visible','on'); 
end; 
return
% END HANDLE FUNCTION showtags

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showemprops(vhandle,eventdata) 
%%%%%   show object properties 
global bodies handles nbody inpmod backup model_was_saved
for ib=1:nbody; 
    set(handles{ib,3},'visible','off'); 
    set(handles{ib,4},'visible','on'); 
end; 
return
% END HANDLE FUNCTION showemprops

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hideemprops(vhandle,eventdata) 
%%%%%   hide object properties
global bodies handles nbody inpmod backup model_was_saved
for ib=1:nbody; 
    set(handles{ib,4},'visible','off'); 
    set(handles{ib,3},'visible','on'); 
end;
return
% END HANDLE FUNCTION hideemprops
