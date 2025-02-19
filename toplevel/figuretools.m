function figuretools()
%
% FIGURETOOLS : A collection of the most useful tools and figure
%               manipulation utilities that everybody would like to have
%               available, but without the full (and not allways useful)
%               complement of the MATLAB figure menu.  
%           ==> However, if one must access some utility from the MATLAB
%               Figure Menu, it is always possible to toggle it on / off.
%               FIGURETOOLS provides for:
%               1. Colour map editing (specific to the current figure).
%               2. Zooming, panning and data inspection.
%               3. Exporting the figure in various graphics formats.
%               4. Printing the current figure. 
%               5. Editing the current figure, axes and their children. 
%           >>> This program is part of the MATGPR suite.
%
%    Usage : figuretools;
% 
%   Author : Andreas Tzanis
%            Department of Geophysics, 
%            University of Athens
%            atzanis@geol.uoa.gr
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
%

matlabrelease = sscanf(version('-release'),'%d %f');

%%% Set up the parent menu
tools = uimenu('Label','绘图工具' );

%%% Insert Colorbar
uimenu(tools, 'Label', '插入颜色栏','Tag', 'cbar', ...
    'Callback', @mycbar)
%%% Edit colormap
uimenu(tools, 'Label', '编辑颜色图','Callback', 'colormapeditor;');
 
%%% ------ Zoom, Pan and Inspect data -------------------------------------
uimenu(tools, 'Label', '缩放', 'Tag', 'zoomon', 'Separator','on', ...
    'Callback', @myzoom);
uimenu(tools, 'Label', '取消缩放', ...
    'Callback', @myzoomout);

%%% Panning possible as from R14 and upward
if matlabrelease >= 14,
    uimenu(tools, 'Label', '平移', 'Tag', 'panon', 'Callback', @mypan);
end

%%% Inspection is possible with MATLAB's "datacursor" as from R14 an
%%% upward. For earlier releases use "getpoint.m' (which probably is more
%%% efficient).
if matlabrelease < 14,
    uimenu(tools, 'Label', '检查数据', 'Callback', 'getpoint;'); 
else
    uimenu(tools, 'Label', '检查数据', 'Tag', 'inspon', ...
        'Callback', @mydatacursor);
end

%%% ------ Save Figure, or copy to (Windows) clipboard --------------------
uimenu(tools, 'Label' , '保存为 ...', 'Tag', 'savefig', ... 
    'Separator', 'on', 'Callback', @mysaveas);
%%% Copy to windows clipboard
cpfig=uimenu(tools, 'Label', '复制图像');
uimenu(cpfig, 'Label', '到剪贴板', 'Callback', 'print -dbitmap');
% uimenu(cpfig, 'Label', 'to Enhanced Metafile', 'Callback', 'print -dmeta');
 
%%% ------ Printing utilities ---------------------------------------------
%%% Set printing options
uimenu(tools, 'Label', '页面设置', 'Separator', 'on', ...
    'callback', 'pagesetupdlg(gcf);'); 
%%% Set up printer
uimenu(tools, 'Label', '打印设置', ...
    'callback', 'printdlg(''-setup'',gcf);'); 
%%% Preview printout
uimenu(tools, 'Label', '打印预览', ...
    'callback', 'printpreview(gcf);'); 
%%% Print figure
uimenu(tools, 'Label', '打印 ', ...
    'callback', 'printdlg(''-crossplatform'',gcf);'); 

%%% ------Edit figure, axes properties / Edit and annotate objects --------
uimenu(tools, 'Label', '图像特征', 'Tag', 'figedt', ...
    'Separator', 'on', 'callback', 'propedit(gcf)'); 
uimenu(tools, 'Label', '轴特征', 'Tag', 'axedt', ...
    'callback', 'propedit(gca)'); 
%%% Edit plot
uimenu(tools, 'Label', '绘图编辑','Tag', 'plotedt', ... 
    'Callback', @myplotedit);

%%% If desired, recover the MATLAB figure menu for additional utilities 
uimenu(tools, 'Label', '隐藏图像菜单', 'Tag', 'figmenu', ...
     'Callback', @myfigmenu);

%%% You can also close the figure from this menu
uimenu(tools, 'Label', '关闭', 'Separator', 'on', ...
    'Callback', 'delete(gcf)');
 
     
return
%%% -----------------------------------------------------------------------
function mycbar(cbar,eventdata )
if strcmp(get(cbar,'Checked'),'off'),
     colorbar;
     set(cbar,'Checked','on');
else,
     colorbar off; 
     set(cbar,'Checked','off');
end
return
%%% -----------------------------------------------------------------------
function myzoom(zoomon,eventdata )
if strcmp(get(zoomon,'Checked'),'off'),
    zoom(gcf,'on');
    set(zoomon,'Checked','on'); 
else, 
    zoom(gcf,'off'); 
    set(zoomon,'Checked','off');
end; 
return
%%% -----------------------------------------------------------------------
function myzoomout(zoomon,eventdata )
zoom(gcf,'out');
if strcmp(get(findobj(gcf,'Tag','zoomon'),'Checked'),'on'),
    zoom(gcf,'off'); 
    set(findobj(gcf,'Tag','zoomon'),'Checked','off');
end; 
return
%%% -----------------------------------------------------------------------
function mypan(panon,eventdata);
if strcmp(get(panon,'Checked'),'off'),
    pan(gcf,'on');
    set(panon,'Checked','on'); 
else, 
    pan(gcf,'off'); 
    set(panon,'Checked','off');
end; 
return
%%% -----------------------------------------------------------------------
function mydatacursor(inspon,eventdata)
if strcmp(get(inspon,'Checked'),'off'),     
    datacursormode(gcf,'on');     
    dcm_obj = datacursormode(gcf);      
    set(dcm_obj,'displaystyle','window');    
    set(gcf,'toolbar','none');     
    set(inspon,'Checked','on');     
else,     
    datacursormode(gcf,'off');     
    dcm_obj = datacursormode(gcf);     
    clear dcm_obj;     
    set(inspon,'Checked','off');     
end; 
return
%%% -----------------------------------------------------------------------
function mysaveas(savefig,eventdata)
[filename, pathname, filterindex] = uiputfile( ...
     {'*.fig', 'MATLAB figure (*.fig)'; ...
     '*.ai','Adobe Illustrator `88 (*.ai)'; ...
     '*.bmp','Windows bitmap *.bmp)'; ...
     '*.emf','Enhanced metafile (*.emf)'; ...
     '*.eps','Encapsulated Postscript Level 1 (*.eps)'; ... 
     '*.jpg','JPEG image (*.jpg)'; ...
     '*.pbm','Portable bitmap (*.pbm)'; ... 
     '*.pcx','Paintbrush 24 bit (*.pcx)'; ...
     '*.pgm','Portable Graymap (*.pgm)'; ...
     '*.png','Portable Network Graphics (*.png)'; ...
     '*.ppm','Portable Pixmap (*.ppm)'; ...
     '*.tif','TIF image, compressed (*.tif)'}, ...
     'Save as','Untitled.fig'); 
     if filename==0,
         return
     end
     saveas(gcf,[pathname,filename]); 
return
%%% -----------------------------------------------------------------------
function myplotedit(plotedt,eventdata)
if strcmp(get(plotedt,'Checked'),'off'),
     plotedit(gcf,'on'); 
     set(plotedt,'Checked','on'); 
else,
     plotedit(gcf,'off'); 
     set(plotedt,'Checked','off'); 
end;
return
%%% -----------------------------------------------------------------------
function myfigmenu(figmenu,eventdata)
if strcmp(get(gcf,'menubar'),'figure'), 
     set(gcf,'menubar','none'); 
     set(figmenu,'Checked','on'); 
elseif strcmp(get(gcf,'menubar'),'none'), 
     set(gcf,'menubar','figure'); 
     set(figmenu,'Checked','off'); 
end
return
