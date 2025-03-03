function isosurf_gui(x,y,z,d3d)
%
% ISOSURF_GUI : GUI to control plots of plots of isometric GPR signal
%               intensity surfaces generated by ISOSURF.M
%         ==>   Used exclusively by ISOSURF.M 
% 
%Author : Andreas Tzanis,
%         Department of Geophysics, 
%         University of Athens
%         atzanis@geol.uoa.gr
%
% Copyright (C) 2008, Andreas Tzanis. All rights reserved.
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

% Create gui figure
figure('name','Iso-Surface Display Controls', ...
    'tag','iso3dcontrols', ...
    'position', [800   550   520   200], ...
    'numbertitle','off', ...
    'MenuBar','none',...
    'Resize', 'off', ...
    'Color',[0.76 0.86 0.69]);         % 'Color',[0 0.5 0.5]);
axis off;
fcolor = [0.76 0.87 0.78];             % fcolor = get(gcf,'color');
% data to create direction arrows
ar = ...
 [0 1 1 1 1
  0 0 1 1 1
  0 0 0 1 1
  0 0 0 0 1
  0 0 0 0 0
  0 0 0 0 1
  0 0 0 1 1
  0 0 1 1 1
  0 1 1 1 1];
ar = repmat(ar,[1 1 3]);
%ar(:,:,1) = min(ar(:,:,1),fcolor(1));
%ar(:,:,2) = min(ar(:,:,2),fcolor(2));
%ar(:,:,3) = min(ar(:,:,3),fcolor(3));

%%%% Slider to select iso-surface value to display
maxvalue = max(max(max(d3d)));
minvalue = min(min(min(d3d)));
isovalue  = 0.6*maxvalue;
uicontrol('Style', 'Slider', 'tag','isovalueslider',...
    'Value', isovalue, 'Min', minvalue, 'Max', maxvalue, ...
    'SliderStep', [0.01 0.01],...
    'Units', 'Normalized', 'Position', [0.1 0.1 0.8 0.075], ...
    'backgroundcolor', fcolor, ...
    'tooltip','Change Isosurface', ...
    'Callback', {@newisosurf, x,y,z,d3d} );...
uicontrol('style','text','units','normalized',...
    'string',num2str(round(minvalue)),'position',[0.03 0.1 0.05 0.075], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center','fontsize',10,'fontweight','bold');
uicontrol('style','text','units','normalized',...
    'string',num2str(round(maxvalue)),'position',[0.90 0.1 0.05 0.075],...
    'backgroundcolor',fcolor, ...
    'horizontalal','center','fontsize',10,'fontweight','bold');
uicontrol('style','text','units','normalized', ...
    'string','Iso-Surface Value' ,...
    'position',[0.05 0.01 0.3 0.075], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Style', 'edit', 'tag','isovaluebox',...
    'Units', 'Normalized', 'Position', [0.425, 0.01 0.15 0.09], ...
    'backgroundcolor', fcolor, ...
    'string', num2str(isovalue,'%6.2f'), ...
    'fontsize', 10, ...
    'tooltip','Change Isosurface', ...
    'Callback', {@newisosurf, x,y,z,d3d} );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Panel to control Viewer Position
vpos = uipanel('Title','Viewer Position', ...
    'FontSize',10, ...
    'Fontweight', 'bold', ...
    'BackgroundColor', fcolor,...
    'Position',[.01 .21 .975 .21]);
    [laz,lel] = view(findobj('tag','isoaxis')); 
    % Checkbox to switch on/off free 3D rotation
    uicontrol('Parent', vpos, 'style','checkbox','tag','rot3donoff', ...
        'units','normalized','pos',[0.78 0.4 0.22 0.5], ...
        'backgroundcolor',fcolor, ...
        'fontsize', 11, 'fontweight','bold', ...
        'string','Free Rotate','value',0,...
        'callback', @freerot); %...
    % Edit box for Viewer Azimuth 
    uicontrol('Parent', vpos, 'Style', 'edit', 'tag','isoaxesazbox',...
        'Units', 'Normalized', 'Position', [0.01 0.1 0.175 0.8], ...
        'backgroundcolor', fcolor, ...
        'fontsize', 10, ...
        'string', num2str(round(laz)), ...
        'tooltip','Enter Viewer Azimuth', ...
        'Callback', ...
       ['[laz,lel] = view(findobj(''tag'',''isoaxis'')); ' ...
        'laz = str2num(get(findobj(''tag'',''isoaxesazbox''),''string'')); ' ...
        'view(findobj(''tag'',''isoaxis''),laz,lel); ' ...
        'clear laz lel; ' ]); 
    uicontrol('Parent',vpos, 'style','text','units','normalized', ...
        'string','Azimuth' ,...
        'position',[0.2 0.1 0.15 0.8], ...
        'backgroundcolor',fcolor, ...
        'horizontalal','left',...
        'fontweight', 'bold', ...
        'fontsize',10); 
     % Edit box for Viewer Elevation 
    uicontrol('Parent', vpos, 'Style', 'edit', 'tag','isoaxeselevbox',...
        'Units', 'Normalized', 'Position', [0.35 0.1 0.175 0.8], ...
        'backgroundcolor', fcolor, ...
        'fontsize', 10, ...
        'string', num2str(round(lel)), ...
        'tooltip','Enter Viewer Elevation', ...
        'Callback', ...
       ['[laz,lel] = view(findobj(''tag'',''isoaxis'')); ' ...
        'lel = str2num(get(findobj(''tag'',''isoaxeselevbox''),''string'')); ' ...
        'view(findobj(''tag'',''isoaxis''),laz,lel); ' ...
        'clear laz lel; ' ]); 
    uicontrol('Parent',vpos, 'style','text','units','normalized', ...
        'string','Elevation' ,...
        'position',[0.54 0.1 0.15 0.8], ...
        'backgroundcolor',fcolor, ...
        'horizontalal','left',...
        'fontweight', 'bold', ...
        'fontsize',10);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Panel to control data aspect ratio
asp = uipanel('Title','Aspect Ratio', ...
    'FontSize',10, ...
    'Fontweight', 'bold', ...
    'BackgroundColor', fcolor,...
    'Position',[.01 .425 .2625 .575]);
    daspr = get(findobj('tag','isoaxis'),'dataaspectratio');
% Z-axis aspect 
    uicontrol('parent',asp, 'style','text', ...
        'string','Z:', 'fontsize',12, 'fontweight','bold', ...
        'backgroundcolor', fcolor, ...
        'position', [10, 10, 20, 20]);
    uicontrol('parent',asp, 'style','edit', 'tag', 'zaspectbox', ...
        'backgroundcolor', fcolor, ...
        'tooltip','Set z-axis aspect', ...
        'position',[35, 6, 60, 28], ...
        'string', num2str(round(daspr(3))), ...
        'fontsize',11, 'fontweight', 'bold', ...
        'Callback', ...
       ['daspr = get(findobj(''tag'',''isoaxis''),''dataaspectratio''); ' ...
        'daspr(3) = str2num(get(findobj(''tag'',''zaspectbox''),''string'')); ' ...
        'set(findobj(''tag'',''isoaxis''),''dataaspectratio'',daspr); ' ...
        'clear daspr; ']);
        % Fine tuning - down pointed arrow
        ar1 = permute(ar,[2 1 3]);
        uicontrol('parent',asp, 'style','push','pos',[96 5 20 15], ...
          'CData',ar1,  ...
          'tooltip','Tune z-axis aspect', ...
          'callback',...
           ['daspr = get(findobj(''tag'',''isoaxis''),''dataaspectratio''); ' ...
            'daspr(3) = str2num(get(findobj(''tag'',''zaspectbox''),''string'')) ' ...
            '- daspr(3)/20; ' ...
            'set(findobj(''tag'',''isoaxis''),''dataaspectratio'',daspr); ' ...
            'set(findobj(''tag'',''zaspectbox''),''string'',num2str(daspr(3),''%6.1f'')); ' ...
            'clear daspr; ']);
    % Fine tuning - up pointed arrow
    ar1 = permute(ar,[2 1 3]); ar1 = ar1(end:-1:1,:,:);
    uicontrol('parent',asp, 'style','push','pos',[96 20 20 15], ...
      'CData',ar1,  ...
      'tooltip','Tune z-axis aspect', ...
      'callback',...
       ['daspr = get(findobj(''tag'',''isoaxis''),''dataaspectratio''); ' ...
        'daspr(3) = str2num(get(findobj(''tag'',''zaspectbox''),''string'')) ' ...
        '+ daspr(3)/20; ' ...
        'set(findobj(''tag'',''isoaxis''),''dataaspectratio'',daspr); ' ...
        'set(findobj(''tag'',''zaspectbox''),''string'',num2str(daspr(3),''%6.1f'')); ' ...
        'clear daspr; ']);
% Y-axis aspect 
    uicontrol('parent',asp, 'style','text', ...
        'string','Y:', 'fontsize',12, 'fontweight','bold', ...
        'backgroundcolor', fcolor, ...
        'position', [10, 40, 20, 20]);
    uicontrol('parent',asp, 'style','edit', 'tag', 'yaspectbox', ...
        'backgroundcolor', fcolor, ...
        'tooltip','Set y-axis aspect', ...
        'position',[35, 36, 60, 28], ...
        'string', num2str(round(daspr(2))), ...
        'fontsize',11, 'fontweight', 'bold', ...
        'Callback', ...
       ['daspr = get(findobj(''tag'',''isoaxis''),''dataaspectratio''); ' ...
        'daspr(2) = str2num(get(findobj(''tag'',''yaspectbox''),''string'')); ' ...
        'set(findobj(''tag'',''isoaxis''),''dataaspectratio'',daspr); ' ...
        'clear daspr; ']);
    % Fine tuning - down pointed arrow
    ar1 = permute(ar,[2 1 3]);
    uicontrol('parent',asp, 'style','push','pos',[96 35 20 15], ...
      'CData',ar1,  ...
      'tooltip','Tune y-axis aspect', ...
       'callback',...
       ['daspr = get(findobj(''tag'',''isoaxis''),''dataaspectratio''); ' ...
        'daspr(2) = str2num(get(findobj(''tag'',''yaspectbox''),''string'')) ' ...
        '- daspr(2)/20; ' ...
        'set(findobj(''tag'',''isoaxis''),''dataaspectratio'',daspr); ' ...
        'set(findobj(''tag'',''yaspectbox''),''string'',num2str(daspr(2),''%7.2f'')); ' ...
        'clear daspr; ']);
    % Fine tuning - up pointed arrow
    ar1 = permute(ar,[2 1 3]); ar1 = ar1(end:-1:1,:,:);
    uicontrol('parent',asp, 'style','push','pos',[96 50 20 15], ...
      'CData',ar1,  ...
      'tooltip','Tune y-axis aspect', ...
      'callback',...
       ['daspr = get(findobj(''tag'',''isoaxis''),''dataaspectratio''); ' ...
        'daspr(2) = str2num(get(findobj(''tag'',''yaspectbox''),''string'')) ' ...
        '+ daspr(2)/20; ' ...
        'set(findobj(''tag'',''isoaxis''),''dataaspectratio'',daspr); ' ...
        'set(findobj(''tag'',''yaspectbox''),''string'',num2str(daspr(2),''%7.2f'')); ' ...
        'clear daspr; ']);
% X-axis aspect 
    uicontrol('parent',asp, 'style','text', ...
        'string','X:', 'fontsize',12, 'fontweight','bold', ...
        'backgroundcolor', fcolor, ...
        'position', [10, 70, 20, 20]);
    uicontrol('parent',asp, 'style','edit', 'tag', 'xaspectbox', ...
        'backgroundcolor', fcolor, ...
        'tooltip','Set x-axis aspect', ...
        'position',[35, 66, 60, 28], ...
        'string', num2str(round(daspr(1))), ...
        'fontsize',11, 'fontweight', 'bold', ...
        'Callback', ...
       ['daspr = get(findobj(''tag'',''isoaxis''),''dataaspectratio''); ' ...
        'daspr(1) = str2num(get(findobj(''tag'',''xaspectbox''),''string'')); ' ...
        'set(findobj(''tag'',''isoaxis''),''dataaspectratio'',daspr); ' ...
        'clear daspr; ']);
    % Fine tuning - down pointed arrow
    ar1 = permute(ar,[2 1 3]);
    uicontrol('parent',asp, 'style','push','pos',[96 65 20 15], ...
      'CData',ar1, ...
      'tooltip','Tune x-axis aspect', ...
      'callback',...
       ['daspr = get(findobj(''tag'',''isoaxis''),''dataaspectratio''); ' ...
        'daspr(1) = str2num(get(findobj(''tag'',''xaspectbox''),''string'')) ' ...
        '- daspr(1)/20; ' ...
        'set(findobj(''tag'',''isoaxis''),''dataaspectratio'',daspr); ' ...
        'set(findobj(''tag'',''xaspectbox''),''string'',num2str(daspr(1),''%7.2f'')); ' ...
        'clear daspr; ']);
    % Fine tuning - up pointed arrow
    ar1 = permute(ar,[2 1 3]); ar1 = ar1(end:-1:1,:,:);
    uicontrol('parent',asp, 'style','push','pos',[96 80 20 15], ...
      'CData',ar1, ...
      'tooltip','Tune x-axis aspect', ...
       'callback',...
      ['daspr = get(findobj(''tag'',''isoaxis''),''dataaspectratio''); ' ...
        'daspr(1) = str2num(get(findobj(''tag'',''xaspectbox''),''string'')) ' ...
        '+ daspr(1)/20; ' ...
        'set(findobj(''tag'',''isoaxis''),''dataaspectratio'',daspr); ' ...
        'set(findobj(''tag'',''xaspectbox''),''string'',num2str(daspr(1),''%7.2f'')); ' ...
        'clear daspr; ']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Panel to control position of light source
lht = uipanel('Title','Light Position', ...
    'FontSize',10, ...
    'Fontweight', 'bold', ...
    'BackgroundColor', fcolor,...
    'Position',[.28 .425 .51 .575]);
    [laz,lel]=lightangle(findobj('tag','isolight'));
    % Slider to control elevation
    uicontrol('Parent',lht,'Style', 'Slider', 'tag','lightelevhandle',...
        'Value', lel, 'Min', -90, 'Max', 90, ...
        'SliderStep', [0.01 0.01],...
        'Units', 'Normalized', 'Position', [0.028 0.13 0.6 0.18], ...
        'backgroundcolor', fcolor, ...
        'tooltip','Change elevation of light source', ...
        'Callback',...
       ['[laz,lel]=lightangle(findobj(''tag'',''isolight'')); ' ...
        'lel = get(findobj(''tag'',''lightelevhandle''),''value''); ' ...
        'lightangle(findobj(''tag'',''isolight''),laz,lel); ' ...
        'set(findobj(''tag'',''lightelevbox''),''string'', num2str(round(lel))); ' ...
        'clear laz lel; ' ]); 
    % Slider Annotation
    uicontrol('parent', lht, 'style','text',...
        'units','normalized','string','Elevation' ,...
        'position',[0.175 0.35 0.3 0.15], ...
        'backgroundcolor',fcolor, ...
        'horizontalal','center',...
        'fontsize',10, 'fontweight','bold');
   uicontrol('parent', lht, 'style','text',...
        'units','normalized','string','-90' ,...
        'position',[0.0 0.35 0.1 0.15], ...
        'backgroundcolor',fcolor, ...
        'horizontalal','center',...
        'fontsize',10, 'fontweight','bold');
    uicontrol('parent' ,lht, 'style','text',...
        'units','normalized', 'string','90', ...
        'position',[0.54 0.35 0.1 0.15],...
        'backgroundcolor',fcolor, ...
        'horizontalal','center','fontsize',10, 'fontweight','bold');
    % Edit box for Elevation 
    uicontrol('Parent', lht, 'Style', 'edit', 'tag','lightelevbox',...
        'Units', 'Normalized', 'Position', [0.7 0.09 0.25 0.25], ...
        'backgroundcolor', fcolor, ...
        'fontsize', 10, ...
        'string', num2str(lel,'%6.2f'), ...
        'tooltip','Enter Viewer Elevation', ...
        'Callback', ...
       ['[laz,lel]=lightangle(findobj(''tag'',''isolight'')); ' ...
        'lel = str2num(get(findobj(''tag'',''lightelevbox''),''string'')); ' ...
        'lightangle(findobj(''tag'',''isolight''),laz,lel); ' ...
        'set(findobj(''tag'',''lightelevhandle''),''value'',lel); ' ...
        'clear laz lel; ' ]); 
    
    % Slider to control azimuth
    uicontrol('Parent',lht,'Style', 'Slider', 'tag','lightazhandle',...
        'Value', laz, 'Min', 0, 'Max', 360, ...
        'SliderStep', [0.01 0.01],...
        'Units', 'Normalized', 'Position', [0.025 0.6 0.6 0.18], ...
        'backgroundcolor', fcolor, ...
        'tooltip','Change azimuth of light source', ...
        'Callback',...
       ['[laz,lel]=lightangle(findobj(''tag'',''isolight'')); ' ...
        'laz = get(findobj(''tag'',''lightazhandle''),''value''); ' ...
        'lightangle(findobj(''tag'',''isolight''),laz,lel); ' ...
        'set(findobj(''tag'',''lightazbox''),''string'', num2str(round(laz))); ' ...
        'clear laz lel; ' ]); 
    %Slider Annotation
    uicontrol('parent', lht, 'style','text',...
        'units','normalized','string','Azimuth' ,...
        'position',[0.175 0.82 0.3 0.15], ...
        'backgroundcolor',fcolor, ...
        'horizontalal','center',...
        'fontsize',10, 'fontweight','bold');
    uicontrol('parent', lht, 'style','text',...
        'units','normalized','string','0' ,...
        'position',[0.0 0.82 0.1 0.15], ...
        'backgroundcolor',fcolor, ...
        'horizontalal','center',...
        'fontsize',10, 'fontweight','bold');
    uicontrol('parent' ,lht, 'style','text',...
        'units','normalized', 'string','360', ...
        'position',[0.54 0.82 0.1 0.15],...
        'backgroundcolor',fcolor, ...
        'horizontalal','center','fontsize',10, 'fontweight','bold');
    % Edit box for Azimuth 
    uicontrol('Parent', lht, 'Style', 'edit', 'tag','lightazbox',...
        'Units', 'Normalized', 'Position', [0.7 0.56 0.25 0.25], ...
        'backgroundcolor', fcolor, ...
        'fontsize', 10, ...
        'string', num2str(laz,'%6.2f'), ...
        'tooltip','Enter azimuth of light source', ...
        'Callback', ...
       ['[laz,lel]=lightangle(findobj(''tag'',''isolight'')); ' ...
        'laz = str2num(get(findobj(''tag'',''lightazbox''),''string'')); ' ...
        'lightangle(findobj(''tag'',''isolight''),laz,lel); ' ...
        'set(findobj(''tag'',''lightazhandle''),''value'',laz); ' ...
        'clear laz lel; ' ]); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Push button to change surface colour
uicontrol('style','push', ...
    'units','normalized','pos',[0.8 0.6 0.15 0.15], ...
    'fontsize', 10, 'fontweight','bold', ...
    'string','ReColor','value',0,...
    'callback',...
   ['newcolor = uisetcolor(get(findobj(''tag'',''isosurface''),''facecolor'')); ' ...
    'set(findobj(''tag'',''isosurface''),''facecolor'',newcolor); ' ...
    'clear newcolor; ']); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Push button to Close Session
uicontrol('style','push', ...
    'units','normalized','pos',[0.8 0.8 0.15 0.15], ...
    'fontsize', 10, 'fontweight','bold', ...
    'string','Close','value',0,...
    'callback', 'delete(findobj(''tag'',''isofigure'')); '); 

% Wrapping up
%clear ar ar1 c daspr fcolor maxvalue minvalue isovalue 
%clear   
return

function newisosurf(obj,evd,x, y, z, d3d) %#ok<INUSL>
fc = get(findobj('tag','isosurface'),'facecolor'); 
[laz,lel]=lightangle(findobj('tag','isolight')); 
if strcmp(get(obj,'style'),'slider')
    isovalue = get(findobj('tag','isovalueslider'),'value'); 
elseif strcmpi(get(obj,'style'),'edit')
    isovalue = str2num(get(findobj('tag','isovaluebox'),'string'));
else
    return
end
figure(findobj('tag','isofigure')); 
[rdf,rdv] = isosurface(x,y,z,d3d,isovalue); 
delete(findobj('tag','isosurface')); 
delete(findobj('tag','isolight')); 
rds = patch('Faces',rdf,'Vertices',rdv); 
set(rds,'tag','isosurface','edgecolor','none','facecolor',fc); 
%delete(findobj('tag','isocaps')); 
%[cf,cv] = isocaps(x,y,z,d3d,isovalue,'above');
%rdc = patch('Faces',cf,'Vertices',cv);
%set(rdc,'tag','isocaps','edgecolor','none','facecolor',fc)
isonormals(x,y,z,d3d,rdf);
lh = camlight('headlight'); 
lighting gouraud; 
lightangle(lh,laz,lel); 
set(lh,'tag','isolight'); 
set(findobj('tag','isoaxis'),'zdir','reverse'); 
if strcmp(get(obj,'style'),'slider')
    set(findobj('tag','isovaluebox'),'string', num2str(isovalue)); 
elseif strcmpi(get(obj,'style'),'edit')
    set(findobj('tag','isovalueslider'),'value', isovalue); 
end
return

function freerot(obj,evd) %#ok<INUSD>
rot3d = get(findobj('tag','rot3donoff'),'value');
if rot3d == 1,
    %rotate3d(findobj('Tag','isofigure'),'on');
    figure(findobj('Tag','isofigure'));
    rh = rotate3d;
    set(rh, 'enable','on', 'ActionPostCallback',@setviewerpos);
elseif rot3d == 0,
    rotate3d(findobj('Tag','isofigure'),'off');
end;
return
function setviewerpos(obj,evd) %#ok<INUSL>
newView = round(get(evd.Axes,'View'));
set(findobj('tag','isoaxesazbox'),'string',num2str(newView(1),'%6.2f')); 
set(findobj('tag','isoaxeselevbox'),'string',num2str(newView(2),'%6.2f')); 
return


