function imagecolors()
%
% IMAGECOLORS: Change Image Display colors. The default colormap if JET
%              (rainbow). Each of the menu choices changes the colormap to:
%              GREY, BONE, JET, HSV, HOT, COOL, COPPER and Blue-white-red
%              color reel. These changes are GLOBAL and permanent. In
%              addition, there's a Blue-White-Red color reel scaled exactly
%              w.r.t. the maximum and minimum values of the displayed data
%              set and white will always correspond to zero. This change is
%              neither global nor permanent, but applies only to the
%              particular figure where it has been invoked and for as long
%              as it is alive. 
%              Colormap manipulation features include: 
%              COLOUR SATURATION : increases or decreases contrast
%                                  (saturation) via a slider GUI. This
%                                  option is available only in 'image'
%                                  display mode. The change applies only to
%                                  current figure.
%              BRIGHTEN          : increases the brightness.
%              DARKEN            : decreases the brightness.
%              FLIPUD            : reverses color order.
%              PERMUTE           : interchanges color indices
%              The changes apply to the current global colormap until they
%              are reset or the colormap changes. 
%
%   Requires : twocolors.m, twocolorsc.m
%
%     Author : Andreas Tzanis, 
%              Department of Geophysics, 
%              University of Athens
%              atzanis@geol.uoa.gr
%              (C) 2005, Andreas Tzanis, all rights reserved
%

global ENVAR
maps = {'gray','bone','jet', 'hot', 'Blue - Red', 'Brown - Black', ...
     'Red - Black'};
% Setup GUI with predefined colours
colormenuh = uimenu(gcf,'label','图像颜色');
% The first menu item sets the handle for the handle function callback  
uimenu(colormenuh,'tag', 'cmgui', 'label', maps{1}, ...
    'callback', {@cmset, maps{1}});
for k = 2:size(maps,2);
   uimenu(colormenuh,'label',maps{k}, 'callback',{@cmset, maps{k}});
end
% Adjust colour saturation
uimenu(colormenuh, 'Label','颜色饱和度', 'tag', 'csmenuitem', ...
    'Separator', 'on', ...
    'Callback', @colorsaturation);
if ~strcmp(ENVAR.DisplayOption,'image'),
    set(findobj('tag', 'csmenuitem'), 'Enable', 'off');
end
% Setup gui for colormap manipulaltion functions
uimenu(colormenuh,'label','提高亮度','Separator','on', ...
    'callback',{@cmadjust,'brighten'});
uimenu(colormenuh,'label','降低亮度', ...
    'callback',{@cmadjust,'darken'});
% uimenu(colormenuh,'label','flipud', ...
%     'callback',{@cmadjust,'flipud'});
% uimenu(colormenuh,'tag','cmadjust','label','permute',...
%     'callback',{@cmadjust,'permute'}); 
return

function cmset(cmgui,eventdata,colormaplabel) %#ok<INUSL>
global ENVAR
if strcmp(colormaplabel,'Blue - Red'),
    % sets the blue-white-red color reel    
    ENVAR.colormap = twocolors(64, [0 0 0.4], [0.4 0 0]); 
elseif strcmp(colormaplabel,'Brown - Black'),
    % sets the Brown-white-black color reel    
    ENVAR.colormap = twocolors(64, [0.5 0.25 0.0], [0 0 0]); 
elseif strcmp(colormaplabel,'Red - Black'),
    % sets the red-white-black color reel    
    ENVAR.colormap = twocolors(64, [0.4 0 0], [0 0 0]); 
else
%     ENVAR.colormap = colormap(colormaplabel(1:strfind(colormaplabel,' ')-1));
     ENVAR.colormap = colormap(colormaplabel);
end
% Set colormap for the "GPR DATA" figure
datafig = findobj('tag','datafigure');
if ~isempty(datafig),
    set(datafig,'colormap',ENVAR.colormap);
    colormap(ENVAR.colormap);
end
% Set colormap for the "PROCESSED GPR DATA" (output data) figure
procdatafig = findobj('tag','procdatafigure');
if ~isempty(procdatafig),
    set(procdatafig,'colormap',ENVAR.colormap);
end
% Set colormap for the "Fit Diffraction Hyperbola" figure
velfig = findobj('tag','velanfigure');
if ~isempty(velfig),
    set(velfig,'colormap',ENVAR.colormap);
end
% Set colormap for the "F-K Filter" figure
fkfig = findobj('tag','fkfiltfig');
if ~isempty(fkfig),
    set(fkfig,'colormap',ENVAR.colormap);
end
% Set colormap for the "TAU-P Filter" figure
fkfig = findobj('tag','taupfiltfig');
if ~isempty(fkfig),
    set(fkfig,'colormap',ENVAR.colormap);
end
return

function colorsaturation(obj,evt) %#ok<INUSD>
% Adjusts color saturation in the current axes
cd = get(findobj(gcf,'type','image'),'cdata');
% if there are more axes in one figure, use the extire value range!
if iscell(cd),
    cx = [min(min(min(cd{1:max(size(cd))}))) max(max(max(cd{1:max(size(cd))})))];
else
    cx = [min(min(cd)) max(max(cd))];
end
maxcolor = max(abs(cx));
cvalue   = get(findobj(gcf,'type','axes'),'clim');
if iscell(cvalue),
    cvalue = max(max(cvalue{1:max(size(cvalue))}));
else
    cvalue = max(get(gca,'clim'));
end
if cvalue > maxcolor,
    cvalue = maxcolor;
end
fcolor = get(gcf,'color');
clear cd;
% Hide uicontrols - make sure this jobs finishes before anything else is
% done
hh = findobj(gcf,'type','uicontrol'); 
set(hh,'visible','off'); 
%set up uipanel 
cs = uipanel('Title','', 'tag', 'datacspanel', ...
    'FontSize',10, ...
    'Fontweight', 'bold', ...
    'BackgroundColor', fcolor,...
    'Position',[.92 .1 .06 0.825]);
uicontrol('Parent',cs, 'Style', 'Slider', 'tag','datacshandle',...
    'Value', cvalue, 'Min', 1e-6*maxcolor, 'Max', maxcolor, ...
    'SliderStep', [0.01 0.1],...
    'Units', 'Normalized', 'Position', [0.2 0.2 0.6 0.75], ...
    'tooltip','Change Colour Saturation', ...
    'Callback',...
   ['cvalue = get(gco,''value''); ' ...
    'maxcolor = get(gco,''Max''); ' ...
    'cx = get(gca,''userdata''); ' ...
    'if isempty(cx), ' ...
    '    cd = get(findobj(gcf,''type'',''image''),''cdata''); ' ...
    '    if iscell(cd), ' ...
    '       cx = [min(min(min(cd{1:max(size(cd))}))) ' ...
    '             max(max(max(cd{1:max(size(cd))})))]; ' ...
    '    else; ' ...
    '       cx = [min(min(cd)) max(max(cd))]; ' ...
    '    end; ' ...
    'end; ' ...
    'cax = cx*cvalue/maxcolor; ' ...
    'hh = findobj(gcf,''type'',''axes''); ' ...
    'set(hh,''Clim'', cax); ' ...
    'clear cvalue maxcolor cax transl cx cd hh; ' ]);
uicontrol('Parent',cs, 'Style', 'push', 'tag','datacsok',...
    'string', 'OK', ...
    'Units', 'Normalized', 'Position', [0.2 0.05 0.6 0.1], ...
    'tooltip','Change Colour Saturation', ...
    'Callback', ...
   ['hh = findobj(gcf,''type'',''uicontrol''); ' ...
    'set(hh,''visible'',''on''); ' ...
    'delete(findobj(''tag'',''datacspanel'')); ' ...
    'clear hh; ']);
return

function cmadjust(cmadjust,eventdata,action) %#ok<INUSL>
% adjusts colormap (brightens, darkens, flips, permutes)
if strcmp(action,'brighten'),
    brighten(.25);
end
if strcmp(action,'darken'),
    brighten(-.25);
end
if strcmp(action,'flipud'),
    colormap(flipud(colormap));
end
if strcmp(action,'permute'),
    c = colormap;
    colormap(c(:,[2 3 1]));
end
return