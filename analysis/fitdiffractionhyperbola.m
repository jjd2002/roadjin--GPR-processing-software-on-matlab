function fitdiffractionhyperbola(DATA)
%
% FITDIFRACTIONHYPERBOLA : Interactively fit a difraction front hyperbola
% to common-offset GPR data. Simple approximation to the problem of fitting
% a difraction front hyperbola, assuming non-dispersive propagation in a
% uniform halfspace and point diffractors, or finite sized targets with 
% circular cross section! Parameters involved in are the halfspace
% velocity, the radius of the target and the coordinates of its centre
% (depth and location in the scan line). The parameters can be changed
% interactively, in three ways: 
% 1. By means of slider UI controls. 
% 2. By pointing the the cursor at the apex of a difraction hyperbola and
%    clicking the MIDDLE mouse button to begin the fitting process, then
%    using the scroll-wheel of the mouse to automatically change velocity
%    and target depth until the hyperbola is fitted - this is a fast easy
%    and  accurate method, but only for MATLAB releases featuring the
%    scroll-wheel callback function ...
% 3. Using any combination of (1) and (2) above.
% The fitting is done with the eye and thus depends on the experience of
% the user. The function returns the halfspace velocity as a global runtime
% variable. 
%
%    Usage : fitdiffractionhyperbola(DATA)
%
%    Input : DATA is the MATGPR data structure (IPD or OPD)
%
%  Returns : vhalfspace - the velocity of the halfspace as a global
%                         variable. 
% Requires : hyperbola.m
%
%  Author  : A. Tzanis
%            Department of Geophysics,
%            University of Athens.
%            atzanis@geol.uoa.gr
%
% Copyright (C) 2005, 2008, Andreas Tzanis, all rights reserved
%

% Initialize global variables
global ENVAR VS vhalfspace 

% setup figure
velfig = findobj('tag','velanfigure');
if ~isempty(velfig),
    figure(velfig);
    clf;
else
    figure('name','速度拟合',...
        'tag','velanfigure', ...
        'numbertitle','off', ...
        'menubar','none', ...
        'position',[80 150 800 500], ...
        'CreateFcn',['figuretools; imagecolors; ' ...
        'matgprwindows(''setup''); ' ...
        'matgprwindows(''updateadd'');'],...
        'DeleteFcn',['matgprwindows(''updateremove'', ' ...
        'findobj(''tag'',''velanfigure''));']);
end

% set up axes
axes('position',[0.1300    0.2100    0.7750    0.7150]);
imagesc(DATA.x,DATA.tt2w,DATA.d);
xlabel(DATA.xlab);  ylabel(DATA.zlab);
colormap(ENVAR.colormap);
%imagecolors
set(gca,'nextplot','add');

% set up parameters for diffraction hyperbola calculations
TxRx = DATA.TxRx;
hstart = min(DATA.x(1),DATA.x(length(DATA.x)));
hstop  = max(DATA.x(1),DATA.x(length(DATA.x)));
target_hloc   = hstart + (hstop-hstart)/2.;
target_depth  = 1.0;
target_radius = 0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(VS.v1d) && size(VS.v1d,1) == 1,
    vhalfspace = VS.v1d(1,1);
else
    vhalfspace    = 0.1;               % this is the default start-up value
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hloc,ttime] = hyperbola(hstart,hstop,TxRx,target_hloc,target_depth,...
    target_radius,vhalfspace);
hyp_handle = plot(hloc,ttime,'k','EraseMode','xor');
set(hyp_handle,'tag','hyperbolahandle');

% Velocity slider
vmin = 0.01;
vmax = 0.2998;
uicontrol('style','text','units','normalized','string',num2str(vmin,3),...
    'position',[0.01 0.085 0.05 0.03],'horizontalal','center',...
    'fontsize',8);
uicontrol('style','text','units','normalized','string',num2str(vmax,3),...
    'position',[0.01+0.18 0.085 0.05 0.03],'horizontalal','center',...
    'fontsize',8);
uicontrol('style','text','units','normalized','string','速度 (m/ns)',...
    'position',[0.01+0.04 0.085 0.15 0.03],'horizontalal','center',...
    'fontsize',9);
uicontrol('style','text','tag','vt4','units','normalized',...
    'string',num2str(vhalfspace,3),'position',[0.07 0.0205 0.1 0.03],...
    'horizontalal','center','fontsize',9);
uicontrol('Style', 'Slider', 'Tag', 'vel_handle', ...
    'Value', vhalfspace, 'Min', vmin, 'Max', vmax, ...
    'SliderStep', [0.001 0.01],'Units', 'Normalized', ...
    'Position', [0.01 0.055 0.23 0.03], ...
    'tooltip','Use this slider to control the velocity of the halfspace', ...
    'Callback',{@hyp1,hstart,hstop,TxRx});

% Depth slider
dmin  = 0;
dmax  = 0.1*(DATA.ns -1)*DATA.dt;
uicontrol('style','text','units','normalized','string',num2str(dmin,3),...
    'position',[0.25 0.085 0.05 0.03],'horizontalal','center',...
    'fontsize',8);
uicontrol('style','text','units','normalized','string',num2str(dmax,3),...
    'position',[0.25+0.18 0.085 0.05 0.03],'horizontalal','center',...
    'fontsize',8);
uicontrol('style','text','units','normalized','string','深度 (m)' ,...
    'position',[0.25+0.04 0.085 0.15 0.03],'horizontalal','center',...
    'fontsize',9);
uicontrol('style','text','tag','dt4','units','normalized',...
    'string',num2str(target_depth,3),'position',[0.31 0.0205 0.1 0.03],...
    'horizontalal','center','fontsize',9);
uicontrol('Style', 'Slider', 'Tag', 'depth_handle', ...
    'Value', target_depth, 'Min', dmin, 'Max', dmax, ...
    'SliderStep', [0.001 0.01], 'Units', 'Normalized',...
    'Position', [0.25 0.055 0.23 0.03],  ...
    'tooltip','Use this slider to control the depth of the target', ...
    'Callback', {@hyp1,hstart,hstop,TxRx});

% Target radius slider
rmin  = 0;
rmax  = dmax/4.;
uicontrol('style','text','units','normalized','string',num2str(rmin,3),...
    'position',[0.49 0.085 0.05 0.03],'horizontalal','center',...
    'fontsize',8);
uicontrol('style','text','units','normalized','string',num2str(rmax,3),...
    'position',[0.49+0.18 0.085 0.05 0.03],'horizontalal','center',...
    'fontsize',8);
uicontrol('style','text','units','normalized','string','半径 (m)' ,...
    'position',[0.49+0.04 0.085 0.15 0.03],'horizontalal','center',...
    'fontsize',9);
uicontrol('style','text','tag','rt4','units','normalized',...
    'string',num2str(target_radius,3),'position',[0.56 0.0205 0.1 0.03],...
    'horizontalal','center','fontsize',9);
uicontrol('Style', 'Slider', 'Tag', 'radius_handle', ...
    'Value', target_radius, 'Min', rmin, 'Max', rmax, ...
    'SliderStep', [0.001 0.01],'Units', 'Normalized', ...
    'Position', [0.49 0.055 0.23 0.03],  ...
    'tooltip','Use this slider to control the radius of the target', ...
    'Callback', {@hyp1,hstart,hstop,TxRx});

% Target location slider
uicontrol('style','text','units','normalized','string',num2str(hstart,3),...
    'position',[0.73 0.085 0.05 0.03],'horizontalal','center',...
    'fontsize',8);
uicontrol('style','text','units','normalized','string',num2str(hstop,3),...
    'position',[0.73+0.18 0.085 0.05 0.03],'horizontalal','center',...
    'fontsize',8);
uicontrol('style','text','units','normalized','string','位置 (m)' ,...
    'position',[0.73+0.04 0.085 0.15 0.03],'horizontalal','center',...
    'fontsize',9);
uicontrol('style','text','tag','ht4','units','normalized',...
    'string',num2str(target_hloc,3),'position',[0.80 0.0205 0.1 0.03],...
    'horizontalal','center','fontsize',9);
uicontrol('Style', 'Slider', 'Tag' , 'hloc_handle', ...
    'Value', target_hloc, 'Min', hstart, 'Max', hstop, ...
    'SliderStep', [0.001 0.001], 'Units', 'Normalized', ...
    'Position', [0.73 0.055 0.23 0.03],  ...
    'tooltip','Use this slider to control the location of the target', ...
    'Callback', {@hyp1,hstart,hstop,TxRx});

% Set up the mouse scroll-wheel callback function - easy find velocity and
% target depth ... But only for MATLAB releases featuring a scroll-wheel
% callback function ...
if isprop(findobj('tag','velanfigure'),'WindowScrollWheelFcn'),
    set(findobj('tag','velanfigure'),...
        'WindowScrollWheelFcn', {@hyp2,hstart,hstop,TxRx,'scroll'});
end
% Set up the window button down callback function - just point the cursor
% at the apex of a difraction hyperbola and click MIDDLE mouse button to
% begin the fitting process ...
set(findobj('tag','velanfigure'),...
        'WindowButtonDownFcn', {@hyp2,hstart,hstop,TxRx,'set'});

% Done button - exit and keep last velocity value
uicontrol('style','pushbutton','string','完成', ...
    'Units', 'Normalized', 'Position', [0.91 0.6 0.08 0.08], ...
    'tooltip','Exit and accept the fitted parameters', 'Callback', ...
    'delete(findobj(''tag'',''velanfigure'')); ');
%    ['vhalfspace = get(findobj(''tag'',''vel_handle''),''value''); ' ...
%     'delete(findobj(''tag'',''velanfigure'')); ']);

% Calcel button - exit and erase velocity value
uicontrol('style','pushbutton','string','取消',...
    'Units', 'Normalized', 'Position', [0.91 0.4 0.08 0.08], ...
    'tooltip','Exit and discard all parameters','Callback', @CancelFcn );

return

function hyp1(obj,evt, hstart, hstop, TxRx) %#ok<INUSL>
% Function-handle callback to calulate diffraction hyperbola with given 
% parameters of the medium and target
vhalfspace = get(findobj('tag','vel_handle'),'value'); 
target_depth = get(findobj('tag','depth_handle'),'value'); 
target_radius = get(findobj('tag','radius_handle'),'value'); 
target_hloc = get(findobj('tag','hloc_handle'),'value'); 
%%%
[hloc,ttime] = hyperbola(hstart,hstop,TxRx,target_hloc,target_depth,...
    target_radius,vhalfspace); 
%%%
set(findobj('tag','vt4'),'string',num2str(vhalfspace,3)); 
set(findobj('tag','dt4'),'string',num2str(target_depth,3));  
set(findobj('tag','rt4'),'string',num2str(target_radius,3)); 
set(findobj('tag','ht4'),'string',num2str(target_hloc,3)); 
set(findobj('tag','hyperbolahandle'),'ydata',ttime);
%%%
assignin('caller','vhalfspace',vhalfspace);
%%%
time_loc = (target_depth - target_radius)*2/vhalfspace;
cp = [target_hloc time_loc];
set(findobj('tag','velanfigure'),'userdata',cp);
return

function hyp2(obj, evt, hstart, hstop, TxRx, action) %#ok<INUSL>
vhalfspace = get(findobj('tag','vel_handle'),'value');
target_radius = get(findobj('tag','radius_handle'),'value');
% Point and click middle mouse button at the apex of a diffraction
% hyperbola (or any other point on the radargram image) - get point
% coordinates and store them for later use in the figure's "userdata"
% buffer.
if strcmp(action,'set'),
    if strcmp(get(gcf,'selectiontype'),'extend'),
        cp = get(gca,'currentpoint');
        target_hloc = cp(1,1);
        time_loc   = cp(1,2);
        set(findobj('tag','velanfigure'),'userdata',cp);
    else
        return
    end
    target_depth = target_radius + time_loc*vhalfspace/2;
    [hloc,ttime] = hyperbola(hstart,hstop,TxRx,target_hloc,target_depth,...
        target_radius,vhalfspace);
end
% At a given location, fit the hyperbola adjusting the velocity and the
% DEPTH TO THE CENTRE of the target using the mouse scroll-wheel.
if strcmp(action,'scroll'),
    target_depth = get(findobj('tag','depth_handle'),'value');
    cp = get(findobj('tag','velanfigure'),'userdata');
    if ~isempty(cp)
        target_hloc = cp(1,1);
        time_loc   = cp(1,2);
    else
        target_hloc = get(findobj('tag','hloc_handle'),'value');
        time_loc = (target_depth - target_radius)*2/vhalfspace;
    end
% Here's were we actually use the scroll wheel...    
    if evt.VerticalScrollCount > 0
        vhalfspace = vhalfspace + 0.0025;
        target_depth = target_radius + time_loc*vhalfspace/2;
    elseif evt.VerticalScrollCount < 0
        vhalfspace = vhalfspace - 0.0025;
        target_depth = target_radius + time_loc*vhalfspace/2;
    end
    [hloc,ttime] = hyperbola(hstart,hstop,TxRx,target_hloc,target_depth,...
        target_radius,vhalfspace);
    target_depth = target_radius + time_loc*vhalfspace/2;
end
% Update all UI controls ...
set(findobj('tag','vt4'),'string',num2str(vhalfspace,4));
set(findobj('tag','vel_handle'),'value',vhalfspace);
set(findobj('tag','dt4'),'string',num2str(target_depth,3));
set(findobj('tag','depth_handle'),'value',target_depth);
% set(findobj('tag','rt4'),'string',num2str(target_radius,3));
% set(findobj('tag','radius_handle'),'value',target_radius);
set(findobj('tag','ht4'),'string',num2str(target_hloc,3));
set(findobj('tag','hloc_handle'),'value',target_hloc);
set(findobj('tag','hyperbolahandle'),'ydata',ttime);
assignin('caller','vhalfspace',vhalfspace);
return

function CancelFcn(obj,evt)
vhalfspace = [];
assignin('caller','vhalfspace',vhalfspace);
delete(findobj('tag','velanfigure')); 
return

