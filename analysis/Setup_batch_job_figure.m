function Setup_batch_job_figure
%
% SETUP_BATCH_JOB_FIGURE :
% Prepares the GUI to be used in setting up the batch job queue
% The GUI must be saved in the "analysis" directory under the name
% "Batch_Job_Figure.fig", whence it is called by "Setup_batch_job.m"
%
%  Author : Andreas Tzanis, 
%           Department of Geophysics, 
%           University of Athens, 
%           atzanis@geol.uoa.gr
%
% Copyright (C) 2010, Andreas Tzanis. All rights reserved.
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
figure('name','Set Up Queue for Batch Job', ...`
    'tag','BatchJobSetupFig', ...
    'position', [530   120   660   770], ...
    'numbertitle','off', ...
    'Resize', 'off', ...
    'Color',[0.76 0.86 0.69]);         % 'Color',[0 0.5 0.5]);
set(gcf, 'MenuBar','none')
axis off;
fcolor = [0.76 0.87 0.78];             % fcolor = get(gcf,'color');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust Signal Position
sigpos = uipanel('Title','调整信号坐标', ...
    'FontSize',10, ...
    'Fontweight', 'bold', ...
    'BackgroundColor', fcolor,...
    'Position',[0.01 0.875 .45 .125]);
uicontrol('parent', sigpos, 'style','text','units','normalized', ...
    'string','步骤 #' ,...
    'position',[0.03 0.75 0.3 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', sigpos, 'Style', 'edit', 'tag','BatchStep#',...
    'Units', 'Normalized', 'Position', [0.07 0.3 0.2 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip', 'Priority of this Processing Step in the Batch Que', ...
    'Callback', ...
    ['if isempty(get(gco,''string'')), ' ...
    '    set(findobj(''tag'',''BatchSigpos1''),''string'',''''); ' ...
    '    set(findobj(''tag'',''BatchSigpos2''),''string'',''''); ' ...
    'end; ']);
% Give Signal Position in ns
uicontrol('parent', sigpos, 'style','text','units','normalized', ...
    'string','ns' ,...
    'position',[0.32 0.75 0.3 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', sigpos, 'Style', 'edit', 'tag','BatchSigpos1',...
    'Units', 'Normalized', 'Position', [0.35 0.3 0.2 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip','Give Signal Position in ns', ...
    'Callback', 'set(findobj(''tag'',''BatchSigpos2''),''string'','''')');
uicontrol('parent', sigpos, 'style','text','units','normalized', ...
    'string','或' ,...
    'position',[0.575 0.475 0.1 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
% Give Signal Position in Number of Samples
uicontrol('parent', sigpos, 'style','text','units','normalized', ...
    'string','# 样例' ,...
    'position',[0.65 0.75 0.3 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', sigpos, 'Style', 'edit', 'tag','BatchSigpos2',...
    'Units', 'Normalized', 'Position', [0.7 0.3 0.2 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip','Give Signal Position in Number of Samples', ...
    'Callback', 'set(findobj(''tag'',''BatchSigpos1''),''string'','''')')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trim Time Window
trimtw = uipanel('Title','修剪时窗', ...
    'FontSize',10, ...
    'Fontweight', 'bold', ...
    'BackgroundColor', fcolor,...
    'Position',[0.51 0.875 .45 .125]);
uicontrol('parent', trimtw, 'style','text','units','normalized', ...
    'string','步骤 #' ,...
    'position',[0.03 0.75 0.3 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', trimtw, 'Style', 'edit', 'tag','BatchStep#',...
    'Units', 'Normalized', 'Position', [0.07 0.3 0.2 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip','Priority of this Processing Step in the Batch Que', ...
    'Callback', ...
    ['if isempty(get(gco,''string'')), ' ...
    '    set(findobj(''tag'',''TrimTW1''),''string'',''''); ' ...
    '    set(findobj(''tag'',''TrimTW2''),''string'',''''); ' ...
    'end; ']);
% Give Cutoff Time in ns
uicontrol('parent', trimtw, 'style','text','units','normalized', ...
    'string','ns' ,...
    'position',[0.32 0.75 0.3 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', trimtw, 'Style', 'edit', 'tag','TrimTW1',...
    'Units', 'Normalized', 'Position', [0.35 0.3 0.2 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip','Give cutoff time in ns', ...
    'Callback', 'set(findobj(''tag'',''TrimTW2''),''string'','''')')
uicontrol('parent', trimtw, 'style','text','units','normalized', ...
    'string','或' ,...
    'position',[0.575 0.475 0.1 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
% Give Cutoff Time in Number of Samples
uicontrol('parent', trimtw, 'style','text','units','normalized', ...
    'string','# 样例' ,...
    'position',[0.65 0.75 0.3 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', trimtw, 'Style', 'edit', 'tag','TrimTW2',...
    'Units', 'Normalized', 'Position', [0.7 0.3 0.2 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip','Give cutoff time in Number of Samples', ...
    'Callback', 'set(findobj(''tag'',''TrimTW1''),''string'','''')')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove Global BackGround
rmgbgr = uipanel('Title','去除背景', ...
    'FontSize',10, ...
    'Fontweight', 'bold', ...
    'BackgroundColor', fcolor,...
    'Position',[0.01 0.74 .25 .125]);
uicontrol('parent', rmgbgr, 'style','text','units','normalized', ...
    'string','步骤 #' ,...
    'position',[0.2 0.7 0.6 0.3], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', rmgbgr, 'Style', 'edit', 'tag','BatchStep#',...
    'Units', 'Normalized', 'Position', [0.2 0.3 0.6 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip','Priority of this Processing Step in the Batch Que');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove DC component from each trace
rmdc = uipanel('Title','去除直流分量', ...
    'FontSize',10, ...
    'Fontweight', 'bold', ...
    'BackgroundColor', fcolor,...
    'Position',[0.285 0.74 .175 .125]);
uicontrol('parent', rmdc, 'style','text','units','normalized', ...
    'string','步骤 #' ,...
    'position',[0.2 0.7 0.6 0.3], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', rmdc, 'Style', 'edit', 'tag','BatchStep#',...
    'Units', 'Normalized', 'Position', [0.2 0.3 0.6 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip','Priority of this Processing Step in the Batch Que');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Dewow
% dwow = uipanel('Title','Dewow', ...
%     'FontSize',10, ...
%     'Fontweight', 'bold', ...
%     'BackgroundColor', fcolor,...
%     'Position',[0.51 0.74 .175 .125]);
% uicontrol('parent', dwow, 'style','text','units','normalized', ...
%     'string','STEP #' ,...
%     'position',[0.2 0.7 0.6 0.3], ...
%     'backgroundcolor',fcolor, ...
%     'horizontalal','center',...
%     'fontweight', 'bold', ...
%     'fontsize',10);
% uicontrol('Parent', dwow, 'Style', 'edit', 'tag','BatchStep#',...
%     'Units', 'Normalized', 'Position', [0.2 0.3 0.6 0.4], ...
%     'backgroundcolor', fcolor, ...
%     'fontsize', 10, ...
%     'string', '', ...
%     'tooltip','Priority of this Processing Step in the Batch Que');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Remove DZT Header Gain
% rmdzthg = uipanel('Title','DZT Header Gain', ...
%     'FontSize',10, ...
%     'Fontweight', 'bold', ...
%     'BackgroundColor', fcolor,...
%     'Position',[0.71 .74 .25 .125]);
% uicontrol('parent', rmdzthg, 'style','text','units','normalized', ...
%     'string','STEP #' ,...
%     'position',[0.2 0.7 0.6 0.3], ...
%     'backgroundcolor',fcolor, ...
%     'horizontalal','center',...
%     'fontweight', 'bold', ...
%     'fontsize',10);
% uicontrol('Parent', rmdzthg, 'Style', 'edit', 'tag','BatchStep#',...
%     'Units', 'Normalized', 'Position', [0.2 0.3 0.6 0.4], ...
%     'backgroundcolor', fcolor, ...
%     'fontsize', 10, ...
%     'string', '', ...
%     'tooltip','Priority of this Processing Step in the Batch Que');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard AGC
sagc = uipanel('Title','标准 AGC 增益', ...
    'FontSize',10, ...
    'Fontweight', 'bold', ...
    'BackgroundColor', fcolor,...
    'Position',[0.0925 0.605 .275 .125]);
uicontrol('parent', sagc, 'style','text','units','normalized', ...
    'string','步骤 #' ,...
    'position',[0.07 0.75 0.35 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', sagc, 'Style', 'edit', 'tag','BatchStep#',...
    'Units', 'Normalized', 'Position', [0.07 0.3 0.35 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip', 'Priority of this Processing Step in the Batch Que', ...
    'Callback', ...
    ['if isempty(get(gco,''string'')), ' ...
    '    set(findobj(''tag'',''StdAGC1''),''string'',''''); ' ...
    'end; ']);
% Length of AGC window in ns
uicontrol('parent', sagc, 'style','text','units','normalized', ...
    'string','时窗 (ns)' ,...
    'position',[0.49 0.75 0.5 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', sagc, 'Style', 'edit', 'tag','StdAGC1',...
    'Units', 'Normalized', 'Position', [0.57 0.3 0.35 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip','Define Length of AGC Window in ns');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian Tapered AGC
gtagc = uipanel('Title','高斯锥形 AGC 增益', ...
    'FontSize',10, ...
    'Fontweight', 'bold', ...
    'BackgroundColor', fcolor,...
    'Position',[0.51 0.605 .45 .125]);
uicontrol('parent', gtagc, 'style','text','units','normalized', ...
    'string','步骤 #' ,...
    'position',[0.03 0.75 0.3 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', gtagc, 'Style', 'edit', 'tag','BatchStep#',...
    'Units', 'Normalized', 'Position', [0.07 0.3 0.2 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip','Priority of this Processing Step in the Batch Que', ...
    'Callback', ...
    ['if isempty(get(gco,''string'')), ' ...
    '    set(findobj(''tag'',''GT-AGC1''),''string'',''''); ' ...
    '    set(findobj(''tag'',''GT-AGC2''),''string'',''''); ' ...
    'end; ']);
% Noise Level
uicontrol('parent', gtagc, 'style','text','units','normalized', ...
    'string','噪声等级' ,...
    'position',[0.34 0.75 0.3 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', gtagc, 'Style', 'edit', 'tag','GT-AGC1',...
    'Units', 'Normalized', 'Position', [0.39 0.3 0.2 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip','Define noise level (e.g. 5.0e-7)')
% Length of AGC window in ns
uicontrol('parent', gtagc, 'style','text','units','normalized', ...
    'string','时窗 (ns)' ,...
    'position',[0.63 0.75 0.35 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', gtagc, 'Style', 'edit', 'tag','GT-AGC2',...
    'Units', 'Normalized', 'Position', [0.7 0.3 0.2 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip','Define Length of AGC Window in ns')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse Power Decay
invpd = uipanel('Title','反向功率衰减', ...
    'FontSize',10, ...
    'Fontweight', 'bold', ...
    'BackgroundColor', fcolor,...
    'Position',[0.0925 0.47 .275 .125]);
uicontrol('parent', invpd, 'style','text','units','normalized', ...
    'string','步骤 #' ,...
    'position',[0.07 0.75 0.35 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', invpd, 'Style', 'edit', 'tag','BatchStep#',...
    'Units', 'Normalized', 'Position', [0.07 0.3 0.35 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip', 'Priority of this Processing Step in the Batch Que', ...
    'Callback', ...
    ['if isempty(get(gco,''string'')), ' ...
    '    set(findobj(''tag'',''Invpdecay1''),''string'',''''); ' ...
    'end; ']);
% Exponent of Gain Function g(t) = s*t^p
uicontrol('parent', invpd, 'style','text','units','normalized', ...
    'string','指数' ,...
    'position',[0.49 0.75 0.5 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', invpd, 'Style', 'edit', 'tag','Invpdecay1',...
    'Units', 'Normalized', 'Position', [0.57 0.3 0.35 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip','Exponent of Gain Function g(t) = s*t^p');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Inverse Amplitude Decay
% invad = uipanel('Title','Inverse Amplitude Decay', ...
%     'FontSize',10, ...
%     'Fontweight', 'bold', ...
%     'BackgroundColor', fcolor,...
%     'Position',[0.51 0.47 .45 .125]);
% uicontrol('parent', invad, 'style','text','units','normalized', ...
%     'string','STEP #' ,...
%     'position',[0.03 0.75 0.3 0.2], ...
%     'backgroundcolor',fcolor, ...
%     'horizontalal','center',...
%     'fontweight', 'bold', ...
%     'fontsize',10);
% uicontrol('Parent', invad, 'Style', 'edit', 'tag','BatchStep#',...
%     'Units', 'Normalized', 'Position', [0.07 0.3 0.2 0.4], ...
%     'backgroundcolor', fcolor, ...
%     'fontsize', 10, ...
%     'string', '', ...
%     'tooltip','Priority of this Processing Step in the Batch Que', ...
%     'Callback', ...
%     ['if isempty(get(gco,''string'')), ' ...
%     '    set(findobj(''tag'',''InvAD_Median''),''value'',1); ' ...
%     '    set(findobj(''tag'',''InvAD_Mean''),''value'',0); ' ...
%     '    set(findobj(''tag'',''InvAD_N1''),''string'',''''); ' ...
%     'end; ']);
% % Attenuation is modelled on the mean or median trace
% uicontrol('Parent', invad, 'style','checkbox','tag', 'InvAD_Mean', ...
%     'units','normalized', 'position', [0.34 0.60 0.3 0.3], ...
%     'backgroundcolor', fcolor, ...
%     'string', 'MEAN', ...
%     'fontsize', 10, 'fontweight','bold', ...
%     'value',0, ...
%     'callback', ...
%     ['set(findobj(''tag'',''InvAD_Median''),''value'',0); '...
%     'if get(findobj(''tag'',''InvAD_Median''),''value'') == 0 && '...
%     '   get(findobj(''tag'',''InvAD_Mean''),''value'') == 0, '...
%     '    set(findobj(''tag'',''InvAD_Median''),''value'',1); '...
%     'end;' ]);
% uicontrol('Parent', invad, 'style','checkbox','tag', 'InvAD_Median', ...
%     'units','normalized', 'position', [0.34 0.15 0.3 0.3], ...
%     'backgroundcolor', fcolor, ...
%     'string', 'MEDIAN', ...
%     'fontsize', 10, 'fontweight','bold', ...
%     'value',1, ...
%     'callback', ...
%     ['set(findobj(''tag'',''InvAD_Mean''),''value'',0); '...
%     'if get(findobj(''tag'',''InvAD_Median''),''value'') == 0 && '...
%     '   get(findobj(''tag'',''InvAD_Mean''),''value'') == 0, '...
%     '    set(findobj(''tag'',''InvAD_Mean''),''value'',1); '...
%     'end;' ]);
% % Order of exponential spectrum
% uicontrol('parent', invad, 'style','text','units','normalized', ...
%     'string','ORDER' ,...
%     'position',[0.63 0.75 0.35 0.2], ...
%     'backgroundcolor',fcolor, ...
%     'horizontalal','center',...
%     'fontweight', 'bold', ...
%     'fontsize',10);
% uicontrol('Parent', invad, 'Style', 'edit', 'tag','InvAD_N1',...
%     'Units', 'Normalized', 'Position', [0.7 0.3 0.2 0.4], ...
%     'backgroundcolor', fcolor, ...
%     'fontsize', 10, ...
%     'string', '', ...
%     'tooltip','Define Order of Exponential Spectrum')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Resample Time Axis
% rstime = uipanel('Title','Resample Time Axis', ...
%     'FontSize',10, ...
%     'Fontweight', 'bold', ...
%     'BackgroundColor', fcolor,...
%     'Position',[0.01 0.34 .275 .125]);
% uicontrol('parent', rstime, 'style','text','units','normalized', ...
%     'string','STEP #' ,...
%     'position',[0.07 0.75 0.35 0.2], ...
%     'backgroundcolor',fcolor, ...
%     'horizontalal','center',...
%     'fontweight', 'bold', ...
%     'fontsize',10);
% uicontrol('Parent', rstime, 'Style', 'edit', 'tag','BatchStep#',...
%     'Units', 'Normalized', 'Position', [0.07 0.3 0.35 0.4], ...
%     'backgroundcolor', fcolor, ...
%     'fontsize', 10, ...
%     'string', '', ...
%     'tooltip', 'Priority of this Processing Step in the Batch Que', ...
%     'Callback', ...
%     ['if isempty(get(gco,''string'')), ' ...
%     '    set(findobj(''tag'',''RsTime1''),''string'',''''); ' ...
%     'end; ']);
% % Output samples per scan
% uicontrol('parent', rstime, 'style','text','units','normalized', ...
%     'string','# Samples' ,...
%     'position',[0.49 0.75 0.5 0.2], ...
%     'backgroundcolor',fcolor, ...
%     'horizontalal','center',...
%     'fontweight', 'bold', ...
%     'fontsize',10);
% uicontrol('Parent', rstime, 'Style', 'edit', 'tag','RsTime1',...
%     'Units', 'Normalized', 'Position', [0.57 0.3 0.35 0.4], ...
%     'backgroundcolor', fcolor, ...
%     'fontsize', 10, ...
%     'string', '', ...
%     'tooltip','Give Number of Output Samples per Scan');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Resample Scan Axis
% rsscan = uipanel('Title','Resample Scan Axis', ...
%     'FontSize',10, ...
%     'Fontweight', 'bold', ...
%     'BackgroundColor', fcolor,...
%     'Position',[0.346 0.34 .275 .125]);
% uicontrol('parent', rsscan, 'style','text','units','normalized', ...
%     'string','STEP #' ,...
%     'position',[0.07 0.75 0.35 0.2], ...
%     'backgroundcolor',fcolor, ...
%     'horizontalal','center',...
%     'fontweight', 'bold', ...
%     'fontsize',10);
% uicontrol('Parent', rsscan, 'Style', 'edit', 'tag','BatchStep#',...
%     'Units', 'Normalized', 'Position', [0.07 0.3 0.35 0.4], ...
%     'backgroundcolor', fcolor, ...
%     'fontsize', 10, ...
%     'string', '', ...
%     'tooltip', 'Priority of this Processing Step in the Batch Que', ...
%     'Callback', ...
%     ['if isempty(get(gco,''string'')), ' ...
%     '    set(findobj(''tag'',''RsScan1''),''string'',''''); ' ...
%     'end; ']);
% % Output samples per scan
% uicontrol('parent', rsscan, 'style','text','units','normalized', ...
%     'string','# Traces' ,...
%     'position',[0.49 0.75 0.5 0.2], ...
%     'backgroundcolor',fcolor, ...
%     'horizontalal','center',...
%     'fontweight', 'bold', ...
%     'fontsize',10);
% uicontrol('Parent', rsscan, 'Style', 'edit', 'tag','RsScan1',...
%     'Units', 'Normalized', 'Position', [0.57 0.3 0.35 0.4], ...
%     'backgroundcolor', fcolor, ...
%     'fontsize', 10, ...
%     'string', '', ...
%     'tooltip','Give Number of Output Traces');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Equalize
% eqlize = uipanel('Title','Equalize', ...
%     'FontSize',10, ...
%     'Fontweight', 'bold', ...
%     'BackgroundColor', fcolor,...
%     'Position',[0.685 0.34 .275 .125]);
% uicontrol('parent', eqlize, 'style','text','units','normalized', ...
%     'string','STEP #' ,...
%     'position',[0.07 0.75 0.35 0.2], ...
%     'backgroundcolor',fcolor, ...
%     'horizontalal','center',...
%     'fontweight', 'bold', ...
%     'fontsize',10);
% uicontrol('Parent', eqlize, 'Style', 'edit', 'tag','BatchStep#',...
%     'Units', 'Normalized', 'Position', [0.07 0.3 0.35 0.4], ...
%     'backgroundcolor', fcolor, ...
%     'fontsize', 10, ...
%     'string', '', ...
%     'tooltip', 'Priority of this Processing Step in the Batch Que', ...
%     'Callback', ...
%     ['if isempty(get(gco,''string'')), ' ...
%     '    set(findobj(''tag'',''EqLize1''),''string'',''''); ' ...
%     'end; ']);
% % Output samples per scan
% uicontrol('parent', eqlize, 'style','text','units','normalized', ...
%     'string','Base Value' ,...
%     'position',[0.49 0.75 0.5 0.2], ...
%     'backgroundcolor',fcolor, ...
%     'horizontalal','center',...
%     'fontweight', 'bold', ...
%     'fontsize',10);
% uicontrol('Parent', eqlize, 'Style', 'edit', 'tag','EqLize1',...
%     'Units', 'Normalized', 'Position', [0.57 0.3 0.35 0.4], ...
%     'backgroundcolor', fcolor, ...
%     'fontsize', 10, ...
%     'string', '', ...
%     'tooltip','Give Number of Output Traces');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency Filters
% Frequency Filters
firff = uipanel('Title','FIR 频率滤波器', ...
    'FontSize',10, ...
    'Fontweight', 'bold', ...
    'BackgroundColor', fcolor,...
    'Position',[0.01 0.21 .95 .125]);
uicontrol('parent', firff, 'style','text','units','normalized', ...
    'string','步骤 #' ,...
    'position',[0.015 0.75 0.2 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', firff, 'Style', 'edit', 'tag','BatchStep#',...
    'Units', 'Normalized', 'Position', [0.04 0.3 0.125 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip','Priority of this Processing Step in the Batch Que', ...
    'Callback', ...
    ['if isempty(get(gco,''string'')), ' ...
    '    set(findobj(''tag'',''FIRFF_LP''),''value'',1); ' ...
    '    set(findobj(''tag'',''FIRFF_HP''),''value'',0); ' ...
    '    set(findobj(''tag'',''FIRFF_BP''),''value'',0); ' ...
    '    set(findobj(''tag'',''FIRFF_BS''),''value'',0); ' ...
    '    set(findobj(''tag'',''FIRFF_F1''),''string'',''''); ' ...
    '    set(findobj(''tag'',''FIRFF_F2''),''string'',''''); ' ...
    '    set(findobj(''tag'',''FIRFF_F2''),''enable'',''off''); ' ...
    'end; ']);
% Filter Type
uicontrol('Parent', firff, 'style','checkbox','tag', 'FIRFF_LP', ...
    'units','normalized', 'position', [0.2 0.60 0.3 0.3], ...
    'backgroundcolor', fcolor, ...
    'string', '低通', ...
    'fontsize', 10, 'fontweight','bold', ...
    'value',1, ...
    'callback', ...
    ['set(findobj(''tag'',''FIRFF_F2''),''enable'',''off''); '...
    'set(findobj(''tag'',''FIRFF_HP''),''value'',0); '...
    'set(findobj(''tag'',''FIRFF_BP''),''value'',0); '...
    'set(findobj(''tag'',''FIRFF_BS''),''value'',0); ' ...
    'if get(findobj(''tag'',''FIRFF_LP''),''value'') == 0 && '...
    '   get(findobj(''tag'',''FIRFF_HP''),''value'') == 0 && '...
    '   get(findobj(''tag'',''FIRFF_BP''),''value'') == 0 && '...
    '   get(findobj(''tag'',''FIRFF_BS''),''value'') == 0, '...
    '    set(findobj(''tag'',''FIRFF_LP''),''value'',1); '...
    '    set(findobj(''tag'',''FIRFF_F2''),''enable'',''off''); '...
    'end;' ]);
uicontrol('Parent', firff, 'style','checkbox','tag', 'FIRFF_HP', ...
    'units','normalized', 'position', [0.2 0.15 0.3 0.3], ...
    'backgroundcolor', fcolor, ...
    'string', '高通', ...
    'fontsize', 10, 'fontweight','bold', ...
    'value',0, ...
    'callback', ...
    ['set(findobj(''tag'',''FIRFF_F2''),''enable'',''off''); '...
    'set(findobj(''tag'',''FIRFF_LP''),''value'',0); '...
    'set(findobj(''tag'',''FIRFF_BP''),''value'',0); '...
    'set(findobj(''tag'',''FIRFF_BS''),''value'',0); ' ...
    'if get(findobj(''tag'',''FIRFF_LP''),''value'') == 0 && '...
    '   get(findobj(''tag'',''FIRFF_HP''),''value'') == 0 && '...
    '   get(findobj(''tag'',''FIRFF_BP''),''value'') == 0 && '...
    '   get(findobj(''tag'',''FIRFF_BS''),''value'') == 0, '...
    '    set(findobj(''tag'',''FIRFF_LP''),''value'',1); '...
    '    set(findobj(''tag'',''FIRFF_F2''),''enable'',''off''); '...
    'end;' ]);
uicontrol('Parent', firff, 'style','checkbox','tag', 'FIRFF_BP', ...
    'units','normalized', 'position', [0.4 0.60 0.3 0.3], ...
    'backgroundcolor', fcolor, ...
    'string', '带通', ...
    'fontsize', 10, 'fontweight','bold', ...
    'value',0, ...
    'callback', ...
    ['set(findobj(''tag'',''FIRFF_F2''),''enable'',''on''); '...
    'set(findobj(''tag'',''FIRFF_LP''),''value'',0); '...
    'set(findobj(''tag'',''FIRFF_HP''),''value'',0); '...
    'set(findobj(''tag'',''FIRFF_BS''),''value'',0); ' ...
    'if get(findobj(''tag'',''FIRFF_LP''),''value'') == 0 && '...
    '   get(findobj(''tag'',''FIRFF_HP''),''value'') == 0 && '...
    '   get(findobj(''tag'',''FIRFF_BP''),''value'') == 0 && '...
    '   get(findobj(''tag'',''FIRFF_BS''),''value'') == 0, '...
    '    set(findobj(''tag'',''FIRFF_LP''),''value'',1); '...
    '    set(findobj(''tag'',''FIRFF_F2''),''enable'',''off''); '...
    'end;' ]);
uicontrol('Parent', firff, 'style','checkbox','tag', 'FIRFF_BS', ...
    'units','normalized', 'position', [0.4 0.15 0.3 0.3], ...
    'backgroundcolor', fcolor, ...
    'string', '带阻', ...
    'fontsize', 10, 'fontweight','bold', ...
    'value',0, ...
    'callback', ...
    ['set(findobj(''tag'',''FIRFF_F2''),''enable'',''on''); '...
    'set(findobj(''tag'',''FIRFF_LP''),''value'',0); '...
    'set(findobj(''tag'',''FIRFF_HP''),''value'',0); '...
    'set(findobj(''tag'',''FIRFF_BP''),''value'',0); '...
    'if get(findobj(''tag'',''FIRFF_LP''),''value'') == 0 && '...
    '   get(findobj(''tag'',''FIRFF_HP''),''value'') == 0 && '...
    '   get(findobj(''tag'',''FIRFF_BP''),''value'') == 0 && '...
    '   get(findobj(''tag'',''FIRFF_BS''),''value'') == 0, '...
    '    set(findobj(''tag'',''FIRFF_LP''),''value'',1); '...
    '    set(findobj(''tag'',''FIRFF_F2''),''enable'',''off''); '...
    'end;' ]);
% Define Cutoff Frequency for Low/High Pass Filters, or Lower Cutoff for
% Band Pass/Stop Filters
uicontrol('parent', firff, 'style','text','units','normalized', ...
    'string','F1' ,...
    'position',[0.6 0.75 0.15 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', firff, 'Style', 'edit', 'tag','FIRFF_F1',...
    'Units', 'Normalized', 'Position', [0.6 0.3 0.15 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'tooltip','Cutoff Frequency for Low/High Pass Filters, or Lower Cutoff for Band Pass/Stop Filters')
% Define Upper Cutoff Frequency for Band Pass/Stop Filters
uicontrol('parent', firff, 'style','text','units','normalized', ...
    'string','F2' ,...
    'position',[0.8 0.75 0.15 0.2], ...
    'backgroundcolor',fcolor, ...
    'horizontalal','center',...
    'fontweight', 'bold', ...
    'fontsize',10);
uicontrol('Parent', firff, 'Style', 'edit', 'tag','FIRFF_F2',...
    'Units', 'Normalized', 'Position', [0.8 0.3 0.15 0.4], ...
    'backgroundcolor', fcolor, ...
    'fontsize', 10, ...
    'string', '', ...
    'enable', 'off', ...
    'tooltip','Upper Cutoff Frequency for Band Pass/Stop Filters')

%  Proceed ...
uicontrol('Style','Pushbutton', 'Units','normalized', ...
    'Position',[.2 .05 .15 .05 ], 'String','设置', ...
    'backgroundcolor', [0.6 1 0], ...
    'Fontweight', 'bold', 'Fontsize',12, ...
    'Callback', ...
   ['Setup_batch_job(''setup'');' ...
    'delete(findobj(''tag'',''BatchJobSetupFig'')); ']);
%  Changed my mind ...
uicontrol('Style','Pushbutton','Units','normalized', ...
    'Position',[.65 .05 .15 .05 ],  'String','关闭', ...
    'backgroundcolor', [1 0.2 0], ...
    'Fontweight', 'bold', 'Fontsize',12, ...
    'Callback', 'delete(findobj(''tag'',''BatchJobSetupFig'')); ');

%%% Finished setting up the figure - save it
disp('===> Save the figure in the "analysis" directory as "Batch_Job_Figure.fig"')
saveas(findobj('tag','BatchJobSetupFig'),'C:\Users\Administrator\Desktop\road\analysis/Batch_Job_Figure.fig')
return