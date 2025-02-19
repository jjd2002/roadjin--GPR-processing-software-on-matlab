function welcome(infotext)
%
% WELCOME : The MATGPR information service. 
%           1. On startup, creates the MATGPR information window on the
%              MATGPR GUI, window and displays the GPL statement and
%              disclaimer stored in the attached subfunction "gpl_statement".
%           2. In runtime is projects on the MATGPR information window the
%              data passed in the (string) cell array "infotext". 
%              ** "Infotext" is normally generated and passed by "showinfo.m" 
%
%    Usage : welcome(infotext)
%
%    Input : infotext – A column cell array of strings. 
%
%   Author : Andreas Tzanis,
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

mess = figure(findobj('tag','fi0'));          % get the MATGPR home window 
c1 = 0.5;    c2 = 0.5;    c3 = 0.5;           % and paint it gray 

if nargin == 0                                % initialize service
    datatitlestring ='Data Info';
%    rand('seed',sum(100*clock));
%        whitebg(mess,[rand rand rand]); 
        whitebg(mess,[c1 c2 c3]); 
    clf
    set(gca,'visible','off','drawmode','fast');
    text1 = text(0.45,1.05,'Welcome to MATGPR R2') ;
    set(text1,'FontSize',13,'Color','k','FontWeight','bold',...
        'horizontalal','center')
% prepare to display info on input data ("GPR DATA")
    top=0.9;
    left=0.05;
    right=0.95;
    bottom=0.05;
    labelheight=0.07;
    spacing=0.02;
% Draw the text window frame
    frameBorder=0.02;
    framePosition=[left-frameBorder bottom-frameBorder ...
            (right-left)+2*frameBorder (top-bottom)+2*frameBorder];
    uicontrol( 'Style','frame', 'Tag', 'indataframe', ...
        'Units','normalized', 'Position',framePosition, ...
        'BackgroundColor',[0.0 0.5 0.5]);
% Draw the text label
    labelPosition=[left top-labelheight (right-left) labelheight];
    uicontrol( 'Style','text', 'Tag', 'indatatitstr', ...
        'Units','normalized', 'Position',labelPosition, ...
        'BackgroundColor',[0.0 0.5 0.5], 'ForegroundColor',[1 1 1], ...
        'String',datatitlestring, 'fontsize', 10, 'fontweight', 'demi');
% Display the info box and display the info text
    textPosition=[left bottom (right-left) top-bottom-labelheight-2*spacing];
% Display the GPL license information
    licinfotext = gpl_statement;
    uicontrol( 'Style','edit', 'tag', 'indatainfobox', ...
        'Units','normalized', 'BackgroundColor',[1 1 1], 'Max',18, ...
        'String',licinfotext, 'horizontalal', 'center', 'fontsize', 10, ...
        'Position',textPosition);
    return
end % if nargin == 0 
% Show info on input data
if nargin == 1,     
    fontn = get(0,'FixedWidthFontName');
    do = ['set(findobj(''tag'',''indatainfobox''),''String'',' ...
        'infotext, ''horizontalal'',''left'',''fontname'',fontn);' ];
    er = 'welcome'; 
    eval(do,er);                % figure(mess)
    
end
return

function licinfo = gpl_statement()
%%% Text of the GPL licence statement
lincinfo = cell(1,1);
licinfo(1,1) = cellstr('Copyright (C) 2005, 2008, 2010  Andreas Tzanis');
licinfo(2,1) = cellstr('    ');
licinfo(3,1) = cellstr('This program is free software; you can redistribute it and/or modify');
licinfo(4,1) = cellstr('it under the terms of the GNU General Public License as published by');
licinfo(5,1) = cellstr('the Free Software Foundation; either version 2 of the License, or');
licinfo(6,1) = cellstr('any later version.');
licinfo(7,1) = cellstr('    ');
licinfo(8,1) = cellstr('This program is distributed in the hope that it will be useful,');
licinfo(9,1) = cellstr('but WITHOUT ANY WARRANTY; without even the implied warranty of');
licinfo(10,1) = cellstr('MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the');
licinfo(11,1) = cellstr('GNU General Public License for more details.');
licinfo(12,1) = cellstr('    ');
licinfo(13,1) = cellstr('You should have received a copy of the GNU General Public License');
licinfo(14,1) = cellstr('along with this program; if not, write to the Free Software');
licinfo(15,1) = cellstr('Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.');
return