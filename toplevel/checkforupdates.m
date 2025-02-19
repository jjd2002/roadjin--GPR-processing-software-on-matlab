function checkforupdates(response)
%
% CHECKFORUPDATES: Checks the MATGPR home page for updates/upgrades and
%                  bug fixes - informs the user.
%
%          Usage : checkforupdates()
%
%         Author : Andreas Tzanis,
%                  Department of Geophysics, 
%                  University of Athens
%                  atzanis@geol.uoa.gr
%
% Copyright (C) 2007, Andreas Tzanis. All rights reserved.
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
try
     new_version = urlread('http://users.uoa.gr/~atzanis/matgpr/version.txt');
     %new_version = urlread('file:///C:/Users/Erric/Desktop/WORK\GPR/MATGPR/Xversion.txt');
catch
% in case user cannot access our site
    new_version = ENVAR.CURRENT_VERSION;
end
if ~strcmp(ENVAR.CURRENT_VERSION,new_version)
    wwandw = which('web'); 
    if ~isempty(wwandw),                  % MATLAB Web browser is available
        msgtext = ['The new matGPR version ' new_version ' has been released!'];
        msg = questdlg(msgtext,'matGPR Updates', 'Learn More', 'Download', 'Cancel', 'Learn More');
        waitfor(msg);
        if strcmp(msg,'Learn More'),
            %web('file:///C:\Users\Erric/Desktop/WORK/GPR/MATGPR/Updates_and_bug_fixes.txt')
            web('http://users.uoa.gr/~atzanis/matgpr/Updates_and_bug_fixes.txt')
            newtext = 'Go to "Check for Updates" to Download!';
            newmsg = msgbox(newtext,'matGPR Updates','help');
            waitfor(newmsg);
        end
        if strcmp(msg, 'Download'),
            web('http://users.uoa.gr/~atzanis/matgpr/matgpr.html')
        end
    else                              % MATLAB Web browser is NOT available
        msgtext = {['New MATGPR version ' new_version ' has been released.'],...
            '    ',...
            'Please visit http://users.uoa.gr/~atzanis/matgpr/matgpr.html ',...
            'to download and check for details.                           '};
        msg = msgbox(msgtext,'Software Update','warn');
        waitfor(msg);
    end
else
    if response == 1,
        msg = msgbox('Sorry, software updates are not presently available',...
            'matGPR Updates','error');
        waitfor(msg);
    end
end
