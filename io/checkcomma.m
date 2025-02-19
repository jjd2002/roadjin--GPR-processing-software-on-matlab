function outstr = checkcomma(instr)
%
% CHECKCOMMA: Replaces commas in the input string with dots. When applied
% to strings or cell arrays of strings returned from MATLAB ui dialogs
% (e.g. inputdlg), and which contain comma delimited decimal numbers 
% (Continental European notation), changes them to dot delimited decimal
% numbers (Anglosaxon, therefore MATLAB notation), so as to avoid
% inadvertent errors and crashes
%
% Usage
% outstr = checkcomma(instr)
%
% Input 
% instr  : The input string or cell array of strings, with comma delimited
%          decimals

% Output 
% outstr : The output string or cell array of strings with dot delimited
%          decimals
%
% Author : Andreas Tzanis
%          Department of Geophysics, 
%          University of Athens
%          atzanis@geol.uoa.gr
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
outstr = char(instr);
for i=1:size(outstr,1);
    comma = findstr(outstr(i,:),',');
    if ~isempty(comma),
        for j=1:length(comma); 
            outstr(i,comma) = '.'; 
        end;
    end
end
return