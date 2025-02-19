function [d_out,t_out,ns_out] = trimtimewindow(d,t,dt,ns)
%
% TRIMTIMEWINDOW: Trims the late times of the time window 
%
%        USAGE  :[d_out,t_out,ns_out] = trimtimewindow(d,t,dt,ns)
%
%        Inputs :  d : The 2-D GPR data matrix
%                  t : The vector of time coordinates (2 way traveltime)
%                 dt : The sampling rate in ns 
%                 ns : The dimension of t, equal to the column dimension
%                      of "d". 

%       Outputs : d_out : The trimmed GPR data matrix
%                 t_out : The trimmed vector of time coordinates
%                 ns_out: The new dimension t_out and d_out
%
%      Requires : checkcomma.m
%
%      Author   : Andreas Tzanis,
%                 Dept. of Geophysics,   
%                 University of Athens
%                 atzanis@geol.uoa.gr
%                 (C) 2005, Andreas Tzanis. All rights reserved.
%
    
if isempty(d), 
    erh = errordlg('No data - the array is empty!','EDITDATA: ERROR');
     uiwait(erh);   clear erh; 
    d_out = []; 
    t_out = [];
    ns_out = [];
    return
end

% Default outputs are : NO ACTION
d_out = d; %#ok<NASGU>
t_out = t;
ns_out = ns;

%%%%%    Input mode
inmode = menu('请选择修剪方式','使用光标','手动输入',...
    '取消');
if inmode == 3,
    d_out = [];
    return
end
if inmode ==1,
%%% Get data figure
    datafig = findobj('tag','datafigure');
    if ~ishandle(datafig),
        erh = errordlg('No data figure to work with',...
            'TRIM TIME WINDOW: ERROR');
        uiwait(erh);
        d_out = [];
        return   
    end
%%% display the help message
    msg = msgbox('Please point and click at the desired time.', ...
        'TRIM TIME WINDOW INFO','help');
    msgpos  = get(msg,'position');
    set(msg,'pos',[20 60 msgpos(3) msgpos(4)]);
    pause(0.5)
    figure(datafig)
    %pause(0.5)
    %[xtr,yt]=ginput(1);
    [xtr,yt]=getpoint('', ' (ns)');   
    nfrom = floor(yt/dt);
    if exist('msg'), %#ok<EXIST>
            close(msg);
    end
        
elseif inmode == 2,
        
    clb(1) = cellstr('修剪至 (ns) :');
    answer = inputdlg(clb,'TRIM TIME WINDOW : Give new range',1 ); 
    if isempty(answer)
        d_out = []; 
        return
    end
    lb = checkcomma(answer);   
    from = str2num(lb(1,:));
    nfrom  = floor(from/dt);
end
d_out  = d(1:nfrom-1,:); 
t_out  = t(1:nfrom-1);
ns_out = size(d_out,1);

return
