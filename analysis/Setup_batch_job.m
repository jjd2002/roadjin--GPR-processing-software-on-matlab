function Setup_batch_job(action)
%
%  SETUP_BATCH_JOB : 
%  Initializes the "Batch Job GUI" with the batch job queue parameters
%  supplied by the ENVAR.QUEUE structure and allows modification/ resetting
%  of the queue (ENVAR.QUEUE) for a new execution of Run_batch_job.m 
%
%  Usage : Setup_batch_job(action)
%  Keyword values : action = 'initialize' will setup the GUI and initialize
%                            the queue. 
%                   action = 'setup' will reset the queue with new values
%                            taken form the GUI.
%
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
global ENVAR

%%% Setup the batch job GUI and initialize it with parameters taken from
%%% ENVAR.QUEUE
%%%
if strcmpi(action,'initialize')
    
    %%% Load the Batch Job GUI (prepared with Setup_batch_job_figure.m)
    openfig('Batch_Job_Figure.fig','reuse');
    
    % Check for the existence of QUEUE in ENVAR. It may not exist when
    % upgrading from Rev. 2.1.3.2 - if so, initialize it
    if ~isfield(ENVAR,'QUEUE'),
        ENVAR.QUEUE.sigpos.step = 0;
        ENVAR.QUEUE.trimtw.step = 0;
        ENVAR.QUEUE.rmbackgr.step = 0;
        ENVAR.QUEUE.rmdc.step = 0;
        ENVAR.QUEUE.dewow.step = 0;
        ENVAR.QUEUE.dzthgain.step = 0;
        ENVAR.QUEUE.stdagc.step = 0;
        ENVAR.QUEUE.gtagc.step = 0;
        ENVAR.QUEUE.invpdecay.step = 0;
        ENVAR.QUEUE.invadecay.step = 0;
        ENVAR.QUEUE.rstime.step = 0;
        ENVAR.QUEUE.rsscan.step = 0;
        ENVAR.QUEUE.equalize.step = 0;
        ENVAR.QUEUE.firff.step = 0;
        perm = [1 2 15 3:14];
        ENVAR = orderfields(ENVAR,perm);
    end

    %%% Initialize Signal Position panel
    if ENVAR.QUEUE.sigpos.step,
        set(findobj('parent',findobj('title','Adjust Signal Position'),...
            'tag','BatchStep#'),...
            'string',num2str(ENVAR.QUEUE.sigpos.step))
        if ~isempty(ENVAR.QUEUE.sigpos.new_t0),
            set(findobj('tag','BatchSigpos1'),'string',num2str(ENVAR.QUEUE.sigpos.new_t0));
        end
        if ~isempty(ENVAR.QUEUE.sigpos.new_zero),
            set(findobj('tag','BatchSigpos2'),'string',num2str(ENVAR.QUEUE.sigpos.new_zero));
        end
    end
    %%% Initialize Time window trimming panel
    if ENVAR.QUEUE.trimtw.step,
        set(findobj('parent',findobj('title','Trim Time Window'),...
            'style','edit','tag','BatchStep#'),...
            'string',num2str(ENVAR.QUEUE.trimtw.step))
        if ~isempty(ENVAR.QUEUE.trimtw.new_tn),
            set(findobj('tag','TrimTW1'),'string',num2str(ENVAR.QUEUE.trimtw.new_tn));
        end
        if ~isempty(ENVAR.QUEUE.trimtw.cutoff),
            set(findobj('tag','TrimTW2'),'string',num2str(ENVAR.QUEUE.trimtw.cutoff));
        end
    end
    %%% Initialize Background Removal panel
    if ENVAR.QUEUE.rmbackgr.step,
        set(findobj('parent',findobj('title','Remove Background'),...
            'style','edit','tag','BatchStep#'),...
            'string',num2str(ENVAR.QUEUE.rmbackgr.step))
    end
    %%% Initialize DC removal panel
    if ENVAR.QUEUE.rmdc.step,
        set(findobj('parent',findobj('title','Remove DC'),...
            'style','edit','tag','BatchStep#'),...
            'string',num2str(ENVAR.QUEUE.rmdc.step ))
    end
    %%% Initialize Dewow panel
    if ENVAR.QUEUE.dewow.step,
        set(findobj('parent',findobj('title','Dewow'),...
            'style','edit','tag','BatchStep#'),...
            'string',num2str(ENVAR.QUEUE.dewow.step))
    end
    %%% Initialize DZT header gain removal panel
    if ENVAR.QUEUE.dzthgain.step,
        set(findobj('parent',findobj('title','DZT Header Gain'),...
            'style','edit','tag','BatchStep#'),...
            'string',num2str(ENVAR.QUEUE.dzthgain.step))
    end
    %%% Initialize Standard AGC panel
    if ENVAR.QUEUE.stdagc.step,
        set(findobj('parent',findobj('title','Standard AGC'),...
            'style','edit','tag','BatchStep#'),...
            'string',num2str(ENVAR.QUEUE.stdagc.step))
        set(findobj('tag','StdAGC1'),'string',num2str(ENVAR.QUEUE.stdagc.wagc))
    end
    %%% Initialize Gaussian-tapered AGC panel
    if ENVAR.QUEUE.gtagc.step,
        set(findobj('parent',findobj('title','Gaussian-tapered AGC'),...
            'style','edit','tag','BatchStep#'),...
            'string',num2str(ENVAR.QUEUE.gtagc.step))
        set(findobj('tag','GT-AGC1'),'string',num2str(ENVAR.QUEUE.gtagc.EPS))
        set(findobj('tag','GT-AGC2'),'string',num2str(ENVAR.QUEUE.gtagc.wagc))
    end
    %%% Initialize Inverse Power Decay panel
    if ENVAR.QUEUE.invpdecay.step,
        set(findobj('parent',findobj('title','Inverse Power Decay'),...
            'style','edit','tag','BatchStep#'),...
            'string',num2str(ENVAR.QUEUE.invpdecay.step))
        set(findobj('tag','Invpdecay1'),'string',num2str(ENVAR.QUEUE.invpdecay.pow))
    end
    %%% Initialize Inverse Amplitude Decay panel
    if ENVAR.QUEUE.invadecay.step,
        set(findobj('parent',findobj('title','Inverse Amplitude Decay'),...
            'style','edit','tag','BatchStep#'),...
            'string',num2str(ENVAR.QUEUE.invadecay.step))
        if ENVAR.QUEUE.invadecay.model == 1,
            set(findobj('tag','InvAD_Median'),'value','max')
        end
        if ENVAR.QUEUE.invadecay.model == 2,
            set(findobj('tag','InvAD_Mean'),'value','max')
        end
        set(findobj('tag','InvAD_N1'),'string',num2str(ENVAR.QUEUE.invadecay.order))
    end
    %%% Initialize Time Resampling panel
    if ENVAR.QUEUE.rstime.step,
        set(findobj('parent',findobj('title','Resample Time Axis'),...
            'style','edit','tag','BatchStep#'),...
            'string',num2str(ENVAR.QUEUE.rstime.step))
        set(findobj('tag','RsTime1'),'string',num2str(ENVAR.QUEUE.rstime.new_ns))
    end
    %%% Initialize Spatial Resampling panel
    if ENVAR.QUEUE.rsscan.step,
        set(findobj('parent',findobj('title','Resample Scan Axis'),...
            'style','edit','tag','BatchStep#'),...
            'string',num2str(ENVAR.QUEUE.rsscan.step))
        set(findobj('tag','RsScan1'),'string',num2str(ENVAR.QUEUE.rsscan.new_ntr))
    end
    %%% Initialize Equalization panel
    if ENVAR.QUEUE.equalize.step,
        set(findobj('parent',findobj('title','Equalize'),...
            'style','edit','tag','BatchStep#'),...
            'string',num2str(ENVAR.QUEUE.equalize.step))
        set(findobj('tag','EqLize1'),'string',num2str(ENVAR.QUEUE.equalize.base))
    end
    %%% Initialize Frequeuency Filtering panel
    if ENVAR.QUEUE.firff.step,
        set(findobj('parent',findobj('title','FIR Frequeuency Filter'),...
            'style','edit','tag','BatchStep#'),...
            'string',num2str(ENVAR.QUEUE.firff.step))
        if ENVAR.QUEUE.firff.filtertype == 1,
            set(findobj('tag','FIRFF_LP'),'value',1)
            set(findobj('tag','FIRFF_HP'),'value',0)
            set(findobj('tag','FIRFF_BP'),'value',0)
            set(findobj('tag','FIRFF_BS'),'value',0)
            set(findobj('tag','FIRFF_F1'),'string',num2str(ENVAR.QUEUE.firff.f1))
        end
        if ENVAR.QUEUE.firff.filtertype == 2,
            set(findobj('tag','FIRFF_LP'),'value',0)
            set(findobj('tag','FIRFF_HP'),'value',1)
            set(findobj('tag','FIRFF_BP'),'value',0)
            set(findobj('tag','FIRFF_BS'),'value',0)
            set(findobj('tag','FIRFF_F1'),'string',num2str(ENVAR.QUEUE.firff.f1))
        end
        if ENVAR.QUEUE.firff.filtertype == 3,
            set(findobj('tag','FIRFF_LP'),'value',0)
            set(findobj('tag','FIRFF_HP'),'value',0)
            set(findobj('tag','FIRFF_BP'),'value',1)
            set(findobj('tag','FIRFF_BS'),'value',0)
            set(findobj('tag','FIRFF_F1'),'string',num2str(ENVAR.QUEUE.firff.f1))
            set(findobj('tag','FIRFF_F2'),'enable','on',...
                'string',num2str(ENVAR.QUEUE.firff.f2))
        end
        if ENVAR.QUEUE.firff.filtertype == 4,
            set(findobj('tag','FIRFF_LP'),'value',0)
            set(findobj('tag','FIRFF_HP'),'value',0)
            set(findobj('tag','FIRFF_BP'),'value',0)
            set(findobj('tag','FIRFF_BS'),'value',1)
            set(findobj('tag','FIRFF_F1'),'string',num2str(ENVAR.QUEUE.firff.f1))
            set(findobj('tag','FIRFF_F2'),'enable','on',...
                'string',num2str(ENVAR.QUEUE.firff.f2))
        end
    end
end

%%% Set up the new batch job queue ------------------------------------------
%%%
if strcmpi(action,'setup')
    
    %%% Set Signal Position priority and parameters
    istep = str2num(get(findobj('parent',...
        findobj('title','Adjust Signal Position'),'style',...
        'edit','tag','BatchStep#'),'string'));
    if ~isempty(istep),
        ENVAR.QUEUE.sigpos.step = istep;
        if ~isempty(get(findobj('tag','BatchSigpos1'),'string')),
            ENVAR.QUEUE.sigpos.new_t0  = str2num(get(findobj('tag','BatchSigpos1'),'string'));
            ENVAR.QUEUE.sigpos.new_zero   = [];
        end
        if ~isempty(get(findobj('tag','BatchSigpos2'),'string')),
            ENVAR.QUEUE.sigpos.new_t0  = [];
            ENVAR.QUEUE.sigpos.new_zero = str2num(get(findobj('tag','BatchSigpos2'),'string'));
        end
    else
        set(findobj('tag','BatchSigpos1'),'string','')
        set(findobj('tag','BatchSigpos2'),'string','')
        ENVAR.QUEUE.sigpos.step = 0;
        ENVAR.QUEUE.sigpos.new_t0  = [];
        ENVAR.QUEUE.sigpos.new_zero   = [];

    end
    %%% Set Time Window trimming priority and parameters
    istep = str2num(get(findobj('parent',...
        findobj('title','Trim Time Window'),'style',...
        'edit','tag','BatchStep#'),'string'));
    if ~isempty(istep),
        ENVAR.QUEUE.trimtw.step = istep;
        if ~isempty(get(findobj('tag','TrimTW1'),'string')),
            ENVAR.QUEUE.trimtw.new_tn  = str2num(get(findobj('tag','TrimTW1'),'string'));
            ENVAR.QUEUE.trimtw.cutoff = [];
        end
        if ~isempty(get(findobj('tag','TrimTW2'),'string')),
            ENVAR.QUEUE.trimtw.new_tn = [];
            ENVAR.QUEUE.trimtw.cutoff = str2num(get(findobj('tag','TrimTW2'),'string'));
        end
    else
        set(findobj('tag','TrimTW1'),'string','')
        set(findobj('tag','TrimTW2'),'string','')
        ENVAR.QUEUE.trimtw.step = 0;
        ENVAR.QUEUE.trimtw.new_tn = [];
        ENVAR.QUEUE.trimtw.cutoff = [];
    end
    %%% Set Background removal priority
    istep = str2num(get(findobj('parent',...
        findobj('title','Remove Background'),'style',...
        'edit','tag','BatchStep#'),'string'));
    if ~isempty(istep),
        ENVAR.QUEUE.rmbackgr.step = istep;
    else
        ENVAR.QUEUE.rmbackgr.step = 0;
    end
    %%% Set DC removal priority
    istep = str2num(get(findobj('parent',...
        findobj('title','Remove DC'),'style',...
        'edit','tag','BatchStep#'),'string'));
    if ~isempty(istep),
        ENVAR.QUEUE.rmdc.step = istep;
    else
        ENVAR.QUEUE.rmdc.step = 0;
    end
    %%% Set Dewow priority
    istep = str2num(get(findobj('parent',...
        findobj('title','Dewow'),'style',...
        'edit','tag','BatchStep#'),'string'));
    if ~isempty(istep),
        ENVAR.QUEUE.dewow.step = istep;
    else
        ENVAR.QUEUE.dewow.step = 0;
    end
    %%% Set DZT header gain removal priority
    istep = str2num(get(findobj('parent',...
        findobj('title','DZT Header Gain'),'style',...
        'edit','tag','BatchStep#'),'string'));
    if~isempty(istep),
        ENVAR.QUEUE.dzthgain.step = istep;
    else
        ENVAR.QUEUE.dzthgain.step = 0;
    end
    %%% Set Standard AGC priority and parameters
    istep = str2num(get(findobj('parent',...
        findobj('title','Standard AGC'),'style',...
        'edit','tag','BatchStep#'),'string'));
    if ~isempty(istep),
        ENVAR.QUEUE.stdagc.step = istep;
        ENVAR.QUEUE.stdagc.wagc = str2num(get(findobj('tag','StdAGC1'),'string'));
    else
        set(findobj('tag','StdAGC1'),'string','')
        ENVAR.QUEUE.stdagc.step = 0;
        ENVAR.QUEUE.stdagc.wagc = [];
    end
    %%% Set Gaussian-tapered AGC priority and parameters
    istep = str2num(get(findobj('parent',...
        findobj('title','Gaussian-tapered AGC'),'style',...
        'edit','tag','BatchStep#'),'string'));
    if ~isempty(istep),
        ENVAR.QUEUE.gtagc.step = istep;
        ENVAR.QUEUE.gtagc.EPS  = str2num(get(findobj('tag','GT-AGC1'),'string'));
        ENVAR.QUEUE.gtagc.wagc = str2num(get(findobj('tag','GT-AGC2'),'string'));
    else
        set(findobj('tag','GT-AGC1'),'string','')
        set(findobj('tag','GT-AGC2'),'string','')
        ENVAR.QUEUE.gtagc.step = 0;
        ENVAR.QUEUE.gtagc.EPS  = [];
        ENVAR.QUEUE.gtagc.wagc = [];
    end
    %%% Set Inv. Power Decay priority and parameters
    istep = str2num(get(findobj('parent',...
        findobj('title','Inverse Power Decay'),'style',...
        'edit','tag','BatchStep#'),'string'));
    if ~isempty(istep),
        ENVAR.QUEUE.invpdecay.step = istep;
        ENVAR.QUEUE.invpdecay.pow = str2num(get(findobj('tag','Invpdecay1'),'string'));
    else
        ENVAR.QUEUE.invpdecay.step = 0;
        ENVAR.QUEUE.invpdecay.pow = [];
    end
    %%% Set Inv. Amplitude Decay priority and parameters
    istep = str2num(get(findobj('parent',...
        findobj('title','Inverse Amplitude Decay'),'style',...
        'edit','tag','BatchStep#'),'string'));
    if ~isempty(istep),
        ENVAR.QUEUE.invadecay.step = istep;
        if get(findobj('tag','InvAD_Median'),'value') == 1,
            ENVAR.QUEUE.invadecay.model = 1;
        end
        if get(findobj('tag','InvAD_Mean'),'value') == 1,
            ENVAR.QUEUE.invadecay.model = 2;
        end
        ENVAR.QUEUE.invadecay.order = str2num(get(findobj('tag','InvAD_N1'),'string'));
    else
        ENVAR.QUEUE.invadecay.step = 0;
        set(findobj('tag','InvAD_Median'),'value',1)
        set(findobj('tag','InvAD_Mean'),'value',0)
        ENVAR.QUEUE.invadecay.order = [];
    end
    %%% Set tome resampling priority and parameters
    istep = str2num(get(findobj('parent',...
        findobj('title','Resample Time Axis'),'style',...
        'edit','tag','BatchStep#'),'string'));
    if ~isempty(istep),
        ENVAR.QUEUE.rstime.step = istep;
        ENVAR.QUEUE.rstime.new_ns = str2num(get(findobj('tag','RsTime1'),'string'));
    else
        ENVAR.QUEUE.rstime.step = 0;
        ENVAR.QUEUE.rstime.new_ns = [];
    end
    %%% Set spatial resampling priority and parameters
    istep = str2num(get(findobj('parent',...
        findobj('title','Resample Scan Axis'),'style',...
        'edit','tag','BatchStep#'),'string'));
    if ~isempty(istep),
        ENVAR.QUEUE.rsscan.step = istep;
        ENVAR.QUEUE.rsscan.new_ntr = str2num(get(findobj('tag','RsScan1'),'string'));
    else
        ENVAR.QUEUE.rsscan.step = 0;
        ENVAR.QUEUE.rsscan.new_ntr = [];
    end
    %%% Set equalization priority and parameters
    istep = str2num(get(findobj('parent',...
        findobj('title','Equalize'),'style',...
        'edit','tag','BatchStep#'),'string'));
    if ~isempty(istep),
        ENVAR.QUEUE.equalize.step = istep;
        ENVAR.QUEUE.equalize.base = str2num(get(findobj('tag','EqLize1'),'string'));
    else
        ENVAR.QUEUE.equalize.step = 0;
        ENVAR.QUEUE.equalize.base = [];
    end
    %%% Set frequency filtering priority and parameters
    istep = str2num(get(findobj('parent',...
        findobj('title','FIR Frequency Filter'),'style',...
        'edit','tag','BatchStep#'),'string'));
    if ~isempty(istep),
        ENVAR.QUEUE.firff.step = istep;
        if get(findobj('tag','FIRFF_LP'),'value') == 1,
            ENVAR.QUEUE.firff.filtertype = 1;
            ENVAR.QUEUE.firff.f1 = str2num(get(findobj('tag','FIRFF_F1'),'string'));
            ENVAR.QUEUE.firff.f2 = [];
        end
        if get(findobj('tag','FIRFF_HP'),'value') == 1,
            ENVAR.QUEUE.firff.filtertype = 2;
            ENVAR.QUEUE.firff.f1 = str2num(get(findobj('tag','FIRFF_F1'),'string'));
            ENVAR.QUEUE.firff.f2 = [];
        end
        if get(findobj('tag','FIRFF_BP'),'value') == 1,
            ENVAR.QUEUE.firff.filtertype = 3;
            ENVAR.QUEUE.firff.f1 = str2num(get(findobj('tag','FIRFF_F1'),'string'));
            ENVAR.QUEUE.firff.f2 = str2num(get(findobj('tag','FIRFF_F2'),'string'));
        end
        if get(findobj('tag','FIRFF_BS'),'value') == 1,
            ENVAR.QUEUE.firff.filtertype = 4;
            ENVAR.QUEUE.firff.f1 = str2num(get(findobj('tag','FIRFF_F1'),'string'));
            ENVAR.QUEUE.firff.f2 = str2num(get(findobj('tag','FIRFF_F2'),'string'));
        end
    else
        set(findobj('tag','FIRFF_LP'),'value',1);
        set(findobj('tag','FIRFF_HP'),'value',0);
        set(findobj('tag','FIRFF_BP'),'value',0);
        set(findobj('tag','FIRFF_BS'),'value',0);
        set(findobj('tag','FIRFF_F1'),'string','')
        set(findobj('tag','FIRFF_F2'),'string','')
        ENVAR.QUEUE.firff.step = 0;
        ENVAR.QUEUE.firff.filtertype = [];
        ENVAR.QUEUE.firff.f1 = [];
        ENVAR.QUEUE.firff.f2 = [];
    end
end

return