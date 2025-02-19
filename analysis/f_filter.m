function  [dff, Ws] = f_filter(scanaxis, d, dt, filtertype, design_option)
%
%    F_FILTER : Design and apply FIR frequency filters
%
%       Usage : [dff, Ws] = freq_filter(scanaxis, d, dt, filtertype, ...
%                               design_option) 
%
%      Inputs : 
%   scanaxis  : Vector of the space coordinates of traces along the
%               scanaxis (need not be equally spaced);
%          d  : 2-D array of distance vs time input data;  
%         dt  : Sampling interval; 
% filtertype  : flag determining the type of the filter
%               =1  Low-pass filter
%               =2  High-pass filter
%               =3  Band-pass filter
%               =4  Band-stop filter
%design_option: keyword determining how the filter will be designed.
%               = 'TESTTR' : Use a test trace. Point and click on some test
%                            trace, then define the cutoff frequecy(ies)
%                            using the spectrum of a test trace as a
%                            template, also by pointing and clicking. 
%               = 'MEANTR' : Compute the mean (average) trace, then define
%                            the cutoff frequency(ies) by pointing and
%                            clicking, using the spectrum of the mean trace
%                            as template. 
%               = 'TYPEIN' : Enter predefined filter characteristics
%                            through dialog boxes.
%
%     Outputs : 
%        dff  : 2-D array of filtered data
%         Ws  : The cutoff frequency (filtertype = 1 or 2), or, the pass /
%               stop bands (filtertype = 3 or 4).
%
%   Requires  :  fir_f1.m, fir_f2.m and checkcomma.m
%
%     Author  :  Andreas Tzanis,
%                Department of Geophysics, 
%                University of Athens
%                atzanis@geol.uoa.gr
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%

if filtertype < 1 || filtertype >4 || ischar(filtertype), 
    erh = errordlg('Wrong filter type. Aborting!','FREQ_FILTER : ERROR');
    uiwait(erh);
    dff = [];  Ws = [];
    return
end;
%%% ns and ntr is Samples_per_scan and No_traces respectively
[ns, ntr] = size(d);
nf = floor(0.75*ns);                      % filter length    
fn = 1/(2*dt*0.001);                      % Nyquist multiplied by 1000 to 
df = fn/(ns);                             % give frequency axis in fake MHz
f  = df:df:fn;

if strcmp(design_option,'TESTTR') || strcmp(design_option,'MEANTR'),
    
    if strcmp(design_option,'TESTTR'), 
%%%%%%%%   Get a test trace to design filter 
        datafig = findobj('tag','datafigure');
        if ~ishandle(datafig),
            erh = errordlg('No data figure for the filter design interface',...
                'F_FILTER : ERROR');
            uiwait(erh);
            dff = [];  Ws = [];
            return
         end
%%%% display the help message
         message     = msgbox('Please click on a test trace',...
             'F_FILTER : HELP');
         messagepos  = get(message,'position');
         set(message,'pos',[20 60 150 messagepos(4)]);
         pause(0.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         figure(datafig)
         pause(0.5)
         [xtr,ytr]=getpoint(' (m)', ' (ns)');        %#ok<NASGU> 
         trnum = round((xtr-scanaxis(1))/(scanaxis(2)-scanaxis(1)));
         tr = d(:,trnum);
         tr = [tr-mean(tr); zeros(ns,1)];
         
    elseif strcmp(design_option,'MEANTR'),
         tr = mean(d, 2);
         tr = [ tr'; zeros(ns,1)];
         trnum = 0;
    end
    fttr=fft(tr); 
%%%%%%%%  Design the filter graphically
    reply = 'No';
    while strcmp(reply,'Yes')~=1,
        if exist('trspfig','var') && ishandle(trspfig),
            figure(trspfig)
        else
            trspfig = figure('name','Select Cutoff Frequecies',...
                'numbertitle','off', 'menubar', 'none', 'toolbar', 'none');
        end
        plot(f,abs(fttr(1:ns))); set(gca,'nextplot','add');
        axis tight
        xlabel('Frequency (MHz)');
        ylabel(['Amplitude spectrum of trace ' num2str(trnum) ]); 
        if filtertype == 1,
            [fcutoff,a]=getpoint(' (MHz)',' ');  %#ok<NASGU>
            fcutoff(1) = checkf(fcutoff(1),df,fn);
            Ws=fcutoff/fn;
            txt0 = ['Cutoff frequency = ' num2str(fcutoff,5) ' MHz'];
         elseif filtertype == 2,
            [fcutoff,a]=getpoint(' (MHz)',' ');   %#ok<NASGU>
            fcutoff(1) = checkf(fcutoff(1),df,fn);
            Ws=fcutoff(1)/fn;
            txt0 = ['Cutoff frequency = ' num2str(fcutoff,5) ' MHz'];
        else
            for i=1:2 
                [fcutoff,a]=getpoint(' (MHz)',' '); %#ok<NASGU>
                fcutoff(1) = checkf(fcutoff(1),df,fn);
                Ws(i) = fcutoff(1)/fn;
                ylimit = get(gca,'ylim');
                plot([fcutoff(1) fcutoff(1)],[0 ylimit(2)], ':r', ...
                    'linewidth',2)
                if i==1,
                    text1 = ['f_1=' num2str(fcutoff(1),5) ];
                    text(fcutoff,0.9*ylimit(2),text1,'fontsize',9,...
                        'fontwei','bold','horizontalal','right')
                end
                if i==2,
                    text2 = ['f_2=' num2str(fcutoff(1),5) ];
                    text(fcutoff,0.9*ylimit(2),text2,'fontsize',9,...
                        'fontwei','bold','horizontalal','left')
                end
                drawnow
            end
            txt0 = ['Cutoff frequencies : ' text1 ' MHz - ' text2 ' MHz'];
        end
% Compute filter coefficients
        if filtertype == 1,
            b = fir_f1(nf, Ws);
        elseif filtertype == 2,
            b = fir_f1(nf, Ws, 'high');
        elseif filtertype ==3,
            b = fir_f1(nf, Ws);
        elseif filtertype == 4,
            b = fir_f1(nf, Ws, 'stop' );
        end
% Test the filter
       if ~isempty(b),           %  An error in FIR_F1/2 may have 
            rb = size(b,1);      %  produced empty  b 
           if rb==1, 
               b=b'; 
           end;
           b = [b; zeros(ns-length(b),1); zeros(ns,1)]; %#ok<AGROW>
           H      = fft(b); 
           fttrf  = fttr.*H;
           fttrf  = fttrf.*conj(H);
           figure(trspfig); 
           cla
           plot(f ,abs(fttr(1:ns))/max(abs(fttr)),...
               f, abs(fttrf(1:ns))/max(abs(fttr)), f,abs(H(1:ns)),'r');
           xlabel('Frequency (MHz)'); 
           ylabel(['Amplitude spectrum of trace ' num2str(trnum) ]); 
           title(txt0);
           axis tight;
           drawnow;
           reply = questdlg('Is this OK? ');
           if strcmp(reply,'Cancel')==1,
               delete(trspfig); 
               dff=[];  Ws = [];  
               if exist('message','var') && ishandle(message),
                   delete(message)
               end
               return
           end
           if strcmp(reply,'No'),
               cla
           end
       else
           reply = 'No';
           cla
       end
    end
    delete(trspfig); 
    if exist('message','var') && ishandle(message),
        delete(message)
    end
end

%%% Cutoff frequencies are predefined and directly typed in 
if strcmp(design_option,'TYPEIN'),
    if filtertype == 1 || filtertype == 2, 
        answer = inputdlg('请给出截止频率(MHz)',...
            'F_FILTER: REQUEST',1);
        if isempty(answer), 
            dff = [];  Ws = [];
            return, 
        end;
        lb = checkcomma(answer);         % decimals will be correct even if , is typed instead of .
        Ws = str2num(lb)/fn;
    end;
    if filtertype == 3 || filtertype == 4, 
        answer = inputdlg({'请给出较低截止频率(MHz)',...
            '请给出较高截止频率(MHz)'}, ...
            'FREQUENCY_FILTER: REQUEST',1);
        if isempty(answer), 
            dff = [];  Ws = [];
            return, 
        end;
        lb = checkcomma(answer);          % decimals will be correct even if , typed instead of .
        Ws = [str2num(lb(1,:)) str2num(lb(2,:))]/fn;
    end;
    if filtertype == 1,
        b = fir_f1(nf, Ws);
    elseif filtertype == 2,
        b = fir_f1(nf, Ws, 'high');
    elseif filtertype == 3,
        b = fir_f1(nf, Ws);
    elseif filtertype == 4,
        b = fir_f1(nf, Ws, 'stop' );
    end;
    if isempty(b),
        erh = errordlg('Error in filter computations! Please try again!',...
                'FREQUENCY_FILTER : ERROR');
        uiwait(erh);
        dff = []; Ws = [];
        return
    end
    b = [b; zeros(ns-length(b),1); zeros(ns,1)];
end

%%%% Filter in the frequency domain twice to preserve phase
hw = helpdlg('Processing... Please wait ...');
drawnow
H     = fft(b);                     % Filter in frequency domain
H     = H * ones(1, ntr);           % Cast filter into a matrix operator 
fttr  = fft([d; zeros(ns,ntr)]);    % transform data 
fttr  = fttr.*H;                    % filter forwards
fttr  = fttr.*conj(H);              % filter backwards
dff   = real(ifft(fttr));           % invert to time domain
dff   = dff(1:ns,:);

close(hw);
return

function fo = checkf(fi,dx,xmax)
if fi < dx,
    fo = dx;
elseif fi > xmax,
    fo = xmax;
else
    fo = fi;
end
