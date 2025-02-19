function dagc = gainagc(d, dt, wagc)
%
%    AGC : Apply Automatic Gain Control to the traces of the GPR section
%          "d" by scaling the amplitude of the data at the centre of the
%          sliding window with respect to the RMS amplitude of the window.
%
% Usage  : dagc = gainagc( d, dt, wagc )
%
% Inputs :
%  wagc  : The length of the AGC window in ns
%    dt  : The data sampling interval
%     d  : 2-D GPR section
%
% Output :
%   dagc : The agc'ed GPR section
%
% Author : Andreas Tzanis,
%          Department of Geophysics, 
%          University of Athens
%          atzanis@geol.uoa.gr
%          (C) 2005, by Andreas Tzanis, all right reserved.
%
[ns, ntr] = size(d);
iwagc = floor(wagc/dt);                        % agc window in samples
% Double check the time window
if iwagc < 1, 
    erh = errordlg([ 'AGC 时窗太短了'],'GAINAGC : ERROR');
    uiwait(erh)
    dagc = [];
    return
end
if iwagc > ns/2, 
    erh = errordlg([ 'AGC 时窗 = ' num2str(wagc) ' 太长了! '],...
        'GAINAGC : ERROR');
    uiwait(erh)
    dagc = [];
    return
end
iwagc = round(iwagc/2);               % windows are symmetric, 
                                      % so work with half length
dagc = zeros(ns,ntr);
h = waitbar(0,'正在进行增益，请稍后!');
for itr = 1:ntr                       % loop over all traces
    
    tr = d(:,itr);                    % Current trace to process
    agcdata = zeros(ns,1);            % work array for agc'ed data         
% compute initial window for first datum 
    sum = 0.0;
    for i = 1:iwagc,
        val = tr(i);
        sum = sum + val*val;
    end
    nwin = 2*iwagc+1;
    rms = sum/nwin;
    if rms <= 0.0, 
        agcdata(1) = 0.0;
    else
        agcdata(1) = tr(1)/sqrt(rms);
    end
% ramping on 
    for i = 1:iwagc, 
        val = tr(i+iwagc);
        sum = sum + val*val;
        nwin= nwin + 1;
        rms = sum/nwin;
        if rms <= 0.0,
            agcdata(i) = 0.0;
        else 
            agcdata(i) = tr(i)/sqrt(rms);
        end
    end
% middle range -- full rms window 
    for i = iwagc+1 : ns-1-iwagc,
        val = tr(i+iwagc);
        sum = sum + val*val;
        val = tr(i-iwagc);
        sum = sum - val*val; 
        rms = sum/nwin;
        if rms <= 0.0,
            agcdata(i) = 0.0;
        else
            agcdata(i) = tr(i)/sqrt(rms);
        end
    end 
% ramping off 
    for i = ns-iwagc : ns,
        val = tr(i-iwagc);
        sum = sum - val*val;
        nwin = nwin-1;
        rms = sum/nwin;
        if rms <= 0.0,
            agcdata(i) = 0.0;
        else
            agcdata(i) = tr(i)/sqrt(rms);
        end
    end
% Trace finished - load onto output array
    dagc(:,itr) = agcdata;
    waitbar(itr/ntr,h);
end                                           % itr loop over traces
close(h);
return    