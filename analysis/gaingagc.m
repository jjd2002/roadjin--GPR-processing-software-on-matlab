function dagc = gaingagc( d, dt, EPS, wagc )
%
%   GAGC  : Apply Automatic Gain Control to the traces of the GPR section
%          "d" by scaling the amplitude of the data at the centre of a
%          Gaussian bell tapered sliding window with respect to the RMS
%          amplitude of the window. 
%
%   Usage : dagc  = gaingagc( d, dt, EPS, wagc )
%
%  Inputs :    d : The input GPR section
%             dt : Data sampling interval 
%            EPS : Noise level (to define the width of the Gaussian taper)
%           wagc : Length of the AGC window in ns
%
%  Output : dagc,  the aplified section
%
% Credits : This function was transcoded and modified from the SU function
%           do_gagc appearing in program sugain. Credits to CWP, Jack K.
%           Cohen, Brian Sumner, and Dave Hale 
%
% Author  : Andreas Tzanis,
%           Department of Geophysics, 
%           University of Athens
%           atzanis@geol.uoa.gr
%           (C) 2005, 2010, by Andreas Tzanis, all right reserved.
%

[ns, ntr] = size(d);
EPS       = sqrt(abs(log(EPS)));
iwagc     = floor(wagc/dt);                        % agc window in samples
% Double check he window
if iwagc < 1, 
    erh = errordlg([ 'AGC 时窗太短了'],'GAINGAGC : ERROR');
    uiwait(erh)
    dagc = []; 
    return
end
if iwagc > ns/2, 
    erh = errordlg([ 'AGC 时窗 = ' num2str(iwagc) '太长了! '],...
        'GAINGAGC : ERROR');
    uiwait(erh)
    dagc = [];
    return
end
% Compute Gaussian window weights 
w   = zeros(iwagc,1);
u   = EPS/iwagc;
u2  = u*u;
for i=1:iwagc
    w(i) = exp(-(u2*i*i));
end
d2 = zeros(ns,1);                    % Initialize sum of squares 
s  = zeros(ns,1);                    % Initialize weighted sum of squares

h = waitbar(0,'正在进行增益，请稍候');
dagc      = zeros(ns,ntr);
% loop over all traces
for itr = 1:ntr
    tr = d(:,itr);                    % Current trace to process
    agcdata = zeros(ns,1);            % work array for agc'ed data 
% agc itr'th trace
    for i = 1:ns
        val     = tr(i);
        d2(i) = val * val;
        s(i)  = d2(i);
    end
    for j = 1:iwagc-1;
        for i = j:ns
            s(i) = s(i) +( w(j)*d2(i-j+1));
        end
        k = ns - j;
        for i = 1:k
            s(i) = s(i) +( w(j)*d2(i+j));
        end
    end
    for i = 1:ns
        if ~s(i),
            agcdata(i) = 0.0;
        else
            agcdata(i) = tr(i)/sqrt(s(i));
        end
    end
    dagc(:,itr) = agcdata;
    waitbar(itr/ntr,h);
end                                    % itr loop over traces
close(h);
return    
