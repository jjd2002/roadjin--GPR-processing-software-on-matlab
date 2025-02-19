function markertr_out = editmarkers(markertr_in, Inpname, Infname, action)
%
% EDITMARKERS : Utility to edit the marker trace information created on 
%               importing a DZT or RD3 raw data file and to augment it 
%               with the co-ordinates of the marker trace locations.
%               Fot GPR systems wihtout survey wheels, this info may be 
%               used for marker interpolation to equal spacing and for  
%               for topographic /static corrections. It may also be used 
%               for creating 3-D radargram data volumes.
%
%       Usage : markertr_out = editmarkers(markertr_in,Inpname,Infname,action)
% 
%      Inputs : 
% markertr_in : Column vector with marker trace ID numbers, or [ n x 4 ]
%               matrix with each column representing the ID numbers and X,
%               Y, Z coordinates of marker traces in a local frame of
%               reference 
%     Inpname : String, the path of the input data file. Used for exporting
%               the marker trace matrix to a disk file.
%   Infname   : String, the name of the input data file. Used for exporting
%               the marker trace matrix to a disk file. 
%    action   : String, keyword
%               = 'editmode' invokes an in-house basic text editor in the
%                            form of a large edit box. 
%               = 'runmode' uses interactive data input by means of dialog
%                           boxes
%
%      Output : 
% markertr_out: Column vector or [ n x 4 ] matrix with marker coordinates
%               in a local frame of reference
%
%      Author : Andreas Tzanis, 
%               Department of Geophysics, 
%               University of Athens
%               atzanis@geol.uoa.gr
%               (C) 2005, Andreas Tzanis, all rights reserved
%

% global markertr_out result
global result
%%% First check if marker coordinates have been set in a previous run
markertr_out = markertr_in;               %% The dummy default
sf           = size(markertr_in);
%%% If not, create one by giving the number of the first trace to mark and 
%%% number of traces to skip (leave unmarked).
if isempty(markertr_in), 
    ask = questdlg(['No marker traces in this data set! '...
        'Do you wish to create now?'], 'EDITMARKERS: REQUEST','Yes');
    if isequal(ask,'Cancel') || isequal(ask,'No'), 
        return; 
    end;
    if strcmp(ask,'Yes'),
        msg = msgbox({'Forcing "editmode" execution!' ...
                      'Please type your data in the edit area' }, ...
                      'EDITMARKERS: NOTIFICATION', 'help');
        uiwait(msg)
        action = 'editmode';
    end
end

%%% First check the state of markertr_in
if strcmp(lower(action),'runmode') && sf(2) > 1, %#ok<STCI>
    ask = questdlg({'More than 1 columns exist in the marker trace matrix.' ... 
                    'The variable may have been assigned in a previous session !' ...
                    'What do you want to do? '}, 'Warning !!!', ...
                    'Check', 'Replace', 'Cancel', 'Check');
    if strcmp(ask,'Cancel'), 
        return; 
    end;
    if strcmp(ask,'Check'), 
        action = 'editmode';
    end
    if strcmp(ask,'Replace'),
        markertr_in = markertr_in(:,1);
        sf = size(markertr_in);
    end
end

%%% Now set the X,Y and Z coordinates of marker traces. If it exists,
%%% do not overwrite existing X,Y,Z information but only after inspection
if strcmp(lower(action),'runmode') && sf(2)==1, %#ok<STCI>
    ask = inputdlg('Are the markers Regularly, or Irregularly spaced? (r/i)',...
        'INPUT MARKER COORDINATES',1);
    if isempty(ask), 
        return; 
    end;
    mrkmode = ask{1};
    XM=[];
    if strcmp(lower(mrkmode),'i'), %#ok<STCI>
        for i=1:length(markertr_in),
            ask = inputdlg(['Give X-coordinate of Marker Trace ' ...
                num2str(markertr_in(i)) ' (' num2str(i) ' of ' ...
                num2str(sf(1)) ')'], 'MARKER X- COORDINATES',1);
            if isempty(ask), 
                return; 
            end;
            answer = checkcomma(ask);
            XM(i) = str2num(answer(1,:));
        end
        XM  = XM(:);          % force column vector
    elseif strcmp(lower(mrkmode),'r'), %#ok<STCI>
        ask = inputdlg({'Give X-location of first marker in m' ... 
                'Give X-spacing between Marker traces in m'}, ...
                'MARKER X- COORDINATES',1);
        if isempty(ask), 
            return; 
        end;
        answer        = checkcomma(ask);
        XM(1)         = str2num(answer(1,:));
        XM_spacing    = str2num(answer(2,:));
        for i=2:length(markertr_in), 
            XM(i) = XM(i-1) + XM_spacing;
        end;
        XM  =  XM(:);          % force column vector
    else
        erh = errordlg('Wrong answer, please try again!',...
            'EDITMARKERS : ERROR');
        uiwait(erh)
        return
    end

    YM=[];
    if strcmp(lower(mrkmode),'i'), %#ok<STCI>
        for i=1:length(markertr_in),
            ask = inputdlg(['Give Y-coordinate of Marker Trace '...
                num2str(markertr_in(i)), ' (' num2str(i) ' of ' ...
                num2str(sf(1)) ')'], 'MARKER Y- COORDINATES',1);
            if isempty(ask), 
                return; 
            end;
            answer = checkcomma(ask);
            YM(i) = str2num(answer(1,:));
        end
        YM  = YM(:);
    elseif strcmp(lower(mrkmode),'r'), %#ok<STCI>
        ask = inputdlg({'Give Y-location of first marker in m' ...
                'Give Y-spacing between Marker traces in m'}, ...
                'MARKER Y- COORDINATES',1);
        if isempty(ask), 
            return; 
        end;
        answer        = checkcomma(ask);
        YM(1)         = str2num(answer(1,:));
        YM_spacing    = str2num(answer(2,:));
        for i=2:length(markertr_in), 
            YM(i) = YM(i-1) + YM_spacing;
        end;
        YM  =  YM(:);
    else
        erh = errordlg('Wrong answer, please try again!',...
            'EDITMARKFILE : ERROR');
        uiwait(erh)
        return
    end

    ZM=[];
    if strcmp(lower(mrkmode),'i'), %#ok<STCI>
        for i=1:length(markertr_in),
            ask = inputdlg(['Give Elevation of Marker Trace '...
                num2str(markertr_in(i)) ' (' num2str(i) ' of ' ...
                num2str(sf(1)) ')'], 'MARKER ELEVATION',1);
            if isempty(ask), 
                return; 
            end;
            answer = checkcomma(ask);
            ZM(i) = str2num(answer(1,:));
        end
        ZM  =  ZM(:);
    elseif strcmp(lower(mrkmode),'r'), %#ok<STCI>
        ask = inputdlg({'Give Elevation of first marker in m' ... 
                'Give elevation spacing between Marker traces in m'}, ...
                'MARKER ELEVATION',1);
        if isempty(ask), 
            return; 
        end;
        answer        = checkcomma(ask);
        ZM(1)         = str2num(answer(1,:));
        ZM_spacing    = str2num(answer(2,:));
        for i=2:length(markertr_in), 
            ZM(i) = ZM(i-1) + ZM_spacing;
        end;
        ZM  =  ZM(:);
    else
        erh = errordlg('Wrong answer, please try again!',...
            'EDITMARKERS : ERROR');
        uiwait(erh)
        return
    end

%%% Done - Assign markertr_out
    markertr_out = [markertr_in XM YM ZM];
end

%%% Edit mode 
if strcmp(lower(action),'editmode'),  %#ok<STCI>
    MRKfig = findobj('tag','markercheckfigure');
    if ~isempty(MRKfig),
        figure(MRKfig);
        clf;
        matgprwindows('setup');
        matgprwindows('updateadd');
    else
        MRKfig = figure('name','Marker Information',...
            'tag','markercheckfigure','numbertitle','off', ...
            'menubar', 'none', ...
            'position',[500   300   500   600], ...
            'CreateFcn',['matgprwindows(''setup''); ' ...
                         'matgprwindows(''updateadd'');'],...
            'DeleteFcn',['matgprwindows(''updateremove'', ' ...
                         'findobj(''tag'',''markercheckfigure''));']);
    end
    datatitlestring = [Inpname Infname ];
    datainfotext    = '     ';
    rand('seed',sum(100*clock));
    %whitebg(MRKfig,[rand rand rand]); 
    c1 = 0.5;    c2 = 0.5;    c3 = 0.5; 
    whitebg(MRKfig,[c1 c2 c3]); 
    set(gca,'visible','off','drawmode','fast');
    text1 = text(0.45,1.05,'MARKER INFORMATION VIEWER') ;
    set(text1,'FontSize',13,'Color','k','FontWeight','bold',...
        'horizontalal','center')
% prepare to display header info on a GPR data file
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
    uicontrol( 'Style','frame', ...
        'Units','normalized', 'Position',framePosition, ...
        'BackgroundColor',[0.0 0.5 0.5]);
% Draw the text label
    labelPosition=[left top-labelheight (right-left) labelheight];
    uicontrol( 'Style','text', ...
        'Units','normalized', 'Position',labelPosition, ...
        'BackgroundColor',[0.0 0.5 0.5], 'ForegroundColor',[1 1 1], ...
        'String',datatitlestring, 'fontsize', 10, 'fontweight', 'demi');
    textPosition=[left bottom (right-left) top-bottom-labelheight-2*spacing];
    uicontrol( 'Style','edit', 'tag', 'markerinfobox', ...
        'Units','normalized', 'Max',18, 'String',datainfotext, ...
        'BackgroundColor',[1 1 1], 'Position',textPosition);
% Prepare pushbuttons 
    uicontrol('style','pushbutton','string','OK', 'tag','markerokbutton', ...
        'Units', 'Normalized', 'Position', [0.6 0.06 0.1 0.05], ...
        'backgroundcolor',[0 0.5 0.5], ...
        'fontweight','bold', 'fontsize',12,  ...
        'tooltip','Exit and accept changes', ...
        'Callback', {@do_getmarkertr,'OK'});
    uicontrol('style','pushbutton','string','Abort',...
        'tag','markerbadbutton', 'Units', 'Normalized', ...
        'Position', [0.75 0.06 0.12 0.05], ...
        'backgroundcolor',[0 0.5 0.5], ... 
        'fontweight','bold', 'fontsize',12,  ...
        'tooltip','Exit and reject changes', 'fontsize',10, ... 
        'Callback', {@do_getmarkertr,'Abort'});
% Display the marker info box and the marker info
    mrktxt = num2str(markertr_in);
    fontn = get(0,'FixedWidthFontName');
    do = ['set(findobj(''tag'',''markerinfobox''),''String'',mrktxt,'...
            '''horizontalal'',''left'',''fontname'',fontn);' ];
    er = 'editmarkers'; 
    eval(do,er);   
% Wait until the pushbuttons are activated, otherwise changes will vanish        
    waitfor(findobj('tag','markercheckfigure'))
end
% Abandon changes
if strcmp(result,'Abort'),
    return
end
%%% Optionally export to a marker file
if ~isempty(markertr_out),
    ask = questdlg('Save to a Marker File?', 'EDITMARKERS: REQUEST','Yes');
    if isequal(ask,'Yes'), 
        Markfile = [Infname(1:findstr(Infname,'.')) 'mrk'];
        fid=fopen([Inpname Markfile],'w');
        for i=1:length(markertr_in)
            fprintf(fid,' %d  %d  %d  %d\r\n',markertr_out(i,:));
        end
        fclose(fid);
    end
end
return

function do_getmarkertr(vhandle,eventdata,action)
% global markertr_out result
global result
if strcmp(action,'OK'),
    q = get(findobj('tag','markerinfobox'),'string');
    w = str2num(q);
    if isempty(w),
        warndlg({'It appears that some elements of the Marker trace ' ...
            'matrix are missing! Please check again or Abort!'},...
            'EDITMARKERS: WARNING');
        return
    end
    assignin('caller','markertr_out',w);   %markertr_out = w;
    result = 'OK';
end
if strcmp(action,'Abort'),
    result = 'Abort';
end
delete(findobj('tag','markercheckfigure'))
return
