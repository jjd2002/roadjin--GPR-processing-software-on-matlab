function vofh = get1dvelocity()
%
% GET1DVELOCITY : Imports 1-D velocity models for passing to 1-D migration
%                 and modelling routines. The models can be imported from
%                 file or typed in. The structure of the model is:
%                 NLAY  
%                 [Velocity, layer 1          Thickness, layer 1;
%                  Velocity, layer 2          Thickness, layer 2;
%                  ...                          ...  ;
%                  Velocity, layer NLAY-1     Thickness, layer NLAY-1;
%                  Velocity, basal halfspace          0   ];
%               => The thickness of the basal halfspace is always equal to
%                  zero! 
%               => A uniform halfspace is introduced as a 1-layer velocity
%                  model with a thickness equal to zero. 
%               => Basic error checking (e.g. for velocities exceeding the
%                  speed of light in free space) is done before the routine
%                  returns.  
%
%         Usage : vofh = get1dvelocity
%
%       Returns :
%          vofh : The 1-D velocity / thickness model 
%
%      Requires : checkcomma.m
%
%       Author  : Andreas Tzanis,
%                 Dept. of Geophysics,   
%                 University of Athens
%                 atzanis@geol.uoa.gr
%                 (C) 2005, Andreas Tzanis, allrights reserved
%

global ENVAR
%%%%%   Ask what to do 
importmodel = questdlg('请问是导入还是创建速度模型?',...
    'GET1DVELOCITY : REQUEST', '导入','创建','取消','创建');
if strcmp(importmodel,'取消'),            % Oops, abort!
    vofh = [];
    return
end

%%%%%   Structure will be read in from file
if strcmp(importmodel,'导入'),
    [vfname, vpname]= uigetfile('*.dat',...
        'GET1DVELOCITY : 1-D velocity model file name', ...
        ENVAR.currentworkdir);
    if vfname == 0,
        vofh = [];
        return
    end
    fid = fopen([vpname vfname],'r');
    nlay = fscanf(fid,' %d ',1);            % read number of layers
    vofh = zeros(nlay,2);
    for i=1:nlay-1
        dummy = fscanf(fid,' %f %f',2);
        vofh(i,1) = dummy(1);               % velocity of i'th layer
        vofh(i,2) = dummy(2);               % thickness of i'th layer
    end
    vofh(nlay,1)  = fscanf(fid,' %f ',1);   % velocity of basal halfspace
    vofh(nlay,2)  = 0;                      % thickness of basal halfspace
    
%%%%%   Structure will be given interactively! 
elseif strcmp(importmodel,'创建'),
    answer = inputdlg('给出包括基础半空间的层数',...
        'GET1DVELOCITY : REQUEST',1);
    if isempty(answer), 
        vofh = []; 
        return, 
    end;
    nlay = str2num(answer{1}); 
    vofh = zeros(nlay,2);                            
    clb = cell(nlay,1);                                        
    for i=1:nlay-1
        clb(i) = cellstr(['Layer ' num2str(i) ...
            ' Give Velocity in m/ns & Thickness in m']);
    end
    clb(nlay) = cellstr('半空间速度');
    for i=1:nlay-1
        cldef(i) = cellstr([num2str(vofh(i,1)) '  ' num2str(vofh(i,2))]);
    end
    cldef(nlay) = cellstr(num2str(vofh(nlay,1)));
    answer = inputdlg(clb,...
        'GET1DVELOCITY : Give layer velocities and thicknesses',1,cldef ); 
    if isempty(answer)
        vofh = []; 
        return      
    end
    lb = checkcomma(answer);       
    for i=1:nlay-1                 
        dummy = str2num(lb(i,:));
        vofh(i,1) = dummy(1);        % velocity of i'th layer      
        vofh(i,2) = dummy(2);        % thickness of i'th layer
    end
    vofh(nlay,1) = sscanf(lb(nlay,:),'%f'); % velocity of basal halfspace
    vofh(nlay,2) = 0;                       % thickness for basal halfspace
    
%%%%%   Abort ...
elseif strcmp(importmodel,'取消'),
    vofh = [];   
    return
end

%%%%%   Check validity of imported velocities
for i=1:nlay
    if vofh(i,1) <= 0 | vofh(i,1) > 0.2998,
        erh = errordlg('Impossible velocity value found! Please try again!',...
            'GET1DVELOCITY : ERROR')
        uiwait(erh);
        vofh = [];  
        return
    end
end

%%%%%   Save a newlly typed velocity model
if strcmp(importmodel,'创建'),
    answer = questdlg('是否保存速度模型 ?','GET1DVELOCITY : REQUEST',...
        'Yes');
    if strcmp(answer,'Yes'),
        [vfname, vpname]= uiputfile('*.dat',...
            'GET1DVELOCITY : Save 1-D velocity model',ENVAR.currentworkdir);
        if isequal(vfname,0),
            return
        end
        fid = fopen([vpname vfname],'w');
        fprintf(fid,' %d \n',nlay);              % write number of layers
        for i=1:nlay,
            fprintf(fid,' %f  %f\n',vofh(i,:));  % write i'th velocity and thickness
        end
        fclose(fid);
    elseif strcmp(answer,'导入') | strcmp(answer,'取消'),
        return
    end
end

return                   % END FUNCTION GET1DVELOCITY
