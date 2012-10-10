%%
%rot40_irc: analysis of IRC 
%
% Version 1.1
% Last modified  R O Zhurakivsky 2011-10-09
% Created        R O Zhurakivsky 2011-03-21

format compact
global flplot

clear 
atomsind

PLOTTYPE_deltaE=0;
PLOTTYPE_mu=1;
PLOTTYPE_Mulliken_q_H=2;
PLOTTYPE_Mulliken_q_A=3;
PLOTTYPE_rho=4;
PLOTTYPE_deltarho=5;
PLOTTYPE_epsilon=6;
PLOTTYPE_EHB=7;
PLOTTYPE_dAB=8;
PLOTTYPE_dAH=9;
PLOTTYPE_aAHB=10;
PLOTTYPE_RHH=11;
PLOTTYPE_aglyc=12;
PLOTTYPE_rho_CH=13;
PLOTTYPE_deltarho_CH=14;
PLOTTYPE_epsilon_CH=15;
PLOTTYPE_EHB_CH=16;
PLOTTYPE_dAH_CH=17;
PLOTTYPE_max=18;

flags={};
%-------------------------------------------------------------------
if 0
    workname='irc120321'%#ok
elseif 0
%    indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\AT'%#ok
%    indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\GC'%#ok
%    indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\AT\AT_wfn_tight_int'%#ok
%    indir='E:\work\Brovarets\120216_IRC_Hyp' %#ok
%    indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\Ade_H' %#ok
%    workname='irc_b3lyp' %#ok

%    indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\AT_tight'
%    workname='irc_b3lyp_tight'%#ok
    
%    indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\AT_int'
%    workname='irc_b3lyp_int'%#ok
    
%    indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\AT_tight_int'
%    indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\GC_tight_int'%#ok
%    indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\AT_tight_int_step1'
%    indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\AT_tight_int_step5'
%    indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp_Hyp_tight_int'
%    indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp_H_Hyp_O'
%    indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp_H_Hyp_O_tight_int'
%    indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp_O_Thy'
%    indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp-dimer'
%    indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp_O_Thy_tight_int'
    indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp_Hyp_tight_int'
    workname='irc_b3lyp_tight_int'%#ok
%    workname='irc120321_freq'
elseif 0
    indir='E:\work\Brovarets\120411_irc_AT_GC\mp2'%#ok
%    workname='irc120411_mp2'
    workname='irc120425_mp2_tight'%#ok
elseif 0
%    indir='E:\work\Brovarets\1204_irc_dnucl_SN\dAdo'%#ok
%    workname='irc_b3lyp'%#ok
%    indir='E:\work\Brovarets\1204_irc_dnucl_SN\b3lyp_631gpd'%#ok
%    indir='E:\work\Brovarets\1204_irc_dnucl_SN\dCyd_631gpd'%#ok
%    indir='E:\work\Brovarets\1204_irc_dnucl_SN\dThd_631gpd'%#ok
%    indir='E:\work\Brovarets\1204_irc_dnucl_SN\dAdo_631Gdp'%#ok
%    indir='E:\work\Brovarets\1204_irc_dnucl_SN\dUrd_631gdp'%#ok
%    indir='E:\work\Brovarets\1204_irc_dnucl_SN\bohill'%#ok
    workname='irc_b3lyp_631gdp'%#ok
elseif 0
%    indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\AT_eps4_tight_int'%#ok
%    indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\AT_eps4_tight_int_2parts'%#ok
%    indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\GC_eps4_tight_int'%#ok
%    indir='E:\work\Brovarets\120216_IRC_Hyp\tight_int_eps4'%#ok
%    indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp_Hyp_eps4_tight_int'%#ok
    indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp_O_Thy_eps4_tight_int'%#ok
    workname='irc_b3lyp_tight_eps4'%#ok
elseif 0
    complextype='HypHyp';
    if 0
        indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp_Hyp_tight_int'
        workname='irc_b3lyp_tight_int'%#ok
    elseif 0
        indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp_Hyp_eps4_tight_int'
        workname='irc_b3lyp_tight_eps4'%#ok
    elseif 1
        indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp-dimer'
        workname='irc_b3lyp_631gdp'%#ok
        flags.semilogy_for_ellipticity=1;
    end
elseif 1
    complextype='HypCyt';
    indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp_Cyd\tight_int_eps4'
    workname='irc_b3lyp_tight_eps4'%#ok
end

flplot=1    %#ok
fl_tosave = 1 %#ok
fl_savepics = 1 %#ok
fl_recreate_matfile = 0 %#ok
fl_reload_sp = 0   %#ok
fl_reload_freq = 0 %#ok
fl_reload_extout = 0 %#ok

%mode='gjf'
mode='anal' %#ok
%-------------------------------------------------------------------

diaryfname0=[indir filesep 'logfile'];
diaryfname=diaryfname0;
for i=2:10000
  if ~(exist(diaryfname,'file')==2)
    break
  end
  diaryfname = [diaryfname0 int2str(i)];
end
diary(diaryfname)

disp(['Time: ' datestr(now)]);
indir %#ok
diaryfname %#ok

%savemode.gsxyz = 1;
gtemplname=[workname '_templ.gjf']  %#ok
fullgtemplname=[CD.templatesdir filesep gtemplname] %#ok

fl_empty_figures = 0;

odir=indir;
if exist(odir,'dir')~=7
   mkdir(odir);
end

lfiles = dir(strcat(indir,filesep,'*.list'));
num_listfiles = size(lfiles,1);
if ~num_listfiles
  error('No LIST files to analyse found');
end

for l_ind=1:num_listfiles
    
    dlm=strfind(lfiles(l_ind).name,'.');
    fnameshort = lfiles(l_ind).name(1:(dlm-1));
    workdbname = [indir filesep fnameshort '.mat']; %#ok
    dlm=strfind(workdbname,'.');
    workdbnamebkp=[workdbname(1:dlm(end)-1) '~' workdbname(dlm(end):end)];
    xlsfile = [indir filesep fnameshort '.xls'] %#ok
    
    psfile = [indir filesep fnameshort '.ps'];
    psfilebkp=[psfile(1:dlm(end)-1) '~' psfile(dlm(end):end)];
    svgfile = [indir filesep fnameshort '.svg'];
    figfile = [indir filesep fnameshort '.fig'];

%-------- additional parameters ---------------------
%dlm=strfind(indir,filesep);
%basedir = indir(dlm(end)+1:end);
aimbondlist=[];
bondlist=[];
anglist=[];
if strcmp(complextype,'HypCyt')
    %glicosidic parameters
    disp(['HypCyt pair job detected using predefined bondlist, anglist and aimbondlist'])
    bondlist = [12 24]; %R_HH
    anglist  = [ 1 12 24; 12 24 13]; %alpha1, alpha2
    molpind={};
    molpind{6}='C6';
    molpind{7}='O6';
    molpind{8}='N1';
    molpind{9}='C2';
    molpind{19}='N4';
    molpind{21}='N3';
    molpind{23}='O2';
elseif strcmp(complextype,'HypHyp')
    %glicosidic parameters
    disp(['HypHyp pair job detected using predefined bondlist, anglist and aimbondlist'])
    bondlist = [13 25]; %R_HH
    anglist  = [ 3 13 25; 13 25 24]; %alpha1, alpha2
    molpind={};
    molpind{7}='C6';
    molpind{8}='O6';
    molpind{9}='N1';
    molpind{11}='C2';
    molpind{16}='N1';
    molpind{17}='C2';
    molpind{27}='C6';
    molpind{28}='O6';
elseif strcmp(complextype,'HypThy')
    %glicosidic parameters
    disp(['HypThy pair job detected using predefined bondlist, anglist and aimbondlist'])
    bondlist = [13 21]; %R_HH
    anglist  = [ 1 13 21; 13 21 20]; %alpha1, alpha2
    molpind={};
    molpind{6}='C6';
    molpind{7}='O6';
    molpind{8}='N1';
    molpind{14}='O4';
    molpind{16}='C4';
    molpind{17}='N3';
    molpind{18}='C2';
    molpind{19}='O2';
elseif strcmp(complextype,'AT')
    disp(['AT pair job detected using predefined bondlist, anglist and aimbondlist'])
    %glicosidic parameters
    bondlist = [11 24]; %indexes of edge glicosidic atoms
    anglist  = [ 1 11 24; 15 24 11];
    aimbondlist = [8 14];
elseif strcmp(complextype,'GC')
    %glicosidic parameters
    disp(['GC pair job detected using predefined bondlist, anglist and aimbondlist'])
    bondlist = [14 26];
    anglist  = [ 1 14 26; 15 26 14];
end
%-------- additional parameters end -----------------
    
    if fl_recreate_matfile || exist(workdbname,'file')~=2

        disp(['Reading ' lfiles(l_ind).name ]);
        [outfiles,irc_directions,irc_shifts]=textread([indir '\' lfiles(l_ind).name],'%s%f%f');

        if ~numel(outfiles)
            disp('No Gaussian output files found. Skipping.');
            continue
        end
        ircdb={};

        for out_ind=1:numel(outfiles)


            dlm=strfind(outfiles{out_ind},'.');
            out = [];
            out.fnameshort = outfiles{out_ind}(1:(dlm-1));
            out.fnamefull = outfiles{out_ind}(1:(dlm(end)-1));

            out.ffname = strcat(indir,filesep,outfiles{out_ind}); 

            out.worktitle = lower(out.fnamefull);
            irc_direction = irc_directions(out_ind);
            if out_ind<=numel(irc_shifts)
                irc_shift = irc_shifts(out_ind);
            else
                irc_shift = 0;
            end

        %    try
                disp(['Loading file: ' out.ffname])
                %processing single OUT file
                fid=fopen(out.ffname,'r'); % don't use text mode - this is dangerous for Unix files ;)

                fl_ircstart = 0;
                fl_good_geom_found = 0;

                frewind(fid);

                clear ms0
                while 1 
                    tline = fgetl(fid);
                    if ~ischar(tline), break, end

                    if ~isempty(strfind(tline,'Cartesian coordinates read from the checkpoint file')) || ...
                       ~isempty(strfind(tline,'Redundant internal coordinates taken from checkpoint file'))
                        disp('Initial coordinates found'); %point #0

                        for i=1:2, tline=fgets(fid); end  
                        [ms0,status]=extrgeom(fid,1,0);
                        if status
                            disp([worktitle,': ',lastwarn]);
                            fclose(fid);
                            break
                        end
                    
                        ms0.point = 0;
                        ms0.irc = 0 + irc_shift;

                    elseif ~isempty(strfind(tline,'Energy From Chk'))
                        %Energy From Chk =   -921.7421354         
                        [xxx,xxx2]=strread(tline,'%s%f','delimiter','=');
                        ms0.energy = xxx2;
                        
                    end
                    
%                    if ~isempty(strfind(tline,'Start new reaction path calculation'))
                    if ~isempty(strfind(tline,'Calculating another point on the path'))
                        fl_ircstart = 1;
                        disp('found "Start IRC calculation"');
                        
                        ms0.desc = [out.fnameshort '_#' sprintf('%03d',ms0.point)];
                        order=1:ms0.atomnum;
                        %if strcmp(mode,'gjf') 
                        if exist([odir filesep 'gjf' filesep ms0.desc '.gjf'],'file')~=2 %create input Gaussian file
                            if exist([odir filesep 'gjf'],'dir')~=7
                                mkdir(odir,'gjf');
                            end
                            savemolgs([odir filesep 'gjf'],ms0,3,order,fullgtemplname); %Gaussian with XYZ
                        end
                        if isempty(ircdb)
                            ircdb=ms0;
                        else
                            ircdb(end+1) = orderfields(ms0,ircdb(1)); %#ok
                        end
                        addstr = [':' num2str(ms0.irc,4)];
                        disp(['Point ' int2str(ms0.point) addstr ' done']);
                        
                        break
                    end
                end
                if (~fl_ircstart)
                    disp('No IRC found');
                    exit;
                end

                clear ms0
                buf={};
                fl_buf_ready = 0;
                while 1
                    tline = fgetl(fid);
                    if ~ischar(tline), break, end

                    if ~isempty(strfind(tline,'Calculating another point on the path'))
        %                disp(['Calculating another point on the path']);
                        fl_buf_ready = 1; %another point found, so our buffer is filled by actual data
                    elseif    ~isempty(strfind(tline,'Maximum number of steps reached'))
                        fl_buf_ready = 1; 
                    elseif    ~isempty(strfind(tline,'PES minimum detected on this side of the pathway'))
                        fl_buf_ready = 1; 

                    elseif ~isempty(strfind(tline,'Input orientation')) 
                        clear ms0
                        for i=1:2, tline=fgets(fid); end  %skip 2 lines to reach "Center     Atomic" line
                        [ms0,status]=extrgeom(fid,1,0);
                        if status
                            disp([worktitle,': ',lastwarn]);
                            fclose(fid);
                            break
                        end

                    elseif ~isempty(strfind(tline,'SCF Done')) 
                    % SCF Done:  E(RB3LYP) =  -941.576443161     A.U. after   12 cycles
                       [xxx,xxx2]=strread(tline,'%s%f','delimiter','=');
                       ms0.energy = xxx2;

                    elseif ~isempty(strfind(tline,'Delta-x Convergence NOT Met'))
        %                disp(['Delta-x Convergence NOT Met']);
                    elseif ~isempty(strfind(tline,'Delta-x Convergence Met'))
        %                disp(['Delta-x Convergence Met']);
                        fl_good_geom_found = 1;
                    end

                    
                    
                    if fl_ircstart || fl_good_geom_found
                        buf{end+1}=tline; %#ok
                    end
                    
                    if fl_buf_ready
                        fl_ircstart = 0;
                        %%%
                        % Point Number:   2          Path Number:   1
                        pat = '\s+Point\sNumber:\s+(\d+)\s+Path\sNumber:\s+(\d+)';
                        A = regexp(buf, pat, 'tokens','once');
                        A=[A{:}];
                        ms0.point = str2double(A{1});

                        %NET REACTION COORDINATE UP TO THIS POINT =    0.26904
                        pat = '\s*NET\sREACTION\sCOORDINATE\sUP\D+([\d\.]+)';
                        A = regexp(buf, pat, 'tokens','once');
                        A=[A{:}];
                        ms0.irc = str2double(A{1})*irc_direction + irc_shift;

                        if fl_good_geom_found
                            ms0.desc = [out.fnameshort '_#' sprintf('%03d',ms0.point)];
                            order=1:ms0.atomnum;
                            %if strcmp(mode,'gjf') 
                            if exist([odir filesep 'gjf' filesep ms0.desc '.gjf'],'file')~=2 %create input Gaussian file
                                if exist([odir filesep 'gjf'],'dir')~=7
                                    mkdir(odir,'gjf');
                                end
                                savemolgs([odir filesep 'gjf'],ms0,3,order,fullgtemplname); %Gaussian with XYZ
                            end

                            if isempty(ircdb)
                                ircdb=ms0;
                            else
                                ircdb(end+1) = orderfields(ms0,ircdb(1)); %#ok
                            end
                        end

                        addstr = '';
                        if isfield(ms0,'irc')
                            addstr = [':' num2str(ms0.irc,4)];
                        end
                        disp(['Point ' int2str(ms0.point) addstr ' done']);

                        buf={};
                        fl_good_geom_found = 0;
                        clear ms0
                        fl_buf_ready = 0;

                    end
                    
                    
                end %while 1 %cycle over IRC points


                fclose(fid);
        %    catch
        %       disp(['error: ',lasterr]);
        %       fclose(fid);
        %    end
        end %i=1:numel(outfiles)

    else %if MAT file exists , load it

        disp(['Loading file: ' workdbname])
        load(workdbname,'ircdb');
        
    end %if exist(workdbname,'file')~=2

 
%--------------------------------------------------------------------------
    pcolor=[{'r'} {'g'} {'b'} {'k'} {'m'} {'c'} ...
            {[.5412 .1686 .8863]}  ...
            {[0 .5 0]} {[.5 0 0]} {[0 0 .5]} {[.5 .5 0]} {[.5 0 .5]} {[0 .5 .5]} {[.5 .5 .5]} ...
            {[.5 .1 .9]} {[.8 .2 .2]} {[.8 .8 .2]}...
            {[.9 .4 .9]} {[.2 .4 .6]} {[.6 .4 .6]} {[.6 .2 .2]} {[.8 .2 .8]} ...
            {[.2 .8 .8]} {[.2 .8 .2]} {[.2 .2 .8]} {[.4 .9 .1]} {[.1 .3 .6]} {'y'} {[.75 .75 .75]} {[.2745 .5098 .7059]}];
    psign='dhxov^<>p+*.';
    pstyle=[{'-'} {'--'} ];
    lw = 0.5; %line width
    ms = 4; %marker size (points)
    msc = 5; %marker size for crosses (points)

    bdwidth = 2; %figure border (points)
    topbdwidth = 20; %figure top border (points)
    set(0,'Units','pixels') 
    scnsize = get(0,'ScreenSize');
    pos  = [bdwidth,... 
            bdwidth,...
            scnsize(3) - 2*bdwidth,...
            scnsize(4) - (topbdwidth + bdwidth)];

    bwl = 0.02; %border width left
    bwr = 0.005;
    bht = 0.03; %border height top
    bhb = 0.05;
    main_h = 0.75; %height of lower pictures
    fw2 = 0.15; %width of right field used for legends 
    fw1 = 0.5*(1-fw2); %width of ordinary field 

    %figure and axes for plotting to PS file
    f_tmp = figure('PaperOrientation','portrait',...
           'PaperType','A4',...
           'PaperUnits','centimeters',...
           'PaperPositionMode','manual',...
           'PaperPosition', [1 1 16 8],...
           'Position',pos,...
           'Visible','on');
%           'Visible','off');
    a_tmp = axes();

    if exist(psfile,'file')
        movefile(psfile,psfilebkp);
    end
    

    data = {[]};
       
%----------------------Energy and dipole moment----------------------------------------------------
    figure(f_tmp);
    [x,sort_ind]=sort([ircdb.irc]);
    y=[ircdb.energy];
    y=y(sort_ind);
    y=(y-min(y))*CC.encoef;
    
    col = 1;
    data(1,col)={'ind'};                data(2:numel(ircdb)+1,col) = num2cell(sort_ind); col=col+1;
    data(1,col)={'filename'};           data(2:numel(ircdb)+1,col) = {ircdb(sort_ind).desc}; col=col+1;
    data(1,col)={'IRC, a.m.u.^{1/2}*bohr'}; data(2:numel(ircdb)+1,col) = num2cell(x); col=col+1;

    txt = 'energy, kcal/mol';
    data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
    
    h_plot = plot(x,y,'r.-','LineWidth',lw);
    saved_axis(1,:)=[min(x) max(x) 0 1.05*max(y)];
    axis(saved_axis(1,:));

    ylabel('\Delta E, kcal/mol');
    
    set(gca,'Box','on','XMinorTick','on','YMinorTick','on');
    xlabel('IRC, a.m.u.^{1/2}*bohr');
    title('');
    print(f_tmp,'-dpsc2', '-append', '-r300', psfile);
%    clf(f_tmp);
    
    
    %--------------------------------------------------------------------------
    %analyzing Gaussian SP output files for dipole moment and charges
    dirlist = dir([indir filesep 'sp' ]);
    if ~isempty(dirlist)
    for iext=1:numel(ircdb)
        
        if ~fl_reload_sp && isfield(ircdb(iext),'DM') && ~isinf(ircdb(iext).DM) 
            continue;
        end
        
        fl_out_file_found = 0;
        for idir=1:numel(dirlist)
            if strncmp(dirlist(idir).name, ircdb(iext).desc, numel(ircdb(iext).desc)) &&...
                strcmp(dirlist(idir).name(end-3:end), '.out')
            
                fl_out_file_found = 1;
                file = ['sp' filesep dirlist(idir).name ];
                file_fullname = [indir filesep file ];
                
                break;
            end
        end    
        
        if ~fl_out_file_found
            disp(['!!! SP file is not found for : ' ircdb(iext).desc])
            
        else

            disp(['...Loading SP file: ' file])
            fid=fopen(file_fullname,'r'); % don't use text mode - this is dangerous for Unix files ;)
            if fid==-1
              disp(['!!!Can''t open file ' file])
              continue;
            end
            
            while 1
              tline = fgetl(fid);
              if feof(fid), break, end
              if ~isempty(strfind(tline,'Dipole moment (field-independent basis, Debye)')), 
                  tline = fgetl(fid);
                  %    X=    -0.9740    Y=    -0.3414    Z=     0.0074  Tot=     1.0321

                  F='(-?\d*\.*\d+E?[+-]?\d*)';
                  pat = ['X=\s+' F '\s+Y=\s+' F '\s+Z=\s+' F '\s+Tot=\s+' F];
                  A = regexp(tline, pat, 'tokens','once');
                  ircdb(iext).DM = sscanf(A{4},'%f'); %#ok
%                  break, 
              elseif ~isempty(strfind(tline,'Mulliken atomic charges:')), 
                  tline = fgetl(fid);
                  
                  A=fscanf(fid,'%d %s %f',[3,inf]);
                  if isempty(A)
                    warning('error in output file - couldn''n load charges. skipping');%#ok
                    break
                  end
                  ircdb(iext).('mcharge')=A(3,:);%#ok %Mulliken atomic charges
                  
              end
            end

            fclose(fid);
        end
        
    end %iext=1:numel(ircdb)
    end

    %--------------------------------------------------------------------------
    if fl_reload_freq || ~isfield(ircdb,'ZPE')
    %analyzing Gaussian FREQ output files for GEC (obsolete)
    dirlist = dir([indir filesep 'freq' ]);
    for iext=1:numel(ircdb)
        
        fl_out_file_found = 0;
        for idir=1:numel(dirlist)
            if strncmp(dirlist(idir).name, ircdb(iext).desc, numel(ircdb(iext).desc)) &&...
                strcmp(dirlist(idir).name(end-3:end), '.out')
            
                fl_out_file_found = 1;
                file = ['freq' filesep dirlist(idir).name ];
                file_fullname = [indir filesep file ];

                break;
            end
        end    
        
        if fl_out_file_found

            disp(['...Loading FREQ file: ' file])
            fid=fopen(file_fullname,'r'); % don't use text mode - this is dangerous for Unix files ;)
            if fid==-1
              disp(['!!!Can''t open file ' file])
              continue;
            end
            
            while 1
              tline = fgetl(fid);
              if feof(fid), break, end
              if ~isempty(strfind(tline,'Zero-point correction')) 
                 [xxx,ZPE]=strread(tline,'%s%f','delimiter','=');
                 ircdb(iext).ZPE=ZPE*CC.encoef;%#ok  %kcal/mol
              end
              if ~isempty(strfind(tline,'Thermal correction to Gibbs Free Energy')) 
                 [xxx,GEC]=strread(tline,'%s%f','delimiter','='); %Gibbs energy correction
                 ircdb(iext).GEC=GEC*CC.encoef;%#ok  %kcal/mol 
              end

            end

            fclose(fid);
        end
        
    end %iext=1:numel(ircdb)
    end

    %--------------------------------------------------------------------------
    %analyzing AIMAll output files
    %extout_processed = 0;
    for iext=1:numel(ircdb)
        if ~fl_reload_extout && isfield(ircdb(iext),'AIM') && ~isempty(ircdb(iext).AIM)
            continue;
        end

        mol = ircdb(iext);
        aimfile = [ 'wfn' filesep mol.desc '.extout' ];
        aimfile_fullname = [indir filesep aimfile];
    
        if exist(aimfile_fullname,'file')==2

            disp(['...Loading AIMAll file: ' aimfile])
            fid=fopen(aimfile_fullname,'r'); % don't use text mode - this is dangerous for Unix files ;)
            if fid==-1
              disp(['!!!Can''t open file ' aimfile])
              continue;
            end

            ms0=struct();
            maxind=0;
            while 1
              tline = fgetl(fid);
              if feof(fid), break, end
              if ~isempty(strfind(tline,'The nuclear coordinates')), 
                  break, 
              end
            end
            A=[];
            line_ind=0;
            while 1
              tline = fgetl(fid);
              if feof(fid), break, end
              if isempty(tline)
                  break;
              end
              line_ind=line_ind+1;

              %      H29           3.3920596300E+00 -4.2126501200E+00 -2.5792230500E+00
              pat = '^\s*([A-za-z]+)\d+\s+(-?\d*\.*\d*E?[+-]?\d*)\s+(-?\d*\.*\d*E?[+-]?\d*)\s+(-?\d*\.*\d*E?[+-]?\d*)'; %!!changed
              A = regexp(tline, pat, 'tokens','once');
              if isempty(A)
                  warning('Atoms coords not found');%#ok
              else
                  ms0.labels(line_ind)=A(1);
                  ms0.x(line_ind)=sscanf(A{2},'%f'); %!!changed
                  ms0.y(line_ind)=sscanf(A{3},'%f'); %!!changed
                  ms0.z(line_ind)=sscanf(A{4},'%f'); %!!changed
              end
            end


            if ~isempty(strcmpcellar(ms0.labels,''))
              warning('rot28_3_importAIMAll:emptylabels','Empty labels found!');
            end
            %In WFN file coordinates are in Bohrs
            ms0.x=ms0.x'*CC.l*1e10;
            ms0.y=ms0.y'*CC.l*1e10;
            ms0.z=ms0.z'*CC.l*1e10;
            ms0.atomnum = uint16(length(ms0.labels));
            ms0.ind = ((maxind+1):(maxind+ms0.atomnum))';

            ms0 = createbondtable(ms0);

            %TBD check mol and ms0 strucrures identity

            clear('CPs');
            CPind=0;

            while 1
              tline = fgetl(fid);
              if feof(fid), break, end
%              if ~isempty(strfind(tline,'CP#')), break, end
              if ~isempty(strfind(tline,'New critical point found')), break, end
              
            end
            fl_searchended = 0;
            while 1
                if feof(fid), break, end
                buf={};
                CP=struct([]);

                CPind=CPind+1;
                fl_coordsfound = 0;
                while 1
                  buf{end+1}=tline;%#ok

                  tline = fgetl(fid);
                  if feof(fid), break, end
                  if ~isempty(strfind(tline,'Electron Density Critical Point Analysis of Molecular Structure'))
                     fl_searchended = 1;
                     break;
                  end

                  if ~isempty(strfind(tline,'Coordinates of critical point and distance from molecular origin'))
                     fl_coordsfound = 1;
                  end

                  if fl_coordsfound
    %                [a1,CP.x] = strread(tline,'%s%f','delimiter','=');
                    A=fscanf(fid,'%s %s %f',[3,3]);
                    if isempty(A)
                        warning('CP coords not found');%#ok
                    else
                        CP(1).ind=CPind;
                        %In WFN file coordinates are in Bohrs %zhr091208
                        CP.x=A(3,1)*CC.l*1e10;
                        CP.y=A(3,2)*CC.l*1e10;
                        CP.z=A(3,3)*CC.l*1e10;
                    end

                    fl_coordsfound = 0;
                  end


%                  if ~isempty(strfind(tline,'CP#'))
                  if ~isempty(strfind(tline,'New critical point found'))
                      
                     break
                  end
                end %while

                if fl_searchended
                    break;
                end


                %Rho(r)                   2.5659186050E-01
                pat = '\s+Rho\(r\)\s+(-?\d*\.*\d+E?[+-]?\d*)';
                A = regexp(buf, pat, 'tokens','once');
                A=[A{:}];
                if isempty(A)
                    warning(['CP Rho for CP #' int2str(CPind) ' not found']);%#ok
                else
                    CP.rho=sscanf(A{1},'%f');
                end

                %DelSq(Rho(r))           -6.1929738296E-01
                pat = '\s+DelSq\(Rho\(r\)\)\s+(-?\d*\.*\d+E?[+-]?\d*)';
                A = regexp(buf, pat, 'tokens','once');
                A=[A{:}];
                if isempty(A)
                    warning(['CP DelSqRho for CP #' int2str(CPind) ' not found']);%#ok
                else
                    CP.DelSqRho=sscanf(A{1},'%f');
                end

                % Ellipticity:  1.7429996616E+00
                pat = '\s+Ellipticity:\s+(-?\d*\.*\d+E?[+-]?\d*)';
                A = regexp(buf, pat, 'tokens','once');
                A=[A{:}];
                if isempty(A)
                    CP.BondEl = NaN;
                else
                    CP.BondEl=sscanf(A{1},'%f');
                end

                % V(r)                    -2.6699369723E-01
                pat = '\s+V\(r\)\s+(-?\d*\.*\d+E?[+-]?\d*)';
                A = regexp(buf, pat, 'tokens','once');
                A=[A{:}];
                if isempty(A)
                    warning(['CP V(r) for CP #' int2str(CPind) ' not found']);%#ok
                else
                    CP.V=sscanf(A{1},'%f');
                end

                % G(r)                     7.0431800070E-03
                pat = '\s+G\(r\)\s+(-?\d*\.*\d+E?[+-]?\d*)';
                A = regexp(buf, pat, 'tokens','once');
                A=[A{:}];
                if isempty(A)
                    warning(['CP G(r) for CP #' int2str(CPind) ' not found']);%#ok
                else
                    CP.G=sscanf(A{1},'%f');
                end

                % K(r)                    -1.4615504038E-03
                pat = '\s+K\(r\)\s+(-?\d*\.*\d+E?[+-]?\d*)';
                A = regexp(buf, pat, 'tokens','once');
                A=[A{:}];
                if isempty(A)
                    warning(['CP K(r) for CP #' int2str(CPind) ' not found']);%#ok
                else
                    CP.K=sscanf(A{1},'%f');
                end

                % L(r)                    -8.5047304108E-03
                pat = '\s+L\(r\)\s+(-?\d*\.*\d+E?[+-]?\d*)';
                A = regexp(buf, pat, 'tokens','once');
                A=[A{:}];
                if isempty(A)
                    warning(['CP L(r) for CP #' int2str(CPind) ' not found']);%#ok
                else
                    CP.L=sscanf(A{1},'%f');
                end

                %Type
                pat = '(Point is a Bond Critical Point)';
                A = regexp(buf, pat, 'tokens','once');
                A=[A{:}];
                if isempty(A)
    %                warning(['CP Type for CP #' int2str(CPind) ' not found']);
                    CP(1).type='3,-3';
                else
                    CP.type='3,-1';
                end

                %CP atoms
                pat = 'Bond path linked to nuclear attractor\s+[A-Za-z]+([0-9]+)';
                A = regexp(buf, pat, 'tokens');
                A=[A{:}];
                if isempty(A)
    %                warning(['CP atoms for CP #' int2str(CPind) ' not found']);
                    CP.atoms=[];
                else
                    CP.atoms=[A{:}];
                end
                CPs(CPind)=CP;%#ok
            end %while 1
            fclose(fid);
    
            clear AIM
            AIM.is_hbond=[];
            AIM.atoms=[];
%            AIM.desc={};
            AIM.ro=[];
            AIM.DelSqRho=[];
%            AIM.pinds=[];
            AIM.BondEl=[];
            AIM.V=[]; 
            AIM.G=[]; 
            AIM.K=[];
            AIM.L=[];
            for i=1:numel(CPs)
               if  ~strcmp(CPs(i).type,'3,-1')
                   continue;
               end

               atoms=[sscanf(CPs(i).atoms{1},'%d') sscanf(CPs(i).atoms{2},'%d')];
               atoms=sort(atoms);
               aa = ms0.btB(find( ms0.btA==atoms(1)));%#ok %atoms connected to atoms(1) 
               if sum(aa==atoms(2))>0 %exclude covalent/ionic bonds
%                    continue
                   AIM.is_hbond(end+1,:)=0;
               else
                   AIM.is_hbond(end+1,:)=1;
               end

               AIM.atoms(end+1,:)=atoms;
%               AIM.desc(end+1)={bondstr};
               AIM.ro(end+1)=CPs(i).rho;
               AIM.DelSqRho(end+1)=CPs(i).DelSqRho;
               AIM.V(end+1)=CPs(i).V;
               AIM.G(end+1)=CPs(i).G;
               AIM.K(end+1)=CPs(i).K;
               AIM.L(end+1)=CPs(i).L;
               AIM.BondEl(end+1)=CPs(i).BondEl;
%               AIM.pinds(end+1,1:numel(HBatomspinds))=HBatomspinds;

            end
            
            ircdb(iext).AIM = AIM;%#ok
            
        else
            disp(['!!! can''t found: ' 'wfn' filesep mol.desc '.extout'])
        end %if exist(aimfile,'file')==2

    end %for iext=1:numel(ircdb)
%    end %if ~isfield(ircdb,'AIM')

    %creating array with indexes of atoms around 'H-bond like BCPs' for all structures
    %along IRC
    uniqCPind = []; 
    uniqCP_irc_struct_ind = [];
    for ind_irc=1:numel(ircdb) %over IRC points
        mol = ircdb(ind_irc);
        if ~isfield(mol,'AIM') || isempty(mol.AIM)
            continue;
        end
        for iCP=1:size(mol.AIM.atoms,1) %cycle over CPs
            if ~mol.AIM.is_hbond(iCP)
                continue
            end
            atoms = mol.AIM.atoms(iCP,:);
            if isempty(uniqCPind) || ~any(sum(uniqCPind==repmat(atoms,size(uniqCPind,1),1),2)==2)
                uniqCPind(end+1,:) = atoms;%#ok
                uniqCP_irc_struct_ind(end+1) = ind_irc;%#ok
            end
        end
    end %for ind_irc=1:numel(ircdb)
    
    
    allHBind = []; %array with indexes of all found AH...B bonds
    for icp=1:size(uniqCPind,1)

       ms0 = ircdb(uniqCP_irc_struct_ind(icp));
       ms0 = createbondtable(ms0);

       atoms = uniqCPind(icp,:);
       HBatoms = [];
       if any(ms0.labels{atoms(1)}~='H') && any(ms0.labels{atoms(2)}~='H') %none of atoms are H
            disp(['Contact with atoms ' int2str(atoms(1)) ',' int2str(atoms(2)) ' skipped']);
            continue;
       elseif any(ms0.labels{atoms(1)}~='H') %second atom is H
           HBatoms = [ms0.btB(find(ms0.btA==atoms(2))) ms0.btA(find(ms0.btB==atoms(2))) atoms(2) atoms(1)];%#ok
       elseif any(ms0.labels{atoms(2)}~='H') %first atom is H
           HBatoms = [ms0.btB(find(ms0.btA==atoms(1))) ms0.btA(find(ms0.btB==atoms(1))) atoms(1) atoms(2)];%#ok
       else %both atoms are H
           disp(['dihydrogen bond with atoms ' int2str(atoms(1)) ',' int2str(atoms(2)) ' skipped']);
           continue;
       end
       if ~isempty(HBatoms) && numel(HBatoms)==3
           if HBatoms(3) < HBatoms(1)
               HBatoms = [HBatoms(3) HBatoms(2) HBatoms(1)];
           end
           allHBind(end+1,:) = HBatoms;%#ok
       elseif numel(HBatoms)==2
           disp(['Two atoms contact found (' int2str(HBatoms) '), ircdb ind=' int2str(uniqCP_irc_struct_ind(icp)) ', skipped.' ]);
       end
    end %icp=1:size(uniqCPind,1)
    uniqHBind = unique(allHBind,'rows'); %unique AH...B bonds
    

%-----------------------------plotting------------------------------------------  
    
    [x,sort_ind]=sort([ircdb.irc]);
    min_x=min(x);
    max_x=max(x);
    xlabel_str='IRC, a.m.u.^{1/2}*bohr';
    
%-----------------------------Dipole moment plotting------------------------------------------  
    y=repmat(NaN,size(x));

    for iirc=1:numel(sort_ind) % over IRC points
        y(iirc) = ircdb(sort_ind(iirc)).DM;%#ok
    end
    txt = ['Dipole monent, Debay' ];%#ok
    data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;

    cla
    pointsign='.-';
    color='green';
    plot( a_tmp, x,y,pointsign,'Color',color,'LineWidth',lw);    
    my_axis([min_x max_x min(y) max(y)]);
    grid off
    set(a_tmp,'Box','on','XMinorTick','on','YMinorTick','on');
    xlabel(a_tmp,xlabel_str);
    ylabel(a_tmp,'Dipole moment, D');
    print(f_tmp,'-dpsc2', '-append', '-r300', psfile);

            
  
%------------------------Distancies, angles, charges figure-----------------------------------------------  
    
    if ~isempty(uniqHBind)

        %---------------------------- A...B --------------------------------
        uniqAB=unique([uniqHBind(:,1) uniqHBind(:,3)],'rows');
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};
        cla
        hold on
        for iHB=1:size(uniqAB,1) 
            row = uniqAB(iHB,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                y(iirc) = adist(ircdb(sort_ind(iirc)), row(1), row(2));%#ok
            end
            leg=[atomlabel(ircdb(1),row(1)) ',' atomlabel(ircdb(1),row(2))];%#ok
            txt = ['A...B (' leg '), A' ];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            pointsign='.-';
            color=pcolor{mod(iHB-1,numel(pcolor))+1};
            buf.h_plots(end+1)=plot( a_tmp, x,y, pointsign, 'Color',color, 'LineWidth',lw);    
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
        end
        my_axis([min_x max_x min_y max_y]);
        grid off
        set(a_tmp,'Box','on','XMinorTick','on','YMinorTick','on');
        xlabel(a_tmp,xlabel_str);
        ylabel(a_tmp,'dAB, A');
        legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);

        %---------------------------- A...H, H...B ------------------------
        uniqAH=unique([uniqHBind(:,1:2); uniqHBind(:,2:3)],'rows');
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};
        cla
        hold on
        for iHB=1:size(uniqAH,1) 
            row = uniqAH(iHB,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                y(iirc) = adist(ircdb(sort_ind(iirc)), row(1), row(2));%#ok
            end
            leg=[atomlabel(ircdb(1),row(1)) ',' atomlabel(ircdb(1),row(2))];%#ok
            txt = ['A...H (' leg '), A' ];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            pointsign='.-';
            color=pcolor{mod(iHB-1,numel(pcolor))+1};
            buf.h_plots(end+1)=plot( a_tmp, x,y, pointsign, 'Color',color, 'LineWidth',lw);    
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
        end
        my_axis([min_x max_x min_y max_y]);
        grid off
        set(a_tmp,'Box','on','XMinorTick','on','YMinorTick','on');
        xlabel(a_tmp,xlabel_str);
        ylabel(a_tmp,'dAH/HB, A');
        legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);

        %---------------------------- A-H-B -------------------------------
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};
        cla
        hold on
        for iHB=1:size(uniqHBind,1) %over all H-bonds
            row = uniqHBind(iHB,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                y(iirc) = valang(ircdb(sort_ind(iirc)), row(1), row(2), row(3));%#ok
            end
            leg=[atomlabel(ircdb(1),row(1)) ',' atomlabel(ircdb(1),row(2))...
                    ',' atomlabel(ircdb(1),row(3))];%#ok
            txt = ['A-H-B (' leg '), °' ];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            pointsign='.-';
            color=pcolor{mod(iHB-1,numel(pcolor))+1};
            buf.h_plots(end+1)=plot( a_tmp, x,y, pointsign, 'Color',color, 'LineWidth',lw);    
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
        end
        my_axis([min_x max_x min_y max_y]);
        grid off
        set(a_tmp,'Box','on','XMinorTick','on','YMinorTick','on');
        xlabel(a_tmp,xlabel_str);
        ylabel(a_tmp,'\angle AH…B, °');
        legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);
    
        %---------------------------- charges on atoms - hydrogens
        uniqatoms=unique(uniqHBind);
        uniqatoms_H = uniqatoms(strcmpcellar(ircdb(1).labels(uniqatoms),'H'));
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};
        cla
        hold on
        for iHB=1:size(uniqatoms_H,1) 
            row = uniqatoms_H(iHB,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                if isfield(ircdb(sort_ind(iirc)),'mcharge') && ~isempty(ircdb(sort_ind(iirc)).mcharge) && ~isempty( ircdb(sort_ind(iirc)).mcharge(row(1)) )
                    y(iirc) = ircdb(sort_ind(iirc)).mcharge(row(1));%#ok
                end
            end
            leg=[atomlabel(ircdb(1),row(1))];%#ok
            txt = ['Mulliken q ' leg ', a.u.' ];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            pointsign='.-';
            color=pcolor{mod(iHB-1,numel(pcolor))+1};
            buf.h_plots(end+1)=plot( a_tmp, x,y, pointsign, 'Color',color, 'LineWidth',lw);    
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
        end
        my_axis([min_x max_x min_y max_y]);
        grid off
        set(a_tmp,'Box','on','XMinorTick','on','YMinorTick','on');
        xlabel(a_tmp,xlabel_str);
        ylabel(a_tmp,'Mulliken atomic charge, e');
        legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);

        %---------------------------- charges on atoms - non-hydrogens
        uniqatoms=unique(uniqHBind);
        uniqatoms_notH = setdiff(uniqatoms,uniqatoms_H);
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};
        cla
        hold on
        for iHB=1:size(uniqatoms_notH,1) 
            row = uniqatoms_notH(iHB,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                if isfield(ircdb(sort_ind(iirc)),'mcharge') && ~isempty(ircdb(sort_ind(iirc)).mcharge) && ~isempty( ircdb(sort_ind(iirc)).mcharge(row(1)) )
                    y(iirc) = ircdb(sort_ind(iirc)).mcharge(row(1));%#ok
                end
            end
            leg=[atomlabel(ircdb(1),row(1))];%#ok
            txt = ['Mulliken q ' leg ', a.u.' ];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            pointsign='.-';
            color=pcolor{mod(iHB-1,numel(pcolor))+1};
            buf.h_plots(end+1)=plot( a_tmp, x,y, pointsign, 'Color',color, 'LineWidth',lw);    
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
        end
        my_axis([min_x max_x min_y max_y]);
        grid off
        set(a_tmp,'Box','on','XMinorTick','on','YMinorTick','on');
        xlabel(a_tmp,xlabel_str);
        ylabel(a_tmp,'Mulliken atomic charge, e');
        legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);

            
        %---------------------------- additional distances
        if exist('bondlist','var')
            min_y=NaN; max_y=NaN;
            buf={};buf.h_plots=[];buf.legs={};
            cla
            hold on
            for iHB=1:size(bondlist,1) %over all additional distances
                row = bondlist(iHB,:);
                y=repmat(NaN,size(x));
                for iirc=1:numel(sort_ind) % over IRC points
                    y(iirc) = adist(ircdb(sort_ind(iirc)), row(1), row(2));%#ok
                end
                leg=[atomlabel(ircdb(1),row(1)) ',' atomlabel(ircdb(1),row(2))];%#ok
                txt = ['H...H (' leg '), A' ];
                data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;

                pointsign='.-';
                color=pcolor{mod(iHB-1,numel(pcolor))+1};
                buf.h_plots(end+1)=plot( a_tmp, x,y, pointsign, 'Color',color, 'LineWidth',lw);    
                buf.legs{end+1}=leg;%#ok
                min_y = min([min_y min(y)]);
                max_y = max([max_y max(y)]);
            end
            my_axis([min_x max_x min_y max_y]);
            grid off
            set(a_tmp,'Box','on','XMinorTick','on','YMinorTick','on');
            xlabel(a_tmp,xlabel_str);
            ylabel(a_tmp,'R(H-H), A');
            legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
            print(f_tmp,'-dpsc2', '-append', '-r300', psfile);
        end

        %---------------------------- additional angles
        if exist('anglist','var')
            min_y=NaN; max_y=NaN;
            buf={};buf.h_plots=[];buf.legs={};
            cla
            hold on
            for iang=1:size(anglist,1) %over all additional angles
                row = anglist(iang,:);
                y=repmat(NaN,size(x));
                for iirc=1:numel(sort_ind) % over IRC points
                    y(iirc) = valang(ircdb(sort_ind(iirc)), row(1), row(2), row(3));%#ok
                end
                leg=[atomlabel(ircdb(1),row(1)) ',' atomlabel(ircdb(1),row(2)) ',' atomlabel(ircdb(1),row(3))];%#ok
                txt = ['glyc angle (' leg '), °' ];
                data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;

                pointsign='.-';
                color=pcolor{mod(iang-1,numel(pcolor))+1};
                buf.h_plots(end+1)=plot( a_tmp, x,y, pointsign, 'Color',color, 'LineWidth',lw);    
                buf.legs{end+1}=leg;%#ok
                min_y = min([min_y min(y)]);
                max_y = max([max_y max(y)]);
            end
            my_axis([min_x max_x min_y max_y]);
            grid off
            set(a_tmp,'Box','on','XMinorTick','on','YMinorTick','on');
            xlabel(a_tmp,xlabel_str);
            ylabel(a_tmp,'Glycosidic angle, °');
            legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
            print(f_tmp,'-dpsc2', '-append', '-r300', psfile);
        end

    end %if ~isempty(uniqHBind)

    
    
%------------------------AIM figures--------------------------  
    if exist('aimbondlist','var')
        uniqCPind = [uniqCPind; aimbondlist]; % #ok
    end
    if ~isempty(uniqCPind)

%------------------------rho --------------------------  
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};
        cla
        hold on
        for iCP=1:size(uniqCPind,1) % over all H-bond CPs
            row = uniqatoms(iHB,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                ircrec=ircdb(sort_ind(iirc));
                if ~isempty(ircrec.AIM)
                    iii = find( sum( repmat(uniqCPind(iCP,:),size(ircrec.AIM.atoms,1),1)==ircrec.AIM.atoms, 2) == 2 ); %find appropriate CP index among molecule CPs
                    if iii
                        y(iirc) = ircrec.AIM.ro(iii);%#ok
                    end
                end
            end
            leg=[atomlabel(ircdb(1),uniqCPind(iCP,1)) ',' atomlabel(ircdb(1),uniqCPind(iCP,2))];%#ok
            txt = ['\rho,a.u. (' leg ')'];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            pointsign='.-';
            color=pcolor{mod(iCP-1,numel(pcolor))+1};
            buf.h_plots(end+1)=plot( a_tmp, x,y, pointsign, 'Color',color, 'LineWidth',lw);    
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
        end
        my_axis([min_x max_x min_y max_y]);
        grid off
        set(a_tmp,'Box','on','XMinorTick','on','YMinorTick','on');
        xlabel(a_tmp,xlabel_str);
        ylabel(a_tmp,'\rho, a.u.');
        legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);        
        

%------------------------ delta rho --------------------------  
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};
        cla
        hold on
        for iCP=1:size(uniqCPind,1) % over all H-bond CPs
            row = uniqatoms(iHB,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                ircrec=ircdb(sort_ind(iirc));
                if ~isempty(ircrec.AIM)
                    iii = find( sum( repmat(uniqCPind(iCP,:),size(ircrec.AIM.atoms,1),1)==ircrec.AIM.atoms, 2) == 2 ); %find appropriate CP index among molecule CPs
                    if iii
                        y(iirc) = ircrec.AIM.DelSqRho(iii);%#ok
                    end
                end
            end
            leg=[atomlabel(ircdb(1),uniqCPind(iCP,1)) ',' atomlabel(ircdb(1),uniqCPind(iCP,2))];%#ok
            txt = ['\Delta \rho,a.u. (' leg ')'];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            pointsign='.-';
            color=pcolor{mod(iCP-1,numel(pcolor))+1};
            buf.h_plots(end+1)=plot( a_tmp, x,y, pointsign, 'Color',color, 'LineWidth',lw);    
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
        end
        my_axis([min_x max_x min_y max_y]);
        grid off
        set(a_tmp,'Box','on','XMinorTick','on','YMinorTick','on');
        xlabel(a_tmp,xlabel_str);
        ylabel(a_tmp,'\Delta \rho, a.u.');
        legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);        
        
        %------------------------ ellipticity --------------------------  
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};
        cla
        hold on
        for iCP=1:size(uniqCPind,1) % over all H-bond CPs
            row = uniqatoms(iHB,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                ircrec=ircdb(sort_ind(iirc));
                if ~isempty(ircrec.AIM)
                    iii = find( sum( repmat(uniqCPind(iCP,:),size(ircrec.AIM.atoms,1),1)==ircrec.AIM.atoms, 2) == 2 ); %find appropriate CP index among molecule CPs
                    if iii
                        y(iirc) = ircrec.AIM.BondEl(iii);%#ok
                    end
                end
            end
            leg=[atomlabel(ircdb(1),uniqCPind(iCP,1)) ',' atomlabel(ircdb(1),uniqCPind(iCP,2))];%#ok
            txt = ['\epsilon (' leg ')'];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            pointsign='.-';
            color=pcolor{mod(iCP-1,numel(pcolor))+1};
            buf.h_plots(end+1)=plot( a_tmp, x,y, pointsign, 'Color',color, 'LineWidth',lw);    
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
        end
        my_axis([min_x max_x min_y max_y]);
        if isfield(flags,'semilogy_for_ellipticity') && flags.semilogy_for_ellipticity
            set(a_tmp,'YScale','log');
        end
        grid off
        set(a_tmp,'Box','on','XMinorTick','on','YMinorTick','on');
        xlabel(a_tmp,xlabel_str);
        ylabel(a_tmp,'ellipticity');
        legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);        

        %------------------------ E_{HB} --------------------------  
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};
        cla
        hold on
        for iCP=1:size(uniqCPind,1) % over all H-bond CPs
            row = uniqatoms(iHB,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                ircrec=ircdb(sort_ind(iirc));
                if ~isempty(ircrec.AIM)
                    iii = find( sum( repmat(uniqCPind(iCP,:),size(ircrec.AIM.atoms,1),1)==ircrec.AIM.atoms, 2) == 2 ); %find appropriate CP index among molecule CPs
                    if iii
                        y(iirc) = -0.5*CC.encoef*ircrec.AIM.V(iii);%#ok
                    end
                end
            end
            leg=[atomlabel(ircdb(1),uniqCPind(iCP,1)) ',' atomlabel(ircdb(1),uniqCPind(iCP,2))];%#ok
            txt = ['E_{HB},k/m (' leg ')'];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            pointsign='.-';
            color=pcolor{mod(iCP-1,numel(pcolor))+1};
            buf.h_plots(end+1)=plot( a_tmp, x,y, pointsign, 'Color',color, 'LineWidth',lw);    
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
        end
        max_y = min([max_y 100]); %cutting off covalent bonds energy range
        my_axis([min_x max_x min_y max_y]);
        grid off
        set(a_tmp,'Box','on','XMinorTick','on','YMinorTick','on');
        xlabel(a_tmp,xlabel_str);
        ylabel(a_tmp,'E_{HB}, kcal/mol');
        legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);        

        
    else 
        fl_empty_figures = fl_empty_figures+1;
    end %if ~isempty(uniqCPind)

    
    if fl_empty_figures==0 && fl_savepics
%        saveas(gcf,[indir filesep fnameshort '.pdf'],'pdf');
%        plot2svg(svgfile, gcf);
        saveas(gcf,figfile,'fig');
    end
    
    if fl_tosave
        if exist(workdbname,'file')
            copyfile(workdbname,workdbnamebkp);
        end
        save(workdbname,'ircdb');
        
      try
        Excel = actxserver ('Excel.Application');
        File=xlsfile;
        if ~exist(File,'file')
          ExcelWorkbook = Excel.workbooks.Add;
          ExcelWorkbook.SaveAs(File,1);
          ExcelWorkbook.Close(false);
        end
        invoke(Excel.Workbooks,'Open',File);

        params={};
        xlswrite1spec(File,data,'IRCdata','A1',params);
        
        invoke(Excel.ActiveWorkbook,'Save');
        Excel.Quit
        Excel.delete
        clear Excel
      catch ME1
        invoke(Excel.ActiveWorkbook,'Save');
        Excel.Quit
        Excel.delete
        clear Excel
      end
        
    end


%break;

end %l_ind=1:numfiles

diary off
