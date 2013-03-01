%%
%rot40_irc: analysis of IRC 
%
% Version 1.1
% Last modified  R O Zhurakivsky 2011-10-09
% Created        R O Zhurakivsky 2011-03-21

%2012-1010 TBD: 
%+0) rename atoms to standart names
%1) separate CH...B and other bonds dependencies
%+2) add cutoffs
%+3) add molecule shapshots generation for cutoffs
%+4) add pictures scale alignment fot different computation methods

%2012-1012
%+на графіках енергій Н-звязків за Еспінозою треба залишати лише ті значення, для яких лаплас ро більше нуля,

%for snapshots add H-bonds plotting and their distances
%add output of H-bonds properties in TS and stationary states

%+1. позначення координат особливих точок налазять місцями один на одний і по різному округлені (має бути до сотих)
%+2. на рисунку відносної електронної енергії відсутні позначення початкових та кінцевих структур, що ускладнює аналіз графіка
%?3. на рисунках, на яких побудовано декілька графіків, слід біля кожного із них робити підпис, якому звязку він відповідає; а також слід вибрати різні значення, щоб у чорно-білому варіанті графік був читабельним
%?4. зроби окрулення так, як це було у надсіланих раніше приклада для пари АТ (особливо це стосується ро та лаплас ро)
%+5. Мюлікеновські заряди поки не будемо робити, тому не витрачай на них
%часу
%+/-6. для звязку С2Н...О2 дані наведи на окремих графіках

%2012-1015
% Зауваження до графіків (стосуються усіх графіків, що ти мені надіслав останнього разу)
% 1. графік для відносної електронної енергії по осі ординат має починатися від 0, не заходячи у відємну область.
% 2. нумерація ключових точок та легенди налазять на самі графіки. 
% 3. графік для ро не треба у масштабі X10-3, треба в нативному з троьма значеннями після коми (це ж стосується і лаплас ро)!!!!
% 4. для  пари HypThy структури 9 точок вийшли перевернутими догори. Виправ взявши за приклад структури як у статті в JBSD. Підписи під графіками мають  відповідати підписам як  у надісланому шаблоні. Зверни увагу, на різне позначення (різну густину точок) Н-звязків у цих структурах. Також значення довжин для середнього звязку зливаються.
% 5. уникай на одному графіку близьких значків  (н-д, ромб та кружочок, так як у чоронобілому форматі вони будуть виглядати практично однаково)
% 6.не подобаються коми на графіках для HypHyp - краще ввести розділення як N1H, HN1, HO6,O6H - так буде зрозуміліше 
% 7. для HypHyp не треба графіка дипольного моменту та глікозидних кутів 
% 8. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Для усіх графіків: значення ро, лаплас ро мають бути з трома знаками після коми, епс - з двома!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

format compact
global flplot

clear 
atomsind

% PLOTTYPE_deltaE=0;
% PLOTTYPE_mu=1;
% PLOTTYPE_Mulliken_q_H=2;
% PLOTTYPE_Mulliken_q_A=3;
% PLOTTYPE_rho=4;
% PLOTTYPE_deltarho=5;
% PLOTTYPE_epsilon=6;
% PLOTTYPE_EHB=7;
% PLOTTYPE_dAB=8;
% PLOTTYPE_dAH=9;
% PLOTTYPE_aAHB=10;
% PLOTTYPE_RHH=11;
% PLOTTYPE_aglyc=12;
% PLOTTYPE_rho_CH=13;
% PLOTTYPE_deltarho_CH=14;
% PLOTTYPE_epsilon_CH=15;
% PLOTTYPE_EHB_CH=16;
% PLOTTYPE_dAH_CH=17;
% PLOTTYPE_max=18;


%!!!2012-1030: MP2 enerhy detected for all point except IRC=0 is wrong - HF energy is
%taken instead of MP2

ATOMS=35;
Hbond_DelSqRho_crit_value = -0.026 %#ok<NOPTS> %if DelSqRho of (3,-1) CP is less than this value assume this CP as Hbond
%-0.022 is taken from Hyp_Hyp eps==1 deltarho(O6H)==0 : deltarho(O6H) in
%nearest calculated point irc=-0.316 is -0.021
%-0.025 is taken fron Hyp-Cyt eps=1 deltarho(O6H)==0: irc=0.031, nearest
%irc=0.060, deltarho=-0.025

global flags
flags={};
iCPcutoff=[0 0 0 0];
cutoffs=[];

aimbondlist=[];
bondlist=[];
anglist=[];

%-------------------------------------------------------------------
flags.develmode = 1;
flags.irc_tick_interval = NaN; %automatically

flplot=1    %#ok
fl_tosave = 1 %#ok
fl_plotpics = 1 %#ok
fl_savepics = 1 %#ok
fl_recreate_matfile = 0 %#ok
fl_reload_sp = 0   %#ok
fl_reload_freq = 0 %#ok
fl_reload_extout = 0 %#ok

%fl_plot_milliken = 0 %#ok

molpind=cell(1,ATOMS);
for i=1:numel(molpind)
    molpind{i}='';
end

%complextype='A_taut-A'; methods=[{'b3lyp_vac'}];
complextype='Hyp_taut-Hyp'; methods=[{'b3lyp_eps4'}];
%complextype='Hyp_taut-Hyp'; methods=[{'b3lyp_vac'}];
%complextype='AC'; methods=[{'mp2_vac'}];
%complextype='AC'; methods=[{'b3lyp_eps4'}];
%complextype='AC'; methods=[{'b3lyp_vac'}];
%complextype='GT'; methods=[{'b3lyp_vac'}];
%complextype='GT'; methods=[{'b3lyp_eps4'}];

%complextype='GC'; methods=[{'b3lyp_vac'},{'b3lyp_eps4'},{'mp2_vac'}];
%complextype='GC'; methods=[{'mp2_vac'}];
%complextype='Hyp-dimer'; methods=[{'b3lyp_vac'},{'b3lyp_eps4'}];
%complextype='Hyp-dimer'; methods=[{'b3lyp_vac'}];
%complextype='HypHyp'; methods=[{'b3lyp_vac'},{'b3lyp_eps4'}];
%complextype='HypCyt'; methods=[{'b3lyp_vac'},{'b3lyp_eps4'}];
%complextype='HypThy'; methods=[{'b3lyp_vac'},{'b3lyp_eps4'}];
%complextype='AT'; methods=[{'b3lyp_vac'},{'b3lyp_eps4'},{'mp2_vac'}];

common_out_dir = 'E:\work\Brovarets\plots';
if ~(exist(common_out_dir,'dir')==7)
     mkdir(common_out_dir);
end

disp(['Time: ' datestr(now)]);

for l_ind=1:numel(methods)
    cur_method=methods{l_ind};

    disp(['Time: ' datestr(now)]);

    if strcmp(complextype,'HypHyp')
        bondlist = [13 25]; %R_HH
        anglist  = [ 3 13 25; 13 25 24]; %alpha1, alpha2
        limits_filename = 'E:\work\Brovarets\120216_IRC_Hyp\Hyp-Hyp\Hyp-Hyp_limits.mat';
        IRCdesc = [{'(Hyp \bullet Hyp)'},{'(TS_{Hyp \bullet Hyp \leftrightarrow Hyp* \bullet Hyp*})'},{'(Hyp* \bullet Hyp*)'}];
        molpind{7}='C6';
        molpind{8}='O6';
        molpind{9}='N1';
        molpind{11}='C2';
        molpind{10}='H';
        molpind{14}='H';
        molpind{16}='N1';%''
        molpind{17}='C2';%''
        molpind{27}='C6';%''
        molpind{28}='O6';%''
        molpind{15}='H';%''
        molpind{18}='H';%''
        molpind{3}='N9';
        molpind{13}='H9';
        molpind{24}='N9';%''
        molpind{25}='H9';%''
        if strcmp(cur_method,'b3lyp_vac')
            indir = 'E:\work\Brovarets\120216_IRC_Hyp\Hyp-Hyp\Hyp-Hyp_tight_int' %#ok<NOPTS>
            workname='irc_b3lyp_tight_int'%#ok
    %        iCPcutoff=[1 3 0 0];
            iCPcutoff=[1 2 0 0];
            cutoffs=[ -6.64 -0.34 -0.15 0.00 0.20 3.37]; %interpolated
        elseif strcmp(cur_method,'b3lyp_eps4')
            indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp-Hyp\Hyp-Hyp_eps4_tight_int' %#ok<NOPTS>
            workname='irc_b3lyp_tight_eps4'%#ok
    %        iCPcutoff=[1 3 0 0];
            iCPcutoff=[1 2 0 0];
            cutoffs=[-6.92 -0.37 -0.19 0.00 0.17 3.17]; %interpolated
        elseif 0
            indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp-dimer'  %#ok
            workname='irc_b3lyp_631gdp'%#ok
            limits_filename = 'E:\work\Brovarets\120216_IRC_Hyp\Hyp-dimer\Hyp-dimer_limits.mat';
            flags.semilogy_for_ellipticity=1;
    %        iCPcutoff=[1 3 2 4];
            iCPcutoff=1:4;
        elseif 0
            indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp-dimer\Hyp-dimer_tight_int' %#ok<NOPTS>
            workname='irc_b3lyp_tight_int'%#ok
            limits_filename = 'E:\work\Brovarets\120216_IRC_Hyp\Hyp-dimer\Hyp-dimer_limits.mat';
            flags.semilogy_for_ellipticity=1;
    %        iCPcutoff=[1 3 2 4];
        elseif 1
            indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp-dimer\Hyp-dimer_tight_int_eps4' %#ok<NOPTS>
            workname='irc_b3lyp_tight_eps4'%#ok
            limits_filename = 'E:\work\Brovarets\120216_IRC_Hyp\Hyp-dimer\Hyp-dimer_limits.mat';
            flags.semilogy_for_ellipticity=1;
            iCPcutoff=[1 2 3 4];
            cutoffs = [-18.87 -0.73 -0.51 -0.25  0.00  4.27  4.57  4.70 23.42];
        end

    elseif strcmp(complextype,'HypCyt')
        bondlist = [12 24]; %R_HH
        anglist  = [ 1 12 24; 12 24 13]; %alpha1, alpha2
        IRCdesc = [{'(Hyp \bullet Cyt)'},{'(TS_{Hyp \bullet Cyt \leftrightarrow Hyp* \bullet Cyt*})'},{'(Hyp* \bullet Cyt*)'}];
        limits_filename = 'E:\work\Brovarets\120216_IRC_Hyp\Hyp-Cyt\Hyp-Cyt_limits.mat';
        molpind{6}='C6';
        molpind{7}='O6';
        molpind{8}='N1';
        molpind{9}='C2';
        molpind{19}='N4';
        molpind{21}='N3';
        molpind{23}='O2';
        molpind{25}='H';
        molpind{26}='H';
        molpind{27}='H';
        molpind{1}='N9';
        molpind{12}='H9';
        molpind{13}='N1';
        molpind{24}='H1';
        if strcmp(cur_method,'b3lyp_vac')
            indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp-Cyt\Hyp-Cyt_tight_int' %#ok<NOPTS>
            workname='irc_b3lyp_tight_int'%#ok
    %        iCPcutoff=[1 3 2 4];
            iCPcutoff=1:4;
            cutoffs=[-5.15 -1.01 -0.79 -0.51 0.00 0.03 0.19 0.44 3.23]; %interpolated
        elseif strcmp(cur_method,'b3lyp_eps4')
    %        indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp-Cyt\Hyp-Cyt_tight_int_eps4'
            indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp-Cyt\Hyp-Cyt_tight_int_eps4_full' %#ok<NOPTS>
    %        iCPcutoff=[1 4 2 3 ];
            iCPcutoff=1:4;
            cutoffs=[-11.15 -7.12 -6.91 -6.65 -0.20 -0.06  0.00  0.20  1.82]; %interpolated
    %        indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp-Cyt\Hyp-Cyt_tight_int_eps4_p2_corrected'
            workname='irc_b3lyp_tight_eps4'%#ok
            flags.irc_tick_interval = 1;
        end

    elseif strcmp(complextype,'HypThy')
        bondlist = [13 21]; %R_HH
        anglist  = [ 1 13 21; 13 21 20]; %alpha1, alpha2
        IRCdesc = [{'(Hyp\ast \bullet Thy)'},{'(TS_{Hyp* \bullet Thy \leftrightarrow Hyp \bullet Thy*})'},{'(Hyp \bullet Thy*)'}];
        molpind{6}='C6';
        molpind{7}='O6';
        molpind{8}='N1';
        molpind{10}='C2';
        molpind{14}='O4';
        molpind{16}='C4';
        molpind{17}='N3';
        molpind{18}='C2';
        molpind{19}='O2';
        molpind{9}='H';
        molpind{15}='H';
        molpind{29}='H';
        molpind{1}='N9';
        molpind{13}='H9';
        molpind{20}='N1';
        molpind{21}='H1';
        if strcmp(cur_method,'b3lyp_vac')
            indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp-Thy\Hyp-Thy_tight_int' %#ok<NOPTS>
            workname='irc_b3lyp_tight_int'%#ok
            limits_filename = 'E:\work\Brovarets\120216_IRC_Hyp\Hyp-Thy\Hyp-Thy_limits.mat';
    %        iCPcutoff=[1 4 2 3 ];
            iCPcutoff=1:4;
            cutoffs=[-5.94 -0.50 -0.26 0.00 0.01 0.14 0.34 0.51 4.82]; %interpolated
        elseif strcmp(cur_method,'b3lyp_eps4')
            indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp-Thy\Hyp-Thy_tight_int_eps4' %#ok<NOPTS>
            workname='irc_b3lyp_tight_eps4'%#ok
            limits_filename = 'E:\work\Brovarets\120216_IRC_Hyp\Hyp-Thy\Hyp-Thy_limits.mat';
    %        iCPcutoff=[1 4 2 3 ];
            iCPcutoff=1:4;
            cutoffs=[-6.37 -0.40 -0.17 0.00 0.03 0.27 0.45 0.62 5.52]; %interpolated
        elseif 0
            indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp_O_Thy\Hyp_O_Thy_tight_int' %#ok<NOPTS>
            workname='irc_b3lyp_tight'%#ok
            limits_filename = 'E:\work\Brovarets\120216_IRC_Hyp\Hyp_O_Thy\Hyp_O_Thy_limits.mat';
    %        iCPcutoff=[1 4 2 3 ];
            iCPcutoff=1:4;
    %        cutoffs=[-5.52 -0.28 -0.45 -0.62 0.00 -0.06 0.17 0.39 6.37];
        elseif 0
            indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp_O_Thy\Hyp_O_Thy_tight_int_eps4' %#ok<NOPTS>
            workname='irc_b3lyp_tight_eps4'%#ok
            limits_filename = 'E:\work\Brovarets\120216_IRC_Hyp\Hyp_O_Thy\Hyp_O_Thy_limits.mat';
    %        iCPcutoff=[1 4 2 3 ];
            iCPcutoff=1:4;
    %        cutoffs=[-5.52 -0.28 -0.45 -0.62 0.00 -0.06 0.17 0.39 6.37];
        end

    elseif strcmp(complextype,'Hyp-dimer')
        bondlist = [13 25]; %R_HH
        anglist  = [ 3 13 25; 13 25 24]; %alpha1, alpha2
        limits_filename = 'E:\work\Brovarets\120216_IRC_Hyp\Hyp-dimer\Hyp-dimer_limits.mat';
        IRCdesc = [{'(Hyp \bullet Hyp)'},{'(TS_{Hyp \bullet Hyp \leftrightarrow Hyp* \bullet Hyp*})'},{'(Hyp* \bullet Hyp*)'}];
        molpind{7}='C6';
        molpind{8}='O6';
        molpind{9}='N1';
        molpind{11}='C2';
        molpind{10}='H';
        molpind{14}='H';
        molpind{16}='N1';%''
        molpind{17}='C2';%''
        molpind{27}='C6';%''
        molpind{28}='O6';%''
        molpind{15}='H';%''
        molpind{18}='H';%''
        molpind{3}='N9';
        molpind{13}='H9';
        molpind{24}='N9';%''
        molpind{25}='H9';%''
        if strcmp(cur_method,'b3lyp_vac_nontight')
            indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp-dimer\Hyp-dimer'  %#ok
            workname='irc_b3lyp'%#ok
            limits_filename = 'E:\work\Brovarets\120216_IRC_Hyp\Hyp-dimer\Hyp-dimer\Hyp-dimer_limits.mat';
            flags.semilogy_for_ellipticity=1;
    %        iCPcutoff=[1 3 2 4];
            iCPcutoff=1:4;
            cutoffs=[-34.41 -28.58 -2.36  0.00  0.23 34.72];
        elseif strcmp(cur_method,'b3lyp_vac')
            indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp-dimer\Hyp-dimer_tight_int' %#ok<NOPTS>
            workname='irc_b3lyp_tight_int'%#ok
            limits_filename = 'E:\work\Brovarets\120216_IRC_Hyp\Hyp-dimer\Hyp-dimer_limits.mat';
            flags.semilogy_for_ellipticity=1;
            iCPcutoff=[1 2 4 6];
            cutoffs=[-36.63 -2.79 -2.67 -2.36  0.00  0.24  0.51  0.73 25.70];
        elseif strcmp(cur_method,'b3lyp_eps4')
            indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp-dimer\Hyp-dimer_tight_int_eps4' %#ok<NOPTS>
            workname='irc_b3lyp_tight_eps4'%#ok
            flags.semilogy_for_ellipticity=1;
            iCPcutoff=[1 2 3 4];
            cutoffs = [-18.87 -0.73 -0.51 -0.25  0.00  4.27  4.57  4.70 23.42];
        end

    elseif strcmp(complextype,'GC')
        bondlist = [14 26];
        anglist  = [ 1 14 26; 15 26 14];
        IRCdesc = [{'(G \bullet C)'},{'(TS_{G \bullet C \leftrightarrow G* \bullet C*})'},{'(G* \bullet C*)'}];
        limits_filename = 'E:\work\Brovarets\120411_irc_AT_GC\GC_limits.mat';
        molpind{6}='C6';
        molpind{7}='O6';
        molpind{8}='N1';
        molpind{9}='C2';
        molpind{10}='N2';
        molpind{21}='N4';
        molpind{23}='N3';
        molpind{24}='C2';
        molpind{25}='O2';
        molpind{27}='H';
        molpind{28}='H';
        molpind{29}='H';
        molpind{1}='N9';
        molpind{14}='H9';
        molpind{15}='N1';
        molpind{26}='H1';
        if strcmp(cur_method,'b3lyp_vac')
            indir = 'E:\work\Brovarets\120411_irc_AT_GC\b3lyp\GC_tight_int' %#ok<NOPTS>
            workname ='irc_b3lyp_tight_int'%#ok
            iCPcutoff=[1 2 3 4 ];
            cutoffs = [-6.46 -0.42 -0.15 0.00 0.00 0.03 0.28 0.50 8.60];
        elseif strcmp(cur_method,'b3lyp_eps4')
            indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\GC_eps4_tight_int' %#ok<NOPTS>
            workname ='irc_b3lyp_tight_eps4'%#ok
            iCPcutoff=[1 2 3 4 ];
            cutoffs = [-10.11 -0.46 -0.25 -0.00  0.00  5.90  6.04  6.31  9.29];
        elseif strcmp(cur_method,'mp2_vac')
            indir='E:\work\Brovarets\120411_irc_AT_GC\mp2\GC_MP2' %#ok<NOPTS>
            workname = 'irc120411_mp2' %#ok<NOPTS>
    %        workname='irc120425_mp2_tight'%#ok
            iCPcutoff=[1 2 3 4 ];
            cutoffs = [-6.91 -0.37 -0.09 -0.09 0.00 0.12 0.22 0.46 8.49];
        end

    elseif strcmp(complextype,'AT')
        bondlist = [11 24]; %indexes of edge glicosidic atoms
        anglist  = [ 1 11 24; 15 24 11];
    %    aimbondlist = [8 14];
        IRCdesc = [{'(A \bullet T)'},{'(TS_{A \bullet T \leftrightarrow A* \bullet T*})'},{'(A* \bullet T*)'}];
        limits_filename = 'E:\work\Brovarets\120411_irc_AT_GC\AT_limits.mat';
        molpind{6}='N6';
        molpind{7}='N1';
        molpind{8}='C2';
        molpind{20}='O4';
        molpind{19}='C4';
        molpind{21}='N3';
        molpind{22}='C2';
        molpind{23}='O2';
        molpind{14}='H';
        molpind{29}='H';
        molpind{30}='H';
        molpind{1}='N9';
        molpind{11}='H9';
        molpind{15}='N1';
        molpind{24}='H1';
        if strcmp(cur_method,'b3lyp_vac')
            indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\AT_tight_int_step5';
            workname = 'irc_b3lyp_tight_int' %#ok<NOPTS>
        elseif strcmp(cur_method,'b3lyp_eps4')
            indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\AT_eps4_tight_int';
            workname ='irc_b3lyp_tight_eps4'%#ok
        elseif strcmp(cur_method,'mp2_vac')
            indir='E:\work\Brovarets\120411_irc_AT_GC\mp2\AT_MP2';
            workname = 'irc120411_mp2' %#ok<NOPTS>
        end

    elseif strcmp(complextype,'AC')
        bondlist = [15 24]; %indexes of edge glicosidic atoms
        anglist  = [ 1 15 24; 23 24 15];
    %    aimbondlist = [8 14];
        IRCdesc = [{'(A \bullet C)'},{'(TS_{A \bullet C \leftrightarrow A* \bullet C*})'},{'(A* \bullet C*)'}];
        limits_filename = 'E:\work\Brovarets\120411_irc_AT_GC\AC_limits.mat';
         molpind{6}='C6';
         molpind{7}='O6';
         molpind{10}='N1';
         molpind{11}='C2';
         molpind{17}='N4';
         molpind{19}='C4';
         molpind{20}='N3';
         molpind{21}='C2';
         molpind{22}='O2';
         molpind{9}='H_{I}';
         molpind{16}='H_{II}';
         molpind{12}='H_{III}';
         molpind{1}='N9';
         molpind{15}='H9';
         molpind{23}='N1';
         molpind{24}='H1';
        if strcmp(cur_method,'b3lyp_vac')
            indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\AC_tight_int';
            workname = 'irc_b3lyp_tight_int' %#ok<NOPTS>
            iCPcutoff=[1 2 3 4];
            cutoffs=[-5.85 -1.48 -1.26 -0.98 -0.09 0.00 0.16 0.37 4.07];
        elseif strcmp(cur_method,'b3lyp_eps4')
            indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\AC_eps4_tight_int';
            workname ='irc_b3lyp_tight_eps4'%#ok
            iCPcutoff=[1 2 3 4];
            cutoffs=[-9.23 -4.07 -3.86 -3.60 -0.08 0.00 0.15 0.37 4.10];
        elseif strcmp(cur_method,'mp2_vac')
             indir='E:\work\Brovarets\120411_irc_AT_GC\mp2\AC_mp2';
             workname = 'irc120411_mp2' %#ok<NOPTS>
            iCPcutoff=[1 2 3 4];
            cutoffs=[-6.02 -0.99 -0.76 -0.46 -0.11 0.00 0.13 0.35 4.71];
        end

    elseif strcmp(complextype,'GT')
        bondlist = [16 24];
        anglist  = [ 1 16 24; 23 24 16];
        IRCdesc = [{'(G \bullet T)'},{'(TS_{G \bullet T \leftrightarrow G* \bullet T*})'},{'(G* \bullet T*)'}];
        limits_filename = 'E:\work\Brovarets\120411_irc_AT_GC\GT_limits.mat';
         molpind{6}='C6';
         molpind{7}='O6';
         molpind{8}='N1';
         molpind{10}='C2';
         molpind{11}='N2';
         molpind{17}='O4';
         molpind{19}='C4';
         molpind{20}='N3';
         molpind{21}='C2';
         molpind{22}='O2';
         molpind{18}='H_{I}';
         molpind{9}='H_{II}';
         molpind{12}='H_{III}';
         molpind{1}='N9';
         molpind{16}='H9';
         molpind{23}='N1';
         molpind{24}='H1';
        if strcmp(cur_method,'b3lyp_vac')
            indir = 'E:\work\Brovarets\120411_irc_AT_GC\b3lyp\GT_tight_int_1' %#ok<NOPTS>
            workname ='irc_b3lyp_tight_int'%#ok
            iCPcutoff=[1 2 4 5];
            cutoffs = [-5.36 -1.97 -1.80 -1.61 -0.18 0.00 0.06 0.29 4.35];
         elseif strcmp(cur_method,'b3lyp_eps4')
             indir='E:\work\Brovarets\120411_irc_AT_GC\b3lyp\GT_eps4_tight_int' %#ok<NOPTS>
             workname ='irc_b3lyp_tight_eps4'%#ok
%             iCPcutoff=[1 2 3 4 ];
%             cutoffs = [-10.11 -0.46 -0.25 -0.00  0.00  5.90  6.04  6.31  9.29];
%         elseif strcmp(cur_method,'mp2_vac')
%             indir='E:\work\Brovarets\120411_irc_AT_GC\mp2\GC_MP2' %#ok<NOPTS>
%             workname = 'irc120411_mp2' %#ok<NOPTS>
%     %        workname='irc120425_mp2_tight'%#ok
%             iCPcutoff=[1 2 3 4 ];
%             cutoffs = [-6.91 -0.37 -0.09 -0.09 0.00 0.12 0.22 0.46 8.49];
        end

    elseif strcmp(complextype,'Hyp_taut-Hyp')
        bondlist = [13 25]; %R_HH
        anglist  = [ 3 13 25; 13 25 24]; %alpha1, alpha2
        limits_filename = 'E:\work\Brovarets\120216_IRC_Hyp\Hyp_taut-Hyp\TS_Hyp_taut-Hyp.mat';
        IRCdesc = [{'(Hyp* \bullet Hyp)'},{'(TS_{Hyp* \bullet Hyp \leftrightarrow Hyp \bullet Hyp*})'},{'(Hyp \bullet Hyp*)'}];
         molpind{7}='C6';
         molpind{8}='O6';
         molpind{9}='N1';
         molpind{11}='C2';
         molpind{10}='H';
         molpind{14}='H';
         molpind{16}='N1''';%''
         molpind{17}='C2''';%''
         molpind{27}='C6''';%''
         molpind{28}='O6''';%''
         molpind{15}='H''';%''
         molpind{18}='H''';%''
         molpind{3}='N9';
         molpind{13}='H9';
         molpind{24}='N9''';%''
         molpind{25}='H9''';%''
        if strcmp(cur_method,'b3lyp_vac')
            indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp_taut-Hyp\b3lyp_tight_int' %#ok<NOPTS>
            workname='irc_b3lyp_tight_int'%#ok
            iCPcutoff=[1 4 2 3];
            cutoffs=[-5.61 -0.28 -0.04 0.00 0.04 0.28 5.87];
            borders=[-5.61 5.61];
%          elseif strcmp(cur_method,'b3lyp_eps4')
%              indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp_taut-Hyp\b3lyp_eps4' %#ok<NOPTS>
%              workname='irc_b3lyp_tight_eps4'%#ok
%              iCPcutoff=[1 2 3 4];
%              cutoffs = [-5.72 -0.35 -0.12 0.00 0.10 0.14 0.33 0.50 5.98];
         elseif strcmp(cur_method,'b3lyp_eps4')
             indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp_taut-Hyp\b3lyp_eps4_mirror' %#ok<NOPTS>
             workname='irc_b3lyp_tight_eps4'%#ok
%             iCPcutoff=[1 2 3 4];
%             cutoffs = [-5.72 -0.35 -0.12 0.00 0.10 0.14 0.33 0.50 5.98];
        end
    elseif strcmp(complextype,'A_taut-A')
        bondlist = [13 25]; %R_HH
        anglist  = [ 3 13 25; 13 25 24]; %alpha1, alpha2
        limits_filename = 'E:\work\Brovarets\1302_irc_AA\TS_A_taut-A.mat';
        IRCdesc = [{'(A* \bullet A)'},{'(TS_{A* \bullet A \leftrightarrow A \bullet A*})'},{'(A \bullet A*)'}];
if 0
        molpind{7}='C6';
         molpind{8}='O6';
         molpind{9}='N1';
         molpind{11}='C2';
         molpind{10}='H';
         molpind{14}='H';
         molpind{16}='N1''';%''
         molpind{17}='C2''';%''
         molpind{27}='C6''';%''
         molpind{28}='O6''';%''
         molpind{15}='H''';%''
         molpind{18}='H''';%''
         molpind{3}='N9';
         molpind{13}='H9';
         molpind{24}='N9''';%''
         molpind{25}='H9''';%''
end         
        if strcmp(cur_method,'b3lyp_vac')
            indir='E:\work\Brovarets\1302_irc_AA\b3lyp_vac' %#ok<NOPTS>
            workname='irc_b3lyp_tight_int'%#ok
%             iCPcutoff=[1 4 2 3];
%             cutoffs=[-5.61 -0.28 -0.04 0.00 0.04 0.28 5.87];
%         elseif strcmp(cur_method,'b3lyp_eps4')
%             indir='E:\work\Brovarets\120216_IRC_Hyp\Hyp-dimer\Hyp-dimer_tight_int_eps4' %#ok<NOPTS>
%             workname='irc_b3lyp_tight_eps4'%#ok
%             flags.semilogy_for_ellipticity=1;
%             iCPcutoff=[1 2 3 4];
%             cutoffs = [-18.87 -0.73 -0.51 -0.25  0.00  4.27  4.57  4.70 23.42];
        end
        
    end
%-------------------------------------------------------------------

diaryfname0=[indir filesep 'logfile'];
diaryfname=diaryfname0;
for i=2:10000
  if ~(exist(diaryfname,'file')==2), break, end
  diaryfname = [diaryfname0 int2str(i)];
end
diary(diaryfname)

indir %#ok
diaryfname %#ok

gtemplname=[workname '_templ.gjf']  %#ok
fullgtemplname=[CD.templatesdir filesep gtemplname] %#ok

fl_empty_figures = 0;

odir=indir;
if exist(odir,'dir')~=7
   mkdir(odir);
end
dlm=strfind(odir,filesep);
baseodir = odir(dlm(end)+1:end);
gjf_odir = [baseodir '_gjf'];

 lfiles = dir(strcat(indir,filesep,'*.list'));
 num_listfiles = size(lfiles,1);
 if ~num_listfiles
   error('No LIST files to analyse found');
 end
% 
% for l_ind=1:num_listfiles
    listingfile=lfiles(1); %processing only one listing per directory!!!
    

    %forming filenames on the base on list filename
    dlm=strfind(listingfile.name,'.');
    fnameshort = listingfile.name(1:(dlm-1));
    workdbname = [indir filesep fnameshort '.mat']; 
    dlm=strfind(workdbname,'.');
    workdbnamebkp=[workdbname(1:dlm(end)-1) '~' workdbname(dlm(end):end)];
    xlsfile = [indir filesep fnameshort '.xls'] %#ok
    
    psfile = [indir filesep fnameshort '.ps'];
    psfilebkp=[psfile(1:dlm(end)-1) '~' psfile(dlm(end):end)];
    svgfile = [indir filesep fnameshort '.svg'];
    figfile = [indir filesep fnameshort '.fig'];

    
    if fl_recreate_matfile || exist(workdbname,'file')~=2

        disp(['Reading ' listingfile.name ]);
        [outfiles,irc_directions,irc_shifts,irc_types] =
                textread([indir '\' listingfile.name],'%s%d%f%d');
        %irc_types - 1 - refinement, 2 - calculate as mirror
        IRC_TYPE_REFINE = 1;
        IRC_TYPE_MIRROR = 2;
        
        if ~numel(outfiles)
            disp('No Gaussian output files found. Skipping.');
            continue
        end
        ircdb={};

        for out_ind=1:numel(outfiles)
            if (outfiles{out_ind}(1)==';' || outfiles{out_ind}(1)=='%'), continue, end %skip commented lines

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
            irc_isrefinement = 0;
            irc_is_mirror = 0;
            if out_ind<=numel(irc_types)
                if irc_isrefinements(out_ind)==IRC_TYPE_REFINE
                    irc_isrefinement = 1;
                elseif irc_isrefinements(out_ind)==IRC_TYPE_MIRROR
                    irc_is_mirror = 1;
                end
            end

        %    try
                disp(['Loading file: ' out.ffname])
                %processing single OUT file
                fid=fopen(out.ffname,'r'); % don't use text mode - this is dangerous for Unix files ;)
                if fid<0
                    error(['error while opening file ' out.ffname])
                end

                fl_ircstart = 0;
                fl_good_geom_found = 0;
                fl_chkdata_found = 0;

                frewind(fid);

                clear ms0
                while 1 
                    tline = fgetl(fid);
                    if ~ischar(tline), break, end

                    if ~isempty(strfind(tline,'Cartesian coordinates read from the checkpoint file')) 
                        disp('Initial coordinates found'); %point #0
                        fl_chkdata_found = 1;

                        for i=1:2, tline=fgets(fid); end  
                        [ms0,status]=extrgeom(fid,1,0);
                        if status
                            disp([worktitle,': ',lastwarn]);
                            fclose(fid);
                            break
                        end
                    
                        ms0.point = 0;
                        ms0.irc = 0 + irc_shift;

                    elseif ~isempty(strfind(tline,'Redundant internal coordinates taken from checkpoint file'))
                        disp('Initial coordinates found'); %point #0
                        fl_chkdata_found = 1;

                        for i=1:2, tline=fgets(fid); end 
                        
                        [ms0,status]=extrgeom(fid,1,0);
                        if status
                            disp([worktitle,': ',lastwarn]);
                            fclose(fid);
                            break
                        end
                    
                        ms0.point = 0;
                        ms0.irc = 0 + irc_shift;
                    elseif ~isempty(strfind(tline,'INPUT DATA FOR L123'))
                        if ~fl_chkdata_found 
                            disp(['Pure IRC w/o checkpoint data detected']); %#ok<NBRAK>
                            fl_ircstart = 1;
                            break;
                        end

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
                        if exist([odir filesep gjf_odir filesep ms0.desc '.gjf'],'file')~=2 %create input Gaussian file
                            if exist([odir filesep gjf_odir],'dir')~=7
                                mkdir(odir,gjf_odir);
                            end
                            savemolgs([odir filesep gjf_odir],ms0,3,order,fullgtemplname); %Gaussian with XYZ
                        end
                        
                        ms0.irc_isrefinement = irc_isrefinement;
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
                       [v1,v2]=strread(tline,'%s%f','delimiter','=');
                       ms0.energy = v2;
                    elseif ~isempty(strfind(tline,'EUMP2 =')) 
                    % E2 =    -0.2992511380D+01 EUMP2 =    -0.91922697437759D+03
                       [v1,v2,v3,v4]=strread(tline,'%s%f%s%f','delimiter','=');
                       ms0.energy = v4; %rewrote with MP2 energy if found

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

                        if ms0.point==0
                            ms0.irc = irc_shift;
                        else
                            %NET REACTION COORDINATE UP TO THIS POINT =    0.26904
                            pat = '\s*NET\sREACTION\sCOORDINATE\sUP\D+([\d\.]+)';
                            A = regexp(buf, pat, 'tokens','once');
                            A=[A{:}];
                            ms0.irc = str2double(A{1})*irc_direction + irc_shift;
                        end
                        
                        if fl_good_geom_found
                            ms0.desc = [out.fnameshort '_#' sprintf('%03d',ms0.point)];
                            order=1:ms0.atomnum;
                            if exist([odir filesep gjf_odir filesep ms0.desc '.gjf'],'file')~=2 %create input Gaussian file
                                if exist([odir filesep gjf_odir],'dir')~=7
                                    mkdir(odir,gjf_odir);
                                end
                                savemolgs([odir filesep gjf_odir],ms0,3,order,fullgtemplname); %Gaussian with XYZ
                            end

                            ms0.irc_isrefinement = irc_isrefinement;
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
    bdwidth = 2; %figure border (points)
    topbdwidth = 20; %figure top border (points)
    set(0,'Units','pixels') 
%    scnsize = get(0,'ScreenSize');
    scnsize = [0 0 1024 576];
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
    f_tmp = figure('PaperOrientation','landscape',...
           'PaperType','A5',...
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
    label_fontsize = 14;

    data = {[]};
       
%----------------------Energy and dipole moment----------------------------------------------------

    figure(f_tmp);
    hold on
    if exist(limits_filename,'file')==2
        load(limits_filename,'PlLim');
    else
        PlLim=PlotsLimits; %new instance of class PlotsLimits
    end
    
    [x,sort_ind]=sort([ircdb.irc]);
    if flags.develmode || ~isfield(ircdb,'irc_isrefinement')
        points2plot = ones(1,numel(ircdb));
    else
        points2plot = ~[ircdb.irc_isrefinement];
        points2plot = points2plot(sort_ind);
    end
    y=[ircdb.energy];
    y=y(sort_ind);
    y=(y-min(y))*CC.encoef;
    
    col = 1;
    data(1,col)={'ind'};                data(2:numel(ircdb)+1,col) = num2cell(sort_ind); col=col+1;
    data(1,col)={'filename'};           data(2:numel(ircdb)+1,col) = {ircdb(sort_ind).desc}; col=col+1;
    data(1,col)={'IRC, Bohr'}; data(2:numel(ircdb)+1,col) = num2cell(x); col=col+1;

    fig_desc = 'E';
    txt = 'energy, kcal/mol';
    data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
    
    if fl_plotpics
        h_plot = irc_plot(gca,x(points2plot==1),y(points2plot==1),1); %#ok<NASGU>
        params={};
        params.flags=flags;
        params.limits = PlLim.correctlimit(fig_desc,[min(x) max(x) 0 max(y)]);
        my_axis(params,fig_desc,cutoffs,l_ind,[1 1 0 1],1);
        fig_name = [complextype '_' fig_desc '_' cur_method];
        set(gcf, 'Name',fig_name, 'FileName',fig_name);
        orient landscape
        export_fig([common_out_dir filesep fig_name],'-eps','-pdf','-transparent', '-r300');
    %    print(f_tmp,'-dpsc2', '-append', '-r300', psfile);
    end
    
    %----------------differentiating E over IRC
    fig_desc = 'diffE_irc';
    txt = '{\delta}E/{\delta}IRC'; 

    xd = x(points2plot==1); %#ok<FNDSB> %bug with x(points2plot)
    yd = y(points2plot==1); %#ok<FNDSB> %bug with x(points2plot)
    [xd, xd_ind] = unique(xd);
    yd = yd(xd_ind);
    y_deriv = deriv(yd)./deriv(xd);

    data(1,col)={'IRC, Bohr'};   
    data(2:numel(y_deriv)+1,col) = num2cell(xd); col=col+1;
    data(1,col)={txt};   
    data(2:numel(y_deriv)+1,col) = num2cell(y_deriv); col=col+1;
    
    if fl_plotpics
        cla%reset(gca)
        h_plot = irc_plot(gca,xd,y_deriv,1);
        params={};
        params.flags=flags;
        params.limits = PlLim.correctlimit(fig_desc,[min(xd) max(xd) min(y_deriv) max(y_deriv)]);
        my_axis(params,fig_desc,cutoffs,l_ind);
        fig_name = [complextype '_' fig_desc '_' cur_method];
        set(gcf,'Name',fig_name,'FileName',fig_name);
        orient landscape
        export_fig([common_out_dir filesep fig_name],'-eps','-pdf','-transparent', '-r300');
    end
    
%    print(f_tmp,'-dpsc2', '-append', '-r300', psfile);
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
                  
              elseif ~isempty(strfind(tline,'Summary of Natural Population Analysis:'))
                  for ii=1:5, tline = fgetl(fid); end
                  
                  A=fscanf(fid,'%s %d %f %f %f %f %f',[7,double(ircdb(iext).atomnum)]);
                  if isempty(A)
                    warning('error in output file - couldn''n load charges. skipping');%#ok
                    break
                  end
                  ircdb(iext).('nbocharge')=A(3,:);%#ok %NBO Natural charges
                  
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

               %exclude covalent/ionic bonds
%                aa = ms0.btB(ms0.btA==atoms(1));%#ok %atoms connected to atoms(1) 
%                if sum(aa==atoms(2))>0 
               if CPs(i).DelSqRho > Hbond_DelSqRho_crit_value %-0.02
                   AIM.is_hbond(end+1,:)=1;
                   if CPs(i).DelSqRho < 0
                       warning('CAN:irc',['contant between atoms ' int2str(atoms(1)) ' and ' int2str(atoms(2)) ' assumed as H-bond, but DelSqRho='...
                           num2str(CPs(i).DelSqRho,'%0.3f') ]);
                   end
               else
                   AIM.is_hbond(end+1,:)=0;
               end
               
            end
            
            ircdb(iext).AIM = AIM;%#ok
            
        else
            disp(['!!! can''t found: ' 'wfn' filesep mol.desc '.extout'])
        end %if exist(aimfile,'file')==2

    end %for iext=1:numel(ircdb)
%    end %if ~isfield(ircdb,'AIM')

    %creating array with indexes of atoms around 'H-bond like BCPs' for all structures
    %along IRC
    uniqCPind = []; %contains pairs of indexes of bond atoms
    uniqCP_Hatom_ind = []; %what H atom is in each H-bond from uniqCPind
    uniqCP_irc_struct_ind = [];
    uniqCP_Aatom_ind = [];
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
            if isempty(uniqCPind) || ~any(sum(uniqCPind==repmat(atoms,size(uniqCPind,1),1),2)==2) %add if is new one
                uniqCPind(end+1,:) = atoms;%#ok
                
                ii=find([mol.labels{atoms}]=='H',1);
                curHatom = atoms(ii);
                curBatom = atoms(3-ii);
                uniqCP_Hatom_ind(end+1) = curHatom; %#ok<AGROW>
                
                %find A bond atom
                thirdatom=[];
                for jCP=1:size(mol.AIM.atoms,1)
                    if mol.AIM.atoms(jCP,1)==curHatom && mol.AIM.atoms(jCP,2)~=curBatom
                        thirdatom = mol.AIM.atoms(jCP,2);
                        break;
                    elseif mol.AIM.atoms(jCP,2)==curHatom && mol.AIM.atoms(jCP,1)~=curBatom
                        thirdatom = mol.AIM.atoms(jCP,1);
                        break
                    end
                end
                if isempty(thirdatom)
                    warning('CAN:irc',['third atom for CP(' int2str(curHatom) '...' int2str(curbatom) ') is not founf in ind_irc ' int2str(ind_irc) ]);
                end
                uniqCP_Aatom_ind(end+1) = thirdatom;%#ok
                
                uniqCP_irc_struct_ind(end+1) = ind_irc;%#ok
            end
        end
    end %for ind_irc=1:numel(ircdb)

    % sorting uniqCPind so that bonds with same H atom be neighbours
    [xx,sort_ind]=sort(uniqCP_Hatom_ind); 
    uniqCPind = uniqCPind(sort_ind,:);
    uniqCP_Hatom_ind = uniqCP_Hatom_ind(sort_ind);
    uniqCP_irc_struct_ind = uniqCP_irc_struct_ind(sort_ind);
    uniqCP_Aatom_ind = uniqCP_Aatom_ind(sort_ind);
    
    
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
           HBatoms = [uniqCP_Aatom_ind(icp) atoms(2) atoms(1)];
       elseif any(ms0.labels{atoms(2)}~='H') %first atom is H
           HBatoms = [uniqCP_Aatom_ind(icp) atoms(1) atoms(2)];
       else %both atoms are H
           disp(['dihydrogen bond with atoms ' int2str(atoms(1)) ',' int2str(atoms(2)) ' skipped']);
           continue;
       end
       if ~isempty(HBatoms) && numel(HBatoms)==3
           if HBatoms(3) < HBatoms(1) %sorting
               HBatoms = [HBatoms(3) HBatoms(2) HBatoms(1)];
           end
           allHBind(end+1,:) = HBatoms;%#ok
       elseif numel(HBatoms)==2
           disp(['Two atoms contact found (' int2str(HBatoms) '), ircdb ind=' int2str(uniqCP_irc_struct_ind(icp)) ', skipped.' ]);
       end
    end %icp=1:size(uniqCPind,1)
    
    uniqHBind = unique(allHBind,'rows'); %unique AH...B bonds
    if ~isempty(uniqHBind)
        iCHB = [strcmpcellar(ircdb(1).labels(uniqHBind(:,1)),'C') strcmpcellar(ircdb(1).labels(uniqHBind(:,3)),'C')]; %indexes of CHB bonds
        CHBbondHind=uniqHBind(iCHB,2); %indexes of hydrogen atoms involved in CH..B bonds
    else
        CHBbondHind=[];
    end

%-----------------------------plotting------------------------------------------  
    
    [x,sort_ind]=sort([ircdb.irc]);
    min_x=min(x);
    max_x=max(x);
    
  

%-----------------------------Dipole moment plotting------------------------------------------  
    y=repmat(NaN,size(x));

    for iirc=1:numel(sort_ind) % over IRC points
        y(iirc) = ircdb(sort_ind(iirc)).DM;%#ok
    end
    if ~all(isinf(y)|isnan(y))
        fig_desc = 'mu';
        txt = ['Dipole moment, Debay' ];%#ok
        data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;

        if fl_plotpics
            cla%reset(gca)
            y(y==Inf)=NaN;
            h_plot = irc_plot( a_tmp, x(points2plot==1),y(points2plot==1),1);    
            params={};
            params.flags=flags;
            params.limits = PlLim.correctlimit(fig_desc,[min_x max_x min(y) max(y)]);
            my_axis(params,fig_desc,cutoffs);
            fig_name = [complextype '_' fig_desc '_' cur_method];
            set(gcf,'Name',fig_name,'FileName',fig_name);
            orient landscape
            export_fig([common_out_dir filesep fig_name],'-eps','-pdf','-transparent', '-r300');
%        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);
        end
    end
            
  
%------------------------Distancies, angles, charges figure-----------------------------------------------  
    
    if fl_plotpics && ~isempty(uniqHBind)

        AHB_ind=zeros(0,3); %non CHB bonds
        CHB_ind=zeros(0,3);
        for i=1:size(uniqHBind,1)
            if ~isempty(uniqHBind) && (ircdb(1).labels{uniqHBind(i,1)}=='C' || ircdb(1).labels{uniqHBind(i,3)}=='C')
                CHB_ind(end+1,:)=uniqHBind(i,:); %#ok<AGROW>
            else
                AHB_ind(end+1,:)=uniqHBind(i,:); %#ok<AGROW>
            end
        end

for icycle=1:2 %2 loops - one for AH...B and one for CH...B bonds

        if icycle==1         
            HB_ind = AHB_ind;
        else
            HB_ind = CHB_ind;
        end %icycle
        if isempty(HB_ind), continue, end        
        
        %---------------------------- A...B --------------------------------
        uniqAB=unique([HB_ind(:,1) HB_ind(:,3)],'rows');
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};buf.xpos=[];buf.ypos=[];buf.ht=[];
        cla%reset(gca)
        hold on
        for iHB=1:size(uniqAB,1) 
            row = uniqAB(iHB,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                y(iirc) = adist(ircdb(sort_ind(iirc)), row(1), row(2));
            end
%            leg=[atomlabel(ircdb(1),row(1)) ',' atomlabel(ircdb(1),row(2))];%#ok
            leg=[molpind{row(1)} '...' molpind{row(2)}];
            txt = ['A...B (' leg '), A' ];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            if iHB==1
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1),1);
            else
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1));    
            end
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
            if flags.develmode || size(uniqAB,1)>1
                buf.ht(end+1) = text(x(10),min(max_y,y(10)+0.05*(max_y-min_y)),leg,'FontSize',label_fontsize);
            end
        end
        if icycle==1, fig_desc='dAB'; else fig_desc='dAB_CHB'; end
        params={};
        params.flags=flags;
        params.limits = PlLim.correctlimit(fig_desc,[min_x max_x min_y max_y]);
        my_axis(params,fig_desc,cutoffs,l_ind);
        movetexts(buf,flags);
        
        if flags.develmode
            legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        end

        fig_name = [complextype '_' fig_desc '_' cur_method];
        set(gcf,'Name',fig_name,'FileName',fig_name);
        orient landscape
        export_fig([common_out_dir filesep fig_name],'-eps','-pdf','-transparent', '-r300');
%        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);

        %---------------------------- A...H, H...B ------------------------
        uniqAH=unique([HB_ind(:,1:2); HB_ind(:,2:3)],'rows');
        numofCHbonds = numel(strcmpcellar(ircdb(1).labels(uniqAH),'C')); %hack
        
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};buf.ht=[];
        cla%reset(gca)
        hold on
        for iHB=1:size(uniqAH,1) 
            row = uniqAH(iHB,:);
            if ircdb(1).labels{row(1)}=='C' || ircdb(1).labels{row(2)}=='C' %skipping CH bonds plotting
                continue
            end
            
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                y(iirc) = adist(ircdb(sort_ind(iirc)), row(1), row(2));
            end
%            leg=[atomlabel(ircdb(1),row(1)) ',' atomlabel(ircdb(1),row(2))];%#ok
            if ircdb(1).labels{row(1)}=='H'
                leg=[molpind{row(2)} '...' molpind{row(1)}];
            else
                leg=[molpind{row(1)} '...' molpind{row(2)}];
            end
            txt = ['A...H (' leg '), A' ];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            if iHB==1
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1),1);    
            else
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1));    
            end
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
            if flags.develmode || size(uniqAH,1)-numofCHbonds>1
                buf.ht(end+1) = text(x(10),min(max_y,y(10)+0.05*(max_y-min_y)),leg,'FontSize',label_fontsize);
            end
        end
        if icycle==1, fig_desc='dHB'; else fig_desc='dHB_CHB'; end
        params={};
        params.flags=flags;
        params.limits = PlLim.correctlimit(fig_desc,[min_x max_x min_y max_y]);
        my_axis(params,fig_desc,cutoffs,l_ind);
        movetexts(buf,flags);
        if flags.develmode
            legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        end
        
        fig_name = [complextype '_' fig_desc '_' cur_method];
        set(gcf,'Name',fig_name,'FileName',fig_name);
        orient landscape
        export_fig([common_out_dir filesep fig_name],'-eps','-pdf','-transparent', '-r300');
%        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);

        %---------------------------- A-H-B -------------------------------
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};buf.ht=[];
        cla%reset(gca)
        hold on
        for iHB=1:size(HB_ind,1) %over all H-bonds
            row = HB_ind(iHB,:);
%            if row(1)==row(3), continue, end
            
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                y(iirc) = valang(ircdb(sort_ind(iirc)), row(1), row(2), row(3));
            end
%            leg=[atomlabel(ircdb(1),row(1)) ',' atomlabel(ircdb(1),row(2))...
%                    ',' atomlabel(ircdb(1),row(3))];%#ok
            leg=[molpind{row(1)} '' molpind{row(2)} '...' molpind{row(3)}];
            txt = ['A-H-B (' leg '), °' ];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            if iHB==1
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1),1);    
            else
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1));    
            end
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
            if flags.develmode || size(HB_ind,1)>1
                buf.ht(end+1) = text(x(10),min(max_y,y(10)+0.05*(max_y-min_y)),leg,'FontSize',label_fontsize);
            end
       end
        if icycle==1, fig_desc='aAHB'; else fig_desc='aAHB_CHB'; end
        params={};
        params.flags=flags;
        params.limits = PlLim.correctlimit(fig_desc,[min_x max_x min_y max_y]);
        my_axis(params,fig_desc,cutoffs,l_ind);
        movetexts(buf,flags);
        if flags.develmode
            legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        end
        
        fig_name = [complextype '_' fig_desc '_' cur_method];
        set(gcf,'Name',fig_name,'FileName',fig_name);
        orient landscape
        export_fig([common_out_dir filesep fig_name],'-eps','-pdf','-transparent', '-r300');
%        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);

end %for icycle=1:2

        %---------------------------- Milliken charges on atoms - all
        uniqatoms=unique(uniqHBind);
%        uniqatoms_H = uniqatoms(strcmpcellar(ircdb(1).labels(uniqatoms),'H'));
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};buf.ht=[];
        cla%reset(gca)
        hold on
        for iHB=1:size(uniqatoms,1) 
            row = uniqatoms(iHB,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                if isfield(ircdb(sort_ind(iirc)),'mcharge') && ~isempty(ircdb(sort_ind(iirc)).mcharge) && ~isempty( ircdb(sort_ind(iirc)).mcharge(row(1)) )
                    y(iirc) = ircdb(sort_ind(iirc)).mcharge(row(1));%#ok
                end
            end
%            leg=[atomlabel(ircdb(1),row(1))];%#ok
            leg=[molpind{row(1)}];
            txt = ['Mulliken q ' leg ', e' ];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            if iHB==1
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1),1);    
            else
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1));    
            end
            buf.legs{end+1}=leg;%#ok
%            ht = text(x(10),max(y,y(10)+0.05*(max(y)-min(y))),leg,'FontSize',label_fontsize);
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
            buf.ht(end+1)= text(x(10),min(max_y,y(10)+0.05*(max_y-min_y)),leg,'FontSize',label_fontsize);
        end
        fig_desc = 'Milliken';
        params={};
        params.flags=flags;
        params.limits = PlLim.correctlimit(fig_desc,[min_x max_x min_y max_y]);
        my_axis(params,fig_desc,cutoffs,l_ind);
        movetexts(buf,flags);
        if flags.develmode
            legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        end

        fig_name = [complextype '_' fig_desc '_' cur_method];
        set(gcf,'Name',fig_name,'FileName',fig_name);
        orient landscape
        export_fig([common_out_dir filesep fig_name],'-eps','-pdf','-transparent', '-r300');
%        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);

        %---------------------------- NBO charges on atoms - hydrogens
if 0        
        uniqatoms=unique(uniqHBind);
        uniqatoms_H = uniqatoms(strcmpcellar(ircdb(1).labels(uniqatoms),'H'));
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};buf.ht=[];
        cla%reset(gca)
        hold on
        for iHB=1:size(uniqatoms_H,1) 
            row = uniqatoms_H(iHB,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                if isfield(ircdb(sort_ind(iirc)),'nbocharge') && ~isempty(ircdb(sort_ind(iirc)).nbocharge) && ~isempty( ircdb(sort_ind(iirc)).nbocharge(row(1)) )
                    y(iirc) = ircdb(sort_ind(iirc)).nbocharge(row(1));%#ok
                end
            end
%            leg=[atomlabel(ircdb(1),row(1))];%#ok
            leg=[molpind{row(1)}];
            txt = ['NBO q ' leg ', e' ];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            if iHB==1
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1),1);    
            else
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1));    
            end
            buf.legs{end+1}=leg;%#ok
%            text(x(10),max(y,y(10)+0.05*(max(y)-min(y))),leg,'FontSize',label_fontsize);
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
            buf.ht(end+1) = text(x(10),min(max_y,y(10)+0.05*(max_y-min_y)),leg,'FontSize',label_fontsize);
        end
        fig_desc = 'NBO_q';
        params={};
        params.flags=flags;
        params.limits = PlLim.correctlimit(fig_desc,[min_x max_x min_y max_y]);
        my_axis(params,fig_desc,cutoffs,l_ind);
        movetexts(buf,flags);
        if flags.develmode
            legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        end

        fig_name = [complextype '_' fig_desc '_' cur_method];
        set(gcf,'Name',fig_name,'FileName',fig_name);
        orient landscape
        export_fig([common_out_dir filesep fig_name],'-eps','-pdf','-transparent', '-r300');
%        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);
end

        %---------------------------- additional distances
        if exist('bondlist','var')
            min_y=NaN; max_y=NaN;
            buf={};buf.h_plots=[];buf.legs={};buf.ht=[];
            cla%reset(gca)
            hold on
            for iHB=1:size(bondlist,1) %over all additional distances
                row = bondlist(iHB,:);
                y=repmat(NaN,size(x));
                for iirc=1:numel(sort_ind) % over IRC points
                    y(iirc) = adist(ircdb(sort_ind(iirc)), row(1), row(2));%#ok
                end
%                leg=[atomlabel(ircdb(1),row(1)) ',' atomlabel(ircdb(1),row(2))];%#ok
                leg=[molpind{row(1)} '-' molpind{row(2)}];
                txt = ['H...H (' leg '), A' ];
                data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;

                if iHB==1
                    buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1),1);    
                else
                    buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1));
                end
                buf.legs{end+1}=leg;%#ok
                min_y = min([min_y min(y)]);
                max_y = max([max_y max(y)]);
                if flags.develmode
                    buf.ht(end+1) = text(x(10),min(max_y,y(10)+0.05*(max_y-min_y)),leg,'FontSize',label_fontsize);
                end
            end
            fig_desc = 'dHH';
            params={};
            params.flags=flags;
            params.limits = PlLim.correctlimit(fig_desc,[min_x max_x min_y max_y]);
            my_axis(params,fig_desc,cutoffs,l_ind);
            movetexts(buf,flags);
            if flags.develmode
                legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
            end
            
            fig_name = [complextype '_' fig_desc '_' cur_method];
            set(gcf,'Name',fig_name,'FileName',fig_name);
            orient landscape
            export_fig([common_out_dir filesep fig_name],'-eps','-pdf','-transparent', '-r300');
%            print(f_tmp,'-dpsc2', '-append', '-r300', psfile);
        end

        %---------------------------- additional angles
        if exist('anglist','var')
            min_y=NaN; max_y=NaN;
            buf={};buf.h_plots=[];buf.legs={};buf.ht=[];
            cla%reset(gca)
            hold on
            for iang=1:size(anglist,1) %over all additional angles
                row = anglist(iang,:);
                y=repmat(NaN,size(x));
                for iirc=1:numel(sort_ind) % over IRC points
                    y(iirc) = valang(ircdb(sort_ind(iirc)), row(1), row(2), row(3));%#ok
                end
%                leg=[atomlabel(ircdb(1),row(1)) ',' atomlabel(ircdb(1),row(2)) ',' atomlabel(ircdb(1),row(3))];%#ok
                if flags.develmode
                    leg=[molpind{row(1)} '-' molpind{row(2)} '-' molpind{row(3)}];
                else
                    leg=['{\alpha}_' int2str(iang) ];
                end
                txt = ['alpha (' leg '), degree' ];
                data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;

                if iang==1
                    buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1),1);
                else
                    buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1));
                end
                buf.legs{end+1}=leg;%#ok
%                text(x(10),max(y,y(10)+0.05*(max(y)-min(y))),leg,'FontSize',label_fontsize);
                min_y = min([min_y min(y)]);
                max_y = max([max_y max(y)]);
                buf.ht(end+1)= text(x(10),min(max_y,y(10)+0.05*(max_y-min_y)),leg,'FontSize',label_fontsize);
            end
            fig_desc = 'glycang';
            params={};
            params.flags=flags;
            params.limits = PlLim.correctlimit(fig_desc,[min_x max_x min_y max_y]);
            my_axis(params,fig_desc,cutoffs,l_ind);
            movetexts(buf,flags);
            if flags.develmode
                legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
            end
            
            fig_name = [complextype '_' fig_desc '_' cur_method];
            set(gcf,'Name',fig_name,'FileName',fig_name);
            orient landscape
            export_fig([common_out_dir filesep fig_name],'-eps','-pdf','-transparent', '-r300');
%            print(f_tmp,'-dpsc2', '-append', '-r300', psfile);
        end

    end %if ~isempty(uniqHBind)

    
    
%------------------------AIM figures--------------------------  
    if ~isempty(uniqCPind)
        if exist('aimbondlist','var')
            uniqCPind = [uniqCPind; aimbondlist]; %#ok<AGROW> % #ok
        end

        AHB_CPind=zeros(0,2); %non CHB bond CPs
        CHB_CPind=zeros(0,2); %CHB bond CPs
        for i=1:size(uniqCPind,1)
            if ~isempty(CHBbondHind) && sum(uniqCPind(i,:)==CHBbondHind) %errorneous !!! H atom may take place in different H-bonds
                CHB_CPind(end+1,:)=uniqCPind(i,:); %#ok<AGROW>
            else
                AHB_CPind(end+1,:)=uniqCPind(i,:); %#ok<AGROW>
            end
        end

        keypoints={};
        keypoints.irc=repmat(NaN,1,9); %IRC value of the point
        keypoints.inddb=repmat(NaN,1,9); %structure index in ircdb
        keypoints.desc={}; %key point description
        keypoints.atoms(9)={[]};
        
        keypoints.inddb(1)=sort_ind(x==min(x)); 
        keypoints.irc(1)=min([ircdb.irc]); 
        keypoints.desc(1)=IRCdesc(1);
        
        ind = find(x==0,1);
        if ~isempty(ind)
            keypoints.irc(5)=0;
            keypoints.inddb(5)=sort_ind(ind);
            keypoints.desc(5)=IRCdesc(2);
        end
        
        keypoints.irc(9)=max([ircdb.irc]);
        keypoints.inddb(9)=sort_ind(x==max(x));
        keypoints.desc(9)=IRCdesc(3);

for icycle=1:2 %2 loops - one for AH...B and one for CH...B bonds

        if icycle==1         
            CPs = AHB_CPind;
        else
            CPs = CHB_CPind;
        end %icycle
        if isempty(CPs), continue, end
%------------------------rho --------------------------  
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};buf.x=[];buf.y=[];buf.title={};buf.ht=[];buf.atoms=[];
        cla%reset(gca)
        hold on
        for iCP=1:size(CPs,1) % over all H-bond CPs
            row = CPs(iCP,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                ircrec=ircdb(sort_ind(iirc));
                if ~isempty(ircrec.AIM)
                    iii = find( sum( repmat(row,size(ircrec.AIM.atoms,1),1)==ircrec.AIM.atoms, 2) == 2 ); %find appropriate CP index among molecule CPs
                    if iii
                        y(iirc) = ircrec.AIM.ro(iii);%#ok
                    end
                end
            end
%            leg=[atomlabel(ircdb(1),uniqCPind(iCP,1)) ',' atomlabel(ircdb(1),uniqCPind(iCP,2))];%#ok
            if icycle==2, dlm='...'; else dlm=''; end
            if ircdb(1).labels{row(1)}=='H'
                leg=[molpind{row(2)} dlm molpind{row(1)}];%#ok
            else
                leg=[molpind{row(1)} dlm molpind{row(2)}];%#ok
            end
            txt = ['\rho,a.u. (' leg ')'];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;

            buf.x(end+1,:) = x;
            buf.y(end+1,:) = y;
            buf.title(end+1) = {txt};
            buf.atoms(end+1,:) = row;
            
            if iCP==1
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1),1);
            else
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1));
            end
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
            if flags.develmode || size(CPs,1)>1
                buf.ht(end+1) = text(x(10),min(max_y,y(10)+0.05*(max_y-min_y)),leg,'FontSize',label_fontsize);
            end
        end
        if icycle==1, fig_desc='rho'; else fig_desc='rhoCHB'; end
        params={};
        params.flags=flags;
        params.limits = PlLim.correctlimit(fig_desc,[min_x max_x min_y max_y]);
        my_axis(params,fig_desc,cutoffs,l_ind);
        movetexts(buf,flags);
        if flags.develmode
            legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        end
        
        fig_name = [complextype '_' fig_desc '_' cur_method];
        set(gcf,'Name',fig_name,'FileName',fig_name);
        orient landscape
        export_fig([common_out_dir filesep fig_name],'-eps','-pdf','-transparent', '-r300');
%        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);        
        
if icycle==1         
        if iCPcutoff(1) && iCPcutoff(2)
            [X0, i1] = findlinroot( buf.x(iCPcutoff(1),:), buf.y(iCPcutoff(1),:)-buf.y(iCPcutoff(2),:) ); %rho1-rho2 %buf.x is already sorted
disp(['point3: irc=' num2str(X0,'%.3f') ', nearest irc=' num2str(buf.x(iCPcutoff(1),i1),'%.3f') ', rho1-rho2=' num2str(buf.y(iCPcutoff(1),i1)-buf.y(iCPcutoff(2),i1),'%.3f') ]);
            keypoints.irc(3) = X0;
            keypoints.inddb(3) = sort_ind(i1);
            keypoints.atoms(3) = {[buf.atoms(iCPcutoff(1),:); buf.atoms(iCPcutoff(2),:)]};
            
            keypoints.desc{3}=['\rho_{' buf.legs{iCPcutoff(1)} '} = \rho_{' buf.legs{iCPcutoff(2)} '}'];
%            disp(['P3 rho(' buf.legs{iCPcutoff(1)} ')==rho(' buf.legs{iCPcutoff(2)} '): irc=' num2str(keypoints.irc(3),3) ])
        end
        if iCPcutoff(3) && iCPcutoff(4)
            [X0, i1] = findlinroot( buf.x(iCPcutoff(3),:), buf.y(iCPcutoff(3),:)-buf.y(iCPcutoff(4),:) ); %rho1-rho2 %buf.x is already sorted
disp(['point7: irc=' num2str(X0,'%.3f') ', nearest irc=' num2str(buf.x(iCPcutoff(3),i1),'%.3f') ', rho3-rho4=' num2str(buf.y(iCPcutoff(3),i1)-buf.y(iCPcutoff(4),i1),'%.3f') ]);
            keypoints.irc(7) = X0;
            keypoints.inddb(7)=sort_ind(i1);
            keypoints.atoms(7) = {[buf.atoms(iCPcutoff(3),:); buf.atoms(iCPcutoff(4),:)]};
            keypoints.desc{7}=['\rho_{' buf.legs{iCPcutoff(3)} '} = \rho_{' buf.legs{iCPcutoff(4)} '}'];
%            disp(['P3 rho(' buf.legs{iCPcutoff(1)} ')==rho(' buf.legs{iCPcutoff(2)} '): irc=' num2str(keypoints.irc(3),3) ])
        end
end %icycle       
%------------------------ delta rho --------------------------  
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};buf.x=[];buf.y=[];buf.title={};buf.ht=[];
        cla%reset(gca)
        hold on
        for iCP=1:size(CPs,1) % over all H-bond CPs
            row = CPs(iCP,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                ircrec=ircdb(sort_ind(iirc));
                if ~isempty(ircrec.AIM)
                    iii = find( sum( repmat(row,size(ircrec.AIM.atoms,1),1)==ircrec.AIM.atoms, 2) == 2 ); %find appropriate CP index among molecule CPs
                    if iii
                        y(iirc) = ircrec.AIM.DelSqRho(iii);%#ok
                    end
                end
            end
%            leg=[atomlabel(ircdb(1),uniqCPind(iCP,1)) ',' atomlabel(ircdb(1),uniqCPind(iCP,2))];%#ok
            if icycle==2, dlm='...'; else dlm=''; end
            if ircdb(1).labels{row(1)}=='H'
                leg=[molpind{row(2)} dlm molpind{row(1)}];%#ok
            else
                leg=[molpind{row(1)} dlm molpind{row(2)}];%#ok
            end
            txt = ['\Delta \rho,a.u. (' leg ')'];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            buf.x(end+1,:) = x;
            buf.y(end+1,:) = y;
            buf.title(end+1) = {txt};
            
            if iCP==1
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1),1);
            else
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1));
            end
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
            if flags.develmode || size(CPs,1)>1
                buf.ht(end+1) = text(x(10),min(max_y,y(10)+0.05*(max_y-min_y)),leg,'FontSize',label_fontsize);
            end
        end
        if icycle==1, fig_desc='deltarho'; else fig_desc='deltarhoCHB'; end
        params={};
        params.flags=flags;
        params.limits = PlLim.correctlimit(fig_desc,[min_x max_x min_y max_y]);
        my_axis(params,fig_desc,cutoffs,l_ind);
        movetexts(buf,flags);
        if flags.develmode
            legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        end
        
        fig_name = [complextype '_' fig_desc '_' cur_method];
        set(gcf,'Name',fig_name,'FileName',fig_name);
        orient landscape
        export_fig([common_out_dir filesep fig_name],'-eps','-pdf','-transparent', '-r300');
%        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);        
        
if icycle==1         
        if iCPcutoff(1) 
            [X0, i1] = findlinroot( buf.x(iCPcutoff(1),:), buf.y(iCPcutoff(1),:) ); %deltarho %buf.x is already sorted
            keypoints.irc(2) = X0;
disp(['point2: irc=' num2str(X0,'%.3f') ', nearest irc=' num2str(buf.x(iCPcutoff(1),i1),'%.3f') ', deltarho=' num2str(buf.y(iCPcutoff(1),i1),'%.3f') ]);
            keypoints.inddb(2)=sort_ind(i1);
            keypoints.desc{2}=['\Delta\rho_{' buf.legs{iCPcutoff(1)} '} = 0'];
%            disp(['P2 deltarho(' buf.legs{iCPcutoff(1)} ')=0: irc=' num2str(keypoints.irc(2),3) ])
        end
        if iCPcutoff(2)
            [X0, i1] = findlinroot( buf.x(iCPcutoff(1),:), buf.y(iCPcutoff(2),:) ); %deltarho %buf.x is already sorted
            keypoints.irc(4) = X0;
disp(['point4: irc=' num2str(X0,'%.3f') ', nearest irc=' num2str(buf.x(iCPcutoff(2),i1),'%.3f') ', deltarho=' num2str(buf.y(iCPcutoff(2),i1),'%.3f') ]);
            keypoints.inddb(4)=sort_ind(i1);
            keypoints.desc{4}=['\Delta\rho_{' buf.legs{iCPcutoff(2)} '} = 0'];
%            disp(['P4 deltarho(' buf.legs{iCPcutoff(2)} ')=0: irc=' num2str(keypoints.irc(4),3) ])
        end

        if iCPcutoff(3) 
            [X0, i1] = findlinroot( buf.x(iCPcutoff(1),:), buf.y(iCPcutoff(3),:) ); %deltarho %buf.x is already sorted
            keypoints.irc(6) = X0;
disp(['point6: irc=' num2str(X0,'%.3f') ', nearest irc=' num2str(buf.x(iCPcutoff(3),i1),'%.3f') ', deltarho=' num2str(buf.y(iCPcutoff(3),i1),'%.3f') ]);
            keypoints.inddb(6)=sort_ind(i1);
            keypoints.desc{6}=['\Delta\rho_{' buf.legs{iCPcutoff(3)} '} = 0'];
%            disp(['P2 deltarho(' buf.legs{iCPcutoff(1)} ')=0: irc=' num2str(keypoints.irc(2),3) ])
        end
        if iCPcutoff(4)
            [X0, i1] = findlinroot( buf.x(iCPcutoff(1),:), buf.y(iCPcutoff(4),:) ); %deltarho %buf.x is already sorted
            keypoints.irc(8) = X0;
disp(['point8: irc=' num2str(X0,'%.3f') ', nearest irc=' num2str(buf.x(iCPcutoff(4),i1),'%.3f') ', deltarho=' num2str(buf.y(iCPcutoff(4),i1),'%.3f') ]);
            keypoints.inddb(8)=sort_ind(i1);
            keypoints.desc{8}=['\Delta\rho_{' buf.legs{iCPcutoff(4)} '} = 0'];
%            disp(['P4 deltarho(' buf.legs{iCPcutoff(2)} ')=0: irc=' num2str(keypoints.irc(4),3) ])
        end
end %icycle
        %------------------------ ellipticity --------------------------  
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};buf.ht=[];
        cla%reset(gca)
        hold on
        for iCP=1:size(CPs,1) % over all H-bond CPs
            row = CPs(iCP,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                ircrec=ircdb(sort_ind(iirc));
                if ~isempty(ircrec.AIM)
                    iii = find( sum( repmat(row,size(ircrec.AIM.atoms,1),1)==ircrec.AIM.atoms, 2) == 2 ); %find appropriate CP index among molecule CPs
                    if iii
                        y(iirc) = ircrec.AIM.BondEl(iii);%#ok
                    end
                end
            end
%            leg=[atomlabel(ircdb(1),uniqCPind(iCP,1)) ',' atomlabel(ircdb(1),uniqCPind(iCP,2))];%#ok
            if icycle==2, dlm='...'; else dlm=''; end
            if ircdb(1).labels{row(1)}=='H'
                leg=[molpind{row(2)} dlm molpind{row(1)}];%#ok
            else
                leg=[molpind{row(1)} dlm molpind{row(2)}];%#ok
            end
            txt = ['\epsilon (' leg ')'];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            if iCP==1
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1),1);
            else
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1));
            end
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
            if flags.develmode || size(CPs,1)>1
                buf.ht(end+1) = text(x(10),min(max_y,y(10)+0.05*(max_y-min_y)),leg,'FontSize',label_fontsize);
            end
        end
        if icycle==1, fig_desc='epsilon'; else fig_desc='epsilonCHB'; end
        params={};
        params.flags=flags;
        params.limits = PlLim.correctlimit(fig_desc,[min_x max_x min_y max_y]);
        my_axis(params,fig_desc,cutoffs,l_ind);
        movetexts(buf,flags);
        if isfield(flags,'semilogy_for_ellipticity') && flags.semilogy_for_ellipticity
            set(a_tmp,'YScale','log');
        end
        if flags.develmode
            legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        end
        
        fig_name = [complextype '_' fig_desc '_' cur_method];
        set(gcf,'Name',fig_name,'FileName',fig_name);
        orient landscape
        export_fig([common_out_dir filesep fig_name],'-eps','-pdf','-transparent', '-r300');
%        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);        
        set(a_tmp,'YScale','linear');

        %------------------------ E_{HB} --------------------------  
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};buf.ht=[];
        cla%reset(gca)
        hold on
        for iCP=1:size(CPs,1) % over all H-bond CPs
            row = CPs(iCP,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                ircrec=ircdb(sort_ind(iirc));
                if ~isempty(ircrec.AIM)
                    iii = find( sum( repmat(row,size(ircrec.AIM.atoms,1),1)==ircrec.AIM.atoms, 2) == 2 ); %find appropriate CP index among molecule CPs
                    if iii
                        if ircrec.AIM.DelSqRho(iii) > Hbond_DelSqRho_crit_value % -0.020 is for plotting E_HB near critical points
                            y(iirc) = -0.5*CC.encoef*ircrec.AIM.V(iii);%#ok
                        else
                            y(iirc) = NaN;
                        end
                    end
                end
            end
%            leg=[atomlabel(ircdb(1),uniqCPind(iCP,1)) ',' atomlabel(ircdb(1),uniqCPind(iCP,2))];%#ok
            if icycle==2, dlm='...'; else dlm=''; end
            if ircdb(1).labels{row(1)}=='H'
                leg=[molpind{row(2)} dlm molpind{row(1)}];%#ok
            else
                leg=[molpind{row(1)} dlm molpind{row(2)}];%#ok
            end
            txt = ['E_{HB},k/m (' leg ')'];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            if iCP==1
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1),1);
            else
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1));
            end
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
            if flags.develmode || size(CPs,1)>1
                buf.ht(end+1) = text(x(10),min(max_y,y(10)+0.05*(max_y-min_y)),leg,'FontSize',label_fontsize);
            end
        end
%        max_y = min([max_y 100]); %cutting off covalent bonds energy range
        if icycle==1, fig_desc='EHB'; else fig_desc='EHBCHB'; end
        params={};
        params.flags=flags;
        params.limits = PlLim.correctlimit(fig_desc,[min_x max_x min_y max_y]);
        my_axis(params,fig_desc,cutoffs,l_ind);
        movetexts(buf,flags);
        if flags.develmode
            h_leg = legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        end
        
        fig_name = [complextype '_' fig_desc '_' cur_method];
        set(gcf,'Name',fig_name,'FileName',fig_name);
        orient landscape
        export_fig([common_out_dir filesep fig_name],'-eps','-pdf','-transparent', '-r300');
%        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);        

        %------------------------ E_{HB} Nikolayenko--------------------------  
if 0 && icycle==1 %not for CH..B bonds
        min_y=NaN; max_y=NaN;
        buf={};buf.h_plots=[];buf.legs={};buf.ht=[];
        cla%reset(gca)
        hold on
        for iCP=1:size(CPs,1) % over all H-bond CPs
            row = CPs(iCP,:);
            y=repmat(NaN,size(x));
            for iirc=1:numel(sort_ind) % over IRC points
                ircrec=ircdb(sort_ind(iirc));
                if ~isempty(ircrec.AIM)
                    iii = find( sum( repmat(row,size(ircrec.AIM.atoms,1),1)==ircrec.AIM.atoms, 2) == 2 ); %find appropriate CP index among molecule CPs
                    if iii
% E(OH···O) = -3.09 + 239·?(cp), 
% E(OH···N) = 1.72 + 142·?(cp), 
% E(NH···O) = -2.03 + 225·?(cp), 
% E(OH···C) = -0.29 + 288·?(cp)                        
%                        if ircrec.AIM.DelSqRho(iii)>-0.020 % -0.020 is for plotting E_HB near critical points
                            y(iirc) = -2.03 + 225*ircrec.AIM.ro(iii); %E(NH···O)
%                        else
%                            y(iirc) = NaN;
%                        end
                    end
                end
            end
%            leg=[atomlabel(ircdb(1),uniqCPind(iCP,1)) ',' atomlabel(ircdb(1),uniqCPind(iCP,2))];%#ok
            if icycle==2, dlm='...'; else dlm=''; end
            if ircdb(1).labels{row(1)}=='H'
                leg=[molpind{row(2)} dlm molpind{row(1)}];%#ok
            else
                leg=[molpind{row(1)} dlm molpind{row(2)}];%#ok
            end
            txt = ['E_{HB},k/m (' leg ')'];
            data(1,col)={txt};   data(2:numel(ircdb)+1,col) = num2cell(y); col=col+1;
            
            if iCP==1
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1),1);
            else
                buf.h_plots(end+1)=irc_plot( a_tmp, x(points2plot==1),y(points2plot==1));
            end
            buf.legs{end+1}=leg;%#ok
            min_y = min([min_y min(y)]);
            max_y = max([max_y max(y)]);
            if flags.develmode || size(CPs,1)>1
                buf.ht(end+1) = text(x(10),min(max_y,y(10)+0.05*(max_y-min_y)),leg,'FontSize',label_fontsize);
            end
        end
%        max_y = min([max_y 100]); %cutting off covalent bonds energy range
        if icycle==1, fig_desc='EHB_timn'; else fig_desc='illegal'; end
        params={};
        params.flags=flags;
        params.limits = PlLim.correctlimit(fig_desc,[min_x max_x min_y max_y]);
        my_axis(params,fig_desc,cutoffs,l_ind);
        movetexts(buf,flags);
        if flags.develmode
            h_leg = legend(buf.h_plots, buf.legs, 'FontSize',7, 'Location','Best');
        end
        
        fig_name = [complextype '_' fig_desc '_' cur_method];
        set(gcf,'Name',fig_name,'FileName',fig_name);
        orient landscape
        export_fig([common_out_dir filesep fig_name],'-eps','-pdf','-transparent', '-r300');
%        print(f_tmp,'-dpsc2', '-append', '-r300', psfile);       
end

        if exist('h_leg','var'), delete(h_leg), end;
        
end %for icycle


    %plot IRC keypoints extract properties and shapshots  
        [xP,xP_sort_ind]=sort(keypoints.irc);
        ixP=keypoints.inddb(xP_sort_ind);
        xPdesc=keypoints.desc(xP_sort_ind);
        xPatoms=keypoints.atoms(xP_sort_ind);
        table={};
        row=0;

        disp(['cutoff points: ' num2str(xP,'%0.2f')]);

        snap_azimut=0;
        snap_elevation=90;
        snap_xdir=-1;
        snap_ydir=-1;

%         f_mol = figure;
%         ikey=1;
%         ms0 = ircdb(ixP(ikey));
%         ms0 = createbondtable(ms0,1);% create bond table with H bonds
%         plotmol( ms0, 'r', 1, 0, 0, gca ); %plot shapshot
%         axis([-8 8 -4 4]);
%         text( snap_dir*3, snap_dir*4.5, ['point ' int2str(ikey) ', IRC = ' num2str(xP(ikey),'%0.2f') ], 'FontSize',8)
%         text( snap_dir*0, snap_dir*6, xPdesc{ikey}, 'FontSize',7)
%         view(snap_azimut,snap_elevation);
%         print(f_mol,'-depsc2','-tiff','-r300',[indir filesep fnameshort '_snap_p' int2str(ikey)]);
%         close(f_mol)

        %next code tested on 1920x1080 resolution
        scnsize = get(0,'ScreenSize');
        scnxsize = scnsize(3);
        scnysize = scnsize(4);
%         if numel(xP)==6
%             scnysize = 300;
%         elseif numel(xP)==9
%             scnysize = 560;
%         end
        pos  = [bdwidth,... 
                bdwidth,...
                scnxsize - 2*bdwidth,...
                scnysize - (topbdwidth + bdwidth)];
        set(f_tmp, 'Position',pos);
        cla%reset(gca)    
        
        for ikey=1:numel(xP)
            if isnan(xP(ikey)), continue, end

            %Change to plotting by AIM CPs
            ms0 = ircdb(ixP(ikey));
            ms0 = createbondtable(ms0,1);% create bond table with H bonds
            AIM = ms0.AIM;
            if isempty(AIM)
                disp(['No AIM info found: point IRC=' num2str(xP(ikey),3) ' skipped']);
                continue;
            end
            iHbonds = find(AIM.is_hbond);
            if isempty(iHbonds)
                row=row+1;
                table(row).irc = ms0.irc;
                table(row).rho = '-';
                table(row).deltarho = '-';
                table(row).eps = '-';
                table(row).E_HB = '-';
                table(row).Hbond_desc = '-';
                table(row).dAB = '-';
                table(row).dAH = '-';
                table(row).dHB = '-';
                table(row).aAHB = '-';
                table(row).desc = xPdesc(ikey);
            end
            for iH=1:numel(iHbonds)
                row=row+1;
                table(row).irc = ms0.irc;
                table(row).rho = AIM.ro(iHbonds(iH));
                table(row).deltarho = AIM.DelSqRho(iHbonds(iH));
                table(row).eps = 100*AIM.BondEl(iHbonds(iH));
                table(row).E_HB = -0.5*CC.encoef*AIM.V(iHbonds(iH));

                atoms = AIM.atoms(iHbonds(iH),:);
                HBatoms=[];
%here bug lives
                if any(ms0.labels{atoms(1)}~='H') && any(ms0.labels{atoms(2)}~='H') %none of atoms are H
                    disp(['Contact with atoms ' int2str(atoms(1)) ',' int2str(atoms(2)) ' skipped']);
                    continue;
                elseif any(ms0.labels{atoms(1)}~='H') %second atom is H
                    HBatoms = [ms0.btB(ms0.btA==atoms(2)) ms0.btA(ms0.btB==atoms(2)) atoms(2) atoms(1)];%#ok
                elseif any(ms0.labels{atoms(2)}~='H') %first atom is H
                    HBatoms = [ms0.btB(ms0.btA==atoms(1)) ms0.btA(ms0.btB==atoms(1)) atoms(1) atoms(2)];%#ok
                else %both atoms are H
                   disp(['dihydrogen bond with atoms ' int2str(atoms(1)) ',' int2str(atoms(2)) ' skipped']);
                   continue;
                end
                if ms0.labels{HBatoms(3)}=='C' %correct name for CH...B bonds 
                   dummy = HBatoms(1); HBatoms(1)=HBatoms(3); HBatoms(3)=dummy;
                end
                table(row).Hbond_desc = [molpind{HBatoms(1)} molpind{HBatoms(2)} '...' molpind{HBatoms(3)}];
                table(row).dAB = adist(ms0,HBatoms(1),HBatoms(3));
                table(row).dAH = adist(ms0,HBatoms(1),HBatoms(2));
                table(row).dHB = adist(ms0,HBatoms(2),HBatoms(3));
                table(row).aAHB = valang(ms0,HBatoms(1),HBatoms(2),HBatoms(3));
                table(row).desc = xPdesc(ikey);

            end

            if ~isempty(xPatoms{ikey}) %cheating - averaging distances
                d1 = adist(ms0,xPatoms{ikey}(1,1),xPatoms{ikey}(1,2));
                d2 = adist(ms0,xPatoms{ikey}(2,1),xPatoms{ikey}(2,2));
                avgd = (d1+d2)/2;
                
                if strcmp(complextype,'HypHyp') %hack
                    avgd2 = (adist(ms0,8,15)+adist(ms0,15,16))/2;
                end
            end
                
            text2plot={}; %distances to plot on shapshot
            for I=1:size(uniqCPind,1)
                atoms=uniqCPind(I,:);
                

                if ~isempty(xPatoms{ikey}) && (all(xPatoms{ikey}(1,:)==atoms) || all(xPatoms{ikey}(2,:)==atoms))%cheating
                    d = avgd;
                elseif strcmp(complextype,'HypHyp') && ~isempty(xPatoms{ikey}) && (all([8 15]==atoms) || all([15 16]==atoms))
                    d = avgd2;
                else
                    d=adist(ms0,atoms(1),atoms(2));
                end
                text2plot(end+1).text = {num2str(d,'%0.3f')};
                pl_x = sum(ms0.x(atoms))/2;
                if pl_x>0, pl_x = pl_x-0.2; else pl_x = pl_x+0.2; end
                text2plot(end).x = pl_x;
                pl_dy = 1.4*(mod(I,2)-0.5); % displacemant to texts do not overlap
%                 if strncmp(fnameshort,'TS_Hyp-Cyt',10)
%                     pl_dy = -pl_dy;
%                 end
                text2plot(end).y = sum(ms0.y(atoms))/2  - pl_dy;

                text2plot(end).z = sum(ms0.z(atoms))/2;
            end

            ms1 = ircdb(ixP(ikey));
%            ms1 = createbondtable(ms1);% create bond table w/o H bonds
            %using AIM data for plotting bond
            ms1.btA=AIM.atoms(~AIM.is_hbond,1)';
            ms1.btB=AIM.atoms(~AIM.is_hbond,2)';
            ms1.HBlist = AIM.atoms(iHbonds,:);
            ms1.text2plot = text2plot;
%hacking
%             if ikey==3 && strcmp(fnameshort,'TS_Hyp-Cyt')
%                 ms1.btA(end+1:end+2) = [8 21];
%                 ms1.btB(end+1:end+2) = [26 26];
%             end

            subplot(3,3,ikey)    
            plotmol( ms1, 'r', 1, 0, 0, gca ); %plot shapshot
            axis([-8 8 -4 4]);
            if sum(~isnan(xP))<=6
                axis_xsz = 0.162; %0.16
                axis_ysz = 0.16;
                k=1.20;
                set(gca,'Position',[(ikey-1)*axis_xsz-0.01 1-axis_ysz*k axis_xsz*k axis_ysz*k]);
            elseif sum(~isnan(xP))<=9
                axis_xsz = 1/5; %0.2
                axis_ysz = 0.2;
                k=1.24;
                if ikey<6
                    set(gca,'Position',[(ikey-1)*axis_xsz-0.02 1-axis_ysz*k axis_xsz*k axis_ysz*k]);
                else
                    set(gca,'Position',[(ikey-5.5)*axis_xsz-0.02 1-2*axis_ysz*k axis_xsz*k axis_ysz*k]);
                end
            end
            
            %rotation hacks
            if strcmp(fnameshort,'TS_Hyp-Cyt_eps4')
                if ms0.irc<-0.3
                    snap_azimut=0;
                    snap_elevation=-91;
                    snap_xdir=1;
                    snap_ydir=-1;
                else
                    snap_azimut=0;
                    snap_elevation=-89;
                    snap_xdir=-1;
                    snap_ydir=1;
                end
            elseif strcmp(complextype,'HypThy') || strcmp(complextype,'HypCyt')
                snap_azimut=0;
                snap_elevation=92;
                snap_xdir=1;
                snap_ydir=1;
            end
            view(snap_azimut,snap_elevation);

            fontsz=14;
            ht = text( snap_xdir*2.9, snap_ydir*3.8, ['point ' int2str(ikey) ', ' xPdesc{ikey} ], 'FontSize',fontsz,'FontWeight','bold'); %#ok<NASGU>
            ht2 = text( snap_xdir*1, snap_ydir*4.8, num2str(xP(ikey),'%0.2f'), 'FontSize',fontsz,'FontWeight','bold');
%            text( snap_xdir*0, snap_ydir*6, num2str(xP(ikey),'%0.2f'), 'FontSize',8,'FontAngle','italic','FontWeight','bold');

        end


        fig_desc='snaps';
        fig_name = [complextype '_' fig_desc '_' cur_method];
        set(gcf,'Name',fig_name,'FileName',fig_name);
        orient landscape
       % pack
        
       if exist('cutoffs','var') %saving snapshots only if we have something to save
%           print(f_tmp,'-dpsc2', '-append', '-r300', psfile);        
            print(f_tmp,'-dpng','-r600',[indir filesep fnameshort '_snaps']);
%           export_fig([indir filesep fnameshort '_snaps'],'-png','-transparent', '-r200');
%        system(sprintf('ps2pdf -dEPSCrop %s.eps %s.pdf',filo,filo));
       end
        
        tablecell=cell(row+1,11);
        tablecell(1,:)=[{'keypoint'} {'H-bond'} {'IRC'} {'\rho,a.u.'} {'{\Delta}{\rho},a.u.'} {'100*\epsilon'} ...
                        {'d_{AB},A'} {'d_{AH},A'} {'d_{HB},A'} {'{\angle}AHB,\degree'} {'E_{HB},kcal/mol'}];
        for i=1:row
            tablecell(i+1,1)=table(i).desc;
        end
        if ~isempty(table)
            tablecell(2:end,2)={table.Hbond_desc};
            tablecell(2:end,3)={table.irc};
            tablecell(2:end,4)={table.rho};
            tablecell(2:end,5)={table.deltarho};
            tablecell(2:end,6)={table.eps};
            tablecell(2:end,7)={table.dAB};    
            tablecell(2:end,8)={table.dAH};
            tablecell(2:end,9)={table.dHB};
            tablecell(2:end,10)={table.aAHB};    
            tablecell(2:end,11)={table.E_HB};
%            tablecell %#ok
        end
        
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
        save(limits_filename,'PlLim');
        
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

        xlswrite1spec(File,tablecell,'Hbonds_in_keypoints','A1',params);
        
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

    diary off
%break;

end %l_ind=1:numfiles

