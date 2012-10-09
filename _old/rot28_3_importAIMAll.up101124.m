%rot28_3_importAIMAll: import AIMAll .mgp file data to database
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2009-12-08
% Created        R O Zhurakivsky 2009-04-22

%091208 coordinates of CPs are now imported in angstroms

%V = Virial Field = Potential Energy Density = Trace of Stress Tensor
%G = Lagrangian Form of Kinetic Energy Density
%K = Hamiltonian Form of Kinetic Energy Density
%L = K - G = Lagrangian Density = (-1/4)DelSqRho


tic
clear 
format compact

global pind
atomsind
pindsdef

%---------------------------------
moltype=16 %#ok
usedpackage='Gaussian' %#ok
theory='dftV3' %#ok
onlyoriginal=1;  % process db with only original conformations
dbsuffix='' %#ok
%dbsuffix='chilimited' %#ok
flwritefile=1 %#ok %save db file

filestype = 'mgp' %#ok %may be mgp or crt
molconf=''; %fill when analyzing specific conformation dynamics
%molconf='EabcA';
fl_useparam4ident=0; %use some parameters (e.g. chi torsion) to differ conformers

indir='D:\_diplom\CAN2\data\r16\rg_dftV2_aimall' %#ok
%---------------------------------

diaryfname0=[indir filesep 'logfile'];
diaryfname=diaryfname0;
for i=2:inf
  if ~(exist(diaryfname,'file')==2)
    break
  end
  diaryfname = [diaryfname0 int2str(i)];
end
diary(diaryfname)

workdbname0=['r' int2str(moltype)] %#ok
if ~isempty(molconf)
    workdbname0=[workdbname0 '_' molconf];
end
workdbname=workdbname0;
if strcmp(usedpackage,'Gaussian')
  workdbname=[workdbname '_g'];
end
if ~strcmp(theory,'dft')
  workdbname=[workdbname '_' theory];
end
if onlyoriginal
    templ='_or';
    workdbname = [workdbname templ];
end
if ~isempty(dbsuffix)
    workdbname = [workdbname '_' dbsuffix];
end
workdbname=[CD.dbdir filesep workdbname '.mat'] %#ok


load(workdbname,'workdb')
recnum=numel(workdb);
sdesc={};
new=[];
processed=zeros(1,recnum);
for i=1:recnum
   sdesc{i} = workdb(i).prop.sdesc;
end
if isfield(workdb,'new')
	for i=1:recnum
	   new(i) = workdb(i).new;
    end
end

sdesc(new~='Y')={''}; %zhr100405

sfiles = dir(strcat(indir,filesep,['*.' filestype]));
numfiles = size(sfiles,1);
if ~numfiles, error('No files found!'), end

HBtotalind=0;
for f_ind=1:numfiles

    fname=sfiles(f_ind).name;
    dlm=strfind(fname,'.');
    fnameshort = fname(1:(dlm-1));

%try
    fid=fopen([indir filesep fname],'r'); % don't use text mode - this is dangerous for Unix files ;)
    if fid==-1
     error(['Can''t open file ' fname])
    end

    if strcmp(filestype,'mgp')

        ms0=struct();
        maxind=0;
        while 1
          tline = fgetl(fid);
          if feof(fid), break, end
          if ~isempty(strfind(tline,'Nuclear Charges and Cartesian Coordinates')), 
              for i=1:3, tline=fgets(fid); end  %skip 3 lines to reach "-----------------" line
              break, 
          end
        end
    %    tline=fgets(fid);
        A=[];
        line_ind=0;
        while 1
          tline = fgetl(fid);
          if feof(fid), break, end
          if isempty(tline)
              break;
          end
          line_ind=line_ind+1;

          %Br28      35.0          -4.8600825600E+00  -3.4574149200E+00     4.6582665000E-01
          pat = '^([A-za-z]+)\d+\s+(\d*\.*\d*)\s+(-?\d*\.*\d*E?[+-]?\d*)\s+(-?\d*\.*\d*E?[+-]?\d*)\s+(-?\d*\.*\d*E?[+-]?\d*)';
          A = regexp(tline, pat, 'tokens','once');
          ms0.labels(line_ind)=A(1);
          ms0.x(line_ind)=sscanf(A{3},'%f');
          ms0.y(line_ind)=sscanf(A{4},'%f');
          ms0.z(line_ind)=sscanf(A{5},'%f');
        end


    %    A=fscanf(fid,'%s%d%f%f%f%f',[6,inf]);
    %    if A(2,end)==0
    %        A=A(:,1:end-1);
    %    end
    %    labels=cell(size(A,2),1);
    %    for i=1:size(A,2)
    %        labels(i)={char(A(1,i))};
    %    end
    %    ms0.labels=labels;

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
        [ms0,status] = identmol(ms0,moltype);
        if status
          lastwarn
          continue
        end
        if moltype == 910 
            ms0=calcproperties910(ms0,moltype);
        else
            ms0=calcproperties(ms0,moltype);
        end

        if fl_useparam4ident
            ms0_ind=find(abs([workdb.param]-ms0.prop.tchi)<0.01); %!!!need to add comparison based on conformers geometry
        else
            ms0_ind=strcmpcellar(sdesc,ms0.prop.sdesc);
        end
        
        if isempty(ms0_ind)
            fclose(fid);
            warning(['sdesc ' ms0.prop.sdesc ' from file ' fnameshort ' is not found!']);
            continue;
        end
        disp({[int2str(ms0_ind) ':' ms0.prop.sdesc ]} );
        
        clear('CPs');
        CPind=0;

        while 1
          tline = fgetl(fid);
          if feof(fid), break, end
          if ~isempty(strfind(tline,'CP#')), break, end
        end
        while 1
            if feof(fid), break, end
            buf={};
        
            CPind=CPind+1;
            while 1
              buf{end+1}=tline;

              tline = fgetl(fid);
              if feof(fid), break, end
              if ~isempty(strfind(tline,'CP#'))
                 break
              end
            end

    %        A=sscanf(tline,'%s%d%s%s%f%f%f');%Coords
            %CP# 35     Coords = -3.65216489424352E+00 -2.33311379771431E-01 -1.68025020406790E+00
            pat = '\s+Coords\s+=\s+(-?\d*\.*\d+E?[+-]?\d*)\s+(-?\d*\.*\d+E?[+-]?\d*)\s+(-?\d*\.*\d+E?[+-]?\d*)';
            A = regexp(buf, pat, 'tokens','once');
            A=[A{:}];
            CP = struct([]); %zhr101102
            CP(1).ind=CPind;
            %In WFN file coordinates are in Bohrs %zhr091208
            CP.x=sscanf(A{1},'%f')*CC.l*1e10;
            CP.y=sscanf(A{2},'%f')*CC.l*1e10;
            CP.z=sscanf(A{3},'%f')*CC.l*1e10;

            
            %Type = (3,-1) BCP O6 C1
            pat = '\s+Type\s+=\s+\((\d?,[+-]?\d?)\)\s+(\w+)\s+([\w\s]+)';
            A = regexp(buf, pat, 'tokens','once');
            A=[A{:}];
            CP.type=A{1};

            %O6 C1
            pat = '[A-Za-z]+([0-9]+)';
            A = regexp(A{3}, pat, 'tokens');
            A=[A{:}];
            CP.atoms=A;
                    
            %Rho =  2.5209925018E-01
            pat = '\s+Rho\s+=\s+(-?\d*\.*\d+E?[+-]?\d*)';
            A = regexp(buf, pat, 'tokens','once');
            A=[A{:}];
            CP.rho=sscanf(A{1},'%f');
            
    %         tline = fgetl(fid);%GradRho
    %         tline = fgetl(fid);%HessRho_EigVals
    %         tline = fgetl(fid);%HessRho_EigVec1
    %         tline = fgetl(fid);%HessRho_EigVec2
    %         tline = fgetl(fid);%HessRho_EigVec3
    %         tline = fgetl(fid);%DelSqRho
            
            %DelSqRho = -4.7582247986E-01
            pat = '\s+DelSqRho\s+=\s+(-?\d*\.*\d+E?[+-]?\d*)';
            A = regexp(buf, pat, 'tokens','once');
            A=[A{:}];
            CP.DelSqRho=sscanf(A{1},'%f');
            
            %Bond Ellipticity =  2.4291149587E-02
            pat = '\s+Bond\ Ellipticity\s+=\s+(.*)';
%            -?\d*\.*\d+E?[+-]?\d*
            A = regexp(buf, pat, 'tokens','once');
            A=[A{:}];
            CP.BondEl=sscanf(A{1},'%f');
            
            %V = -5.8556505813E-01
            pat = '\s+V\s+=\s+(-?\d*\.*\d+E?[+-]?\d*)';
            A = regexp(buf, pat, 'tokens','once');
            A=[A{:}];
            CP.V=sscanf(A{1},'%f');
            %G
            pat = '\s+G\s+=\s+(-?\d*\.*\d+E?[+-]?\d*)';
            A = regexp(buf, pat, 'tokens','once');
            A=[A{:}];
            CP.G=sscanf(A{1},'%f');
            %K
            pat = '\s+K\s+=\s+(-?\d*\.*\d+E?[+-]?\d*)';
            A = regexp(buf, pat, 'tokens','once');
            A=[A{:}];
            CP.K=sscanf(A{1},'%f');
            %L
            pat = '\s+L\s+=\s+(-?\d*\.*\d+E?[+-]?\d*)';
            A = regexp(buf, pat, 'tokens','once');
            A=[A{:}];
            CP.L=sscanf(A{1},'%f');

    %        tline = fgetl(fid);%Vnuc
    %        tline = fgetl(fid);%DivSigma
             
    %        disp([{int2str(CP.ind)}  {num2str(CP.rho)} {num2str(CP.DelSqRho)}] );
            CPs(CPind)=CP;
        end

    elseif strcmp(filestype,'crt')
   
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
              warning('Atoms coords not found');
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
        [ms0,status] = identmol(ms0,moltype);
        if status
          lastwarn
          continue
        end
        if moltype == 910 
            ms0=calcproperties910(ms0,moltype);
        else
            ms0=calcproperties(ms0,moltype);
        end

        if fl_useparam4ident
            ms0_ind=find(abs([workdb.param]-ms0.prop.tchi)<0.01); %!!!need to add comparison based on conformers geometry
        else
            ms0_ind=strcmpcellar(sdesc,ms0.prop.sdesc);
        end
        
        if isempty(ms0_ind)
            fclose(fid);
            warning(['sdesc ' ms0.prop.sdesc ' from file ' fnameshort ' is not found!']);
            continue;
        end
        disp({[int2str(ms0_ind) ':' ms0.prop.sdesc ]} );
        
        clear('CPs');
        CPind=0;

        while 1
          tline = fgetl(fid);
          if feof(fid), break, end
          if ~isempty(strfind(tline,'CP#')), break, end
        end
        fl_searchended = 0;
        while 1
            if feof(fid), break, end
            buf={};
            CP=struct([]);
        
            CPind=CPind+1;
            fl_coordsfound = 0;
            while 1
              buf{end+1}=tline;

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
                    warning('CP coords not found');
                else
                    CP(1).ind=CPind;
                    %In WFN file coordinates are in Bohrs %zhr091208
                    CP.x=A(3,1)*CC.l*1e10;
                    CP.y=A(3,2)*CC.l*1e10;
                    CP.z=A(3,3)*CC.l*1e10;
                end

                fl_coordsfound = 0;
              end

              
              if ~isempty(strfind(tline,'CP#'))
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
                warning(['CP Rho for CP #' int2str(CPind) ' not found']);
            else
                CP.rho=sscanf(A{1},'%f');
            end
            
            %DelSq(Rho(r))           -6.1929738296E-01
            pat = '\s+DelSq\(Rho\(r\)\)\s+(-?\d*\.*\d+E?[+-]?\d*)';
            A = regexp(buf, pat, 'tokens','once');
            A=[A{:}];
            if isempty(A)
                warning(['CP DelSqRho for CP #' int2str(CPind) ' not found']);
            else
                CP.DelSqRho=sscanf(A{1},'%f');
            end
            
            %Bond Ellipticity =  2.4291149587E-02
%            pat = '\s+Bond\ Ellipticity\s+=\s+(-?\d*\.*\d+E?[+-]?\d*)';
%            A = regexp(buf, pat, 'tokens','once');
%            A=[A{:}];
%            CP.BondEl=sscanf(A{1},'%f');
            
            %V(r)                    -2.6699369723E-01
            pat = '\s+V\(r\)\s+(-?\d*\.*\d+E?[+-]?\d*)';
            A = regexp(buf, pat, 'tokens','once');
            A=[A{:}];
            if isempty(A)
                warning(['CP V(r) for CP #' int2str(CPind) ' not found']);
            else
                CP.V=sscanf(A{1},'%f');
            end
            %G
            pat = '\s+G\(r\)\s+(-?\d*\.*\d+E?[+-]?\d*)';
            A = regexp(buf, pat, 'tokens','once');
            A=[A{:}];
            if isempty(A)
                warning(['CP G(r) for CP #' int2str(CPind) ' not found']);
            else
                CP.G=sscanf(A{1},'%f');
            end
            %K
            pat = '\s+K\(r\)\s+(-?\d*\.*\d+E?[+-]?\d*)';
            A = regexp(buf, pat, 'tokens','once');
            A=[A{:}];
            if isempty(A)
                warning(['CP K(r) for CP #' int2str(CPind) ' not found']);
            else
                CP.K=sscanf(A{1},'%f');
            end
            %L
            pat = '\s+L\(r\)\s+(-?\d*\.*\d+E?[+-]?\d*)';
            A = regexp(buf, pat, 'tokens','once');
            A=[A{:}];
            if isempty(A)
                warning(['CP L(r) for CP #' int2str(CPind) ' not found']);
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
                    
            CPs(CPind)=CP;
        end

    else 
        error('Unknown filestype value');
    end
    fclose(fid);

%catch
%    fclose(fid);
%    error('Some error catched');
%end    

    AIM.desc={};
    AIM.ro=[];
    AIM.DelSqRho=[];
    AIM.pinds=[];
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
       aa = ms0.btB(find( ms0.btA==atoms(1))); %atoms connected to atoms(1) 
       if sum(aa==atoms(2))>0 %exclude covalent/ionic bonds
           continue
       end
 
       HBtotalind=HBtotalind+1;

       if ms0.labels{atoms(1)}~='H' && ms0.labels{atoms(2)}~='H'
           HBatoms = [atoms(1) atoms(2)];
           HBatomspinds = ms0.pind(HBatoms);
           HBatomslabels = pind.labels(HBatomspinds);
           bondstr=[HBatomslabels{1} '...' HBatomslabels{2}];
       elseif ms0.labels{atoms(1)}~='H'
           HBatoms = [ms0.btB(find(ms0.btA==atoms(2))) ms0.btA(find(ms0.btB==atoms(2))) atoms(2) atoms(1)];
           HBatomspinds = ms0.pind(HBatoms);
           HBatomslabels = pind.labels(HBatomspinds);
           bondstr=[HBatomslabels{1} HBatomslabels{2} '...' HBatomslabels{3}];
       elseif ms0.labels{atoms(2)}~='H'
           HBatoms = [ms0.btB(find(ms0.btA==atoms(1))) ms0.btA(find(ms0.btB==atoms(1))) atoms(1) atoms(2)];
           HBatomspinds = ms0.pind(HBatoms);
           HBatomslabels = pind.labels(HBatomspinds);
           bondstr=[HBatomslabels{1} HBatomslabels{2} '...' HBatomslabels{3}];
       else
           HBatoms = [ms0.btB(find(ms0.btA==atoms(1))) ms0.btA(find(ms0.btB==atoms(1))) atoms(1) atoms(2) ...
                      ms0.btB(find(ms0.btA==atoms(2))) ms0.btA(find(ms0.btB==atoms(2)))];
           HBatomspinds = ms0.pind(HBatoms);
           HBatomslabels = pind.labels(HBatomspinds);
           bondstr=[HBatomslabels{1} HBatomslabels{2} '...' HBatomslabels{3} HBatomslabels{4}];
       end

       AIM.desc(end+1)={bondstr};
       AIM.ro(end+1)=CPs(i).rho;
       AIM.DelSqRho(end+1)=CPs(i).DelSqRho;
%       AIM.BondEl(end+1)=CPs(i).BondEl;
       AIM.V(end+1)=CPs(i).V;
       AIM.G(end+1)=CPs(i).G;
       AIM.K(end+1)=CPs(i).K;
       AIM.L(end+1)=CPs(i).L;
       AIM.pinds(end+1,1:numel(HBatomspinds))=HBatomspinds;
       
       disp(bondstr)
    end

    %zhr100713
    if numel(ms0_ind)~=1
        warning([int2str(numel(ms0_ind)) ' structure with same parameter value ' num2str(ms0.prop.tchi,3) ' detected']);
        for ih=1:numel(ms0_ind)
            workdb(ms0_ind(ih)).AIM = AIM; %#ok
        end
    else
        workdb(ms0_ind).AIM = AIM; %#ok
    end
    
    processed(ms0_ind)=1;
%return

end

disp([workdbname ': Total ' int2str(HBtotalind) ' of H-bonds are imported.']);

if sum(processed==0)
    warning([int2str(sum(processed==0)) ' DB record(s) are not processsed: ' sdesc{find(processed==0)} ]);
end

if flwritefile
  dlm=strfind(workdbname,'.');
  workdbnameold=[workdbname(1:dlm(end)-1) '~' workdbname(dlm(end):end)];
  copyfile(workdbname,workdbnameold);

  save(workdbname,'workdb')
end


toc
diary off