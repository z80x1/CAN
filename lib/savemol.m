function savemol(odir,mol,format,order,consttor,fl_Hlast)

%saves molecule mol to output directory odir
%format: 0 - create cartesian XYZ file
%        1 - create Z-matrix ZMT file
%		 2 - create ZMT file  (new version - constructed by bonds)
%		 3 - create ZMT file  (new version - constructed by bonds, with list of variables)
%		 4 - mixed format
%order: ROW vector with atoms indexes that specifies atoms order in Z-matrix (it can
%be not complete than atoms that aren't specified are not sorted)

%consttor: specifies that torsion angle between 4 atoms with indexes are specified
%is constant

%fl_Hlast=1 - will put H atoms at the end
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-08-07
% Created        R O Zhurakivsky 2006-?-?

if nargin<3
    format=0; %XYZ
end
if nargin<4
    order=[];
%   [xxx,order]=sortrows(mol.pind); 
end
if nargin<5
    consttor=[];
    prconst=0;
else
    prconst=1;
end
if nargin<6
   fl_Hlast=0;
end


if isfield(mol,'pind')
  or=createbondchain(mol,order,fl_Hlast);
elseif format==0 %for unknown molecules
  order=reshape(order,1,numel(order));
  or=[order,setdiff(1:mol.atomnum,order)];
else
  error('savemol: Cannot save unknown molecule by ZMT');
end


if format==0

    fname = strcat(odir,filesep,mol.desc,'.xyz');
    fid = fopen(fname,'w');
    if fid==-1
      error('fopen error');
    end

    for i=or
        fprintf(fid,'%s  %-12.8f  %-12.8f  %-12.8f\n',mol.labels{i},mol.x(i),mol.y(i),mol.z(i));
    end

    fclose(fid);
elseif format==1  %create ZMT file
warning('Are you sure you want to use format==1 ?');


    fname = strcat(odir,filesep,mol.desc,'.zmt');
    fid = fopen(fname,'w');


    j=0;
    for i=or
        j=j+1;
        if j==1
          fprintf(fid,'%c\n',mol.labels{i});
        elseif j==2
          R = adist(mol,i,or(1));
          fprintf(fid,'%c\t%d\t%-.8f\n',mol.labels{i},1,R);
        elseif j==3
          R = adist(mol,i,or(2));
          alpha = valang(mol,i,or(2),or(1));
          fprintf(fid,'%c\t%d\t%-.8f\t%d\t%-.8f\t\n',mol.labels{i},2,R,1,alpha);
        else

          R = adist(mol,i,or(j-1));
          alpha = valang(mol,i,or(j-1),or(j-2));
          if prconst==1 & (or(j-3:j)==consttor(1:4))
            T1=torang(mol,or(j),or(j-1),or(j-2),or(j-3));
            fprintf(fid,'%c\t%d\t%-.8f\t%d\t%-.8f\t%d\t%s\n',mol.labels{or(j)},j-1,R,j-2,alpha,j-3,'T1');
          else
            beta = torang(mol,i,or(j-1),or(j-2),or(j-3));
            fprintf(fid,'%c\t%d\t%-.8f\t%d\t%-.8f\t%d\t%-.8f\n',mol.labels{i},or(j-1),R,or(j-2),alpha,or(j-3),beta);
          end
        end
    end
    if prconst==1
       fprintf(fid,'constants\n');
       fprintf(fid,'T1\t%-.8f\n',T1);
    end
    fclose(fid);
elseif format==2  %create ZMT file  (new version - constructed by bonds)
    fname = strcat(odir,'/',mol.desc,'.zmt');
    fid = fopen(fname,'w');


    if isfield(mol,'R')==0
      mol = createzmt(mol,or);
    end
    fprintf(fid,'%c\n',mol.labels{or(1)});
    fprintf(fid,'%c  %2d  %-12.8f\n',mol.labels{or(2)},1,mol.R(or(2)));
    i = or(3);
    fprintf(fid,'%c  %2d  %-12.8f  %2d  %-12.8f\n',mol.labels{i},2,mol.R(i),1,mol.alfa(i));
    for i=or(4:end)

%          if prconst==1 & (or(j-3:j)==consttor(1:4))
%            T1=torang(mol,or(j),mol,or(j-1),mol,or(j-2),mol,or(j-3));
%            fprintf(fid,'%c\t%d\t%-.3f\t%d\t%-.3f\t%d\t%s\n',mol.labels(or(j)),j-1,R,j-2,alpha,j-3,'T1');
%          else
      fprintf(fid,'%c  %2d  %-12.8f  %2d  %-12.8f  %2d  %-12.8f\n',...
         mol.labels{i},mol.iR(i),mol.R(i),mol.ialfa(i),mol.alfa(i),...
         mol.ibeta(i),mol.beta(i));
%          end
    end
%    if prconst==1
%       fprintf(fid,'constants\n');
%       fprintf(fid,'T1\t%-.3f\n',T1);
%    end

    
    fclose(fid);

elseif format==3  %create ZMT file  (new version - constructed by bonds, with list of variables)
    fname = strcat(odir,filesep,mol.desc,'.zmt');
    fid = fopen(fname,'w');

    ivar=1; iR=1; ialfa=1; ibeta=1; 

    if isfield(mol,'R')==0
      mol = createzmt(mol,or);
    end
    fprintf(fid,'%c\n',mol.labels{or(1)});

    
    labR=['R' int2str(iR)]; varlist(ivar,1:3)={1 labR mol.R(or(2))}; ivar=ivar+1; iR=iR+1; 
    fprintf(fid,'%c\t%2d\t%s\n',mol.labels{or(2)},1,labR);

    i = or(3);
    labR=['R' int2str(iR)];    varlist(ivar,1:3)={1 labR mol.R(i)};    ivar=ivar+1; iR=iR+1; 
    labA=['A' int2str(ialfa)]; varlist(ivar,1:3)={2 labA mol.alfa(i)}; ivar=ivar+1; ialfa=ialfa+1; 
    fprintf(fid,'%c\t%2d\t%s\t%2d\t%s\n',mol.labels{i},2,labR,1,labA);

    for i=or(4:end)

      labR=['R' int2str(iR)];    varlist(ivar,1:3)={1 labR mol.R(i)};    ivar=ivar+1; iR=iR+1; 
      labA=['A' int2str(ialfa)]; varlist(ivar,1:3)={2 labA mol.alfa(i)}; ivar=ivar+1; ialfa=ialfa+1; 
      labD=['D' int2str(ibeta)]; varlist(ivar,1:3)={3 labD mol.beta(i)}; ivar=ivar+1; ibeta=ibeta+1; 
      fprintf(fid,'%c\t%2d\t%s\t%2d\t%s\t%2d\t%s\n',...
         mol.labels{i},mol.iR(i),labR,mol.ialfa(i),labA,mol.ibeta(i),labD);
    end

    fprintf(fid,'\tVariables:\n');
    for i=1:size(varlist,1);
      if varlist{i,1}==1
        fprintf(fid,'  %s\t\t%-12.8f\n',varlist{i,2},varlist{i,3});
      end
    end
    for i=1:size(varlist,1);
      if varlist{i,1}==2
        fprintf(fid,'  %s\t\t%-12.8f\n',varlist{i,2},varlist{i,3});
      end
    end
    for i=1:size(varlist,1);
      if varlist{i,1}==3
        fprintf(fid,'  %s\t\t%-12.8f\n',varlist{i,2},varlist{i,3});
      end
    end
    
    fclose(fid);

elseif format==4  %mixed format

    if isfield(mol,'coortype')
      indxyz=find(mol.coortype=='C'); %indexes of cartesian atoms
    else
      indxyz=[];
    end
    
    if numel(indxyz)==mol.atomnum
        fname = strcat(odir,filesep,mol.desc,'.xyz');
    else
        fname = strcat(odir,filesep,mol.desc,'.zmt');
    end
    fid = fopen(fname,'w');
    if fid==-1
      error('fopen error');
    end

    if ~isempty(order)
        order=[];
        warning('With format type = 4 order parameter will be zeroed!');
        or=createbondchain(mol,union(order,indxyz),fl_Hlast);
    end
    if isfield(mol,'R')==0
      mol = createzmt(mol,[]);
    end


    for i=indxyz
        fprintf(fid,'%c  %-12.8f  %-12.8f  %-12.8f\n',mol.labels{i},mol.x(i),mol.y(i),mol.z(i));
    end


    [or,I] = setdiff(or,indxyz);
    [XX,II] = sort(I);
    or = or(II);

    ivar=1; iR=1; ialfa=1; ibeta=1; 

if 0    
    i = or(1);
    fprintf(fid,'%c\n',mol.labels{i});
    i = or(2);
    labR=['R' int2str(iR)]; varlist(ivar,1:3)={1 labR mol.R(i)}; ivar=ivar+1; iR=iR+1; 
    fprintf(fid,'%c\t%2d\t%s\n',mol.labels{i},1,labR);
    i = or(3);
    labR=['R' int2str(iR)];    varlist(ivar,1:3)={1 labR mol.R(i)};    ivar=ivar+1; iR=iR+1; 
    labA=['A' int2str(ialfa)]; varlist(ivar,1:3)={2 labA mol.alfa(i)}; ivar=ivar+1; ialfa=ialfa+1; 
    fprintf(fid,'%c\t%2d\t%s\t%2d\t%s\n',mol.labels{i},2,labR,1,labA);
end
    
if numel(indxyz)<3
  error('Cannot build structure if cartezian atoms are less than 3');
end

%    for i=or(4:end)
    for i=or(1:end)

      labR=['R' int2str(iR)];    varlist(ivar,1:3)={1 labR mol.R(i)};    ivar=ivar+1; iR=iR+1; 
      labA=['A' int2str(ialfa)]; varlist(ivar,1:3)={2 labA mol.alfa(i)}; ivar=ivar+1; ialfa=ialfa+1; 
      labD=['D' int2str(ibeta)]; varlist(ivar,1:3)={3 labD mol.beta(i)}; ivar=ivar+1; ibeta=ibeta+1; 
      fprintf(fid,'%c\t%2d\t%s\t%2d\t%s\t%2d\t%s\n',...
         mol.labels{i},mol.iR(i),labR,mol.ialfa(i),labA,mol.ibeta(i),labD);
    end

    fprintf(fid,'\tVariables:\n');
    for i=1:size(varlist,1);
      if varlist{i,1}==1
        fprintf(fid,'  %s\t\t%-12.8f\n',varlist{i,2},varlist{i,3});
      end
    end
    for i=1:size(varlist,1);
      if varlist{i,1}==2
        fprintf(fid,'  %s\t\t%-12.8f\n',varlist{i,2},varlist{i,3});
      end
    end
    for i=1:size(varlist,1);
      if varlist{i,1}==3
        fprintf(fid,'  %s\t\t%-12.8f\n',varlist{i,2},varlist{i,3});
      end
    end
    
    fclose(fid);

else
  error('undefined format');
end

