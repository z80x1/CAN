function savemolnw(odir,mol,format,order,nwchemtemplfile)

%saves molecule mol to output directory odir
%format: 0 - create cartesian XYZ file
%        1 - create Z-matrix ZMT file
%order: ROW vector with atoms indexes that specifies atoms order in Z-matrix (it can
%be not complete than atoms that aren't specified are not sorted)
%nwchemtemplfile: name of template file
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-03-17
% Created        R O Zhurakivsky 2006-?-?


if nargin<3
    format=3; %XYZ
end
if nargin<4
    order=[];
%   [xxx,order]=sortrows(mol.pind); 
end
%order=reshape(order,1,numel(order));
%or=[order,setdiff(1:mol.atomnum,order)];
or=createbondchain(mol,order);

if nargin<5
   nwchemtemplfile = 'gotempl_02.nw';
end


fidt = fopen(nwchemtemplfile,'r');
fname = strcat(odir,filesep,mol.desc,'.nw');
fid = fopen(fname,'w');

dlm=strfind(mol.desc,'.');
workname=mol.desc(1:(dlm-1));
flgeom=0;
driver_found=0;
while 1

    tline = fgetl(fidt);
    if ~ischar(tline), break, end

    if flgeom
      if format==3 %create NWChem input file from template

        for i=or
            fprintf(fid,'%c  %-12.8f  %-12.8f  %-12.8f\n',mol.labels(i),mol.x(i),mol.y(i),mol.z(i));
        end

      elseif format==4  %create NWChem input file from template with Z-matrix    

        fprintf(fid,'  zmatrix\n');
        if isfield(mol,'R')==0
          mol = createzmt(mol,or);
        end
        fprintf(fid,'    %c\n',mol.labels(or(1)));
        fprintf(fid,'    %c  %2d  %-12.8f\n',mol.labels(or(2)),1,mol.R(or(2)));
        i = or(3);
        fprintf(fid,'    %c  %2d  %-12.8f  %2d  %-12.8f\n',mol.labels(i),2,mol.R(i),1,mol.alfa(i));
        for i=or(4:end)
    
          fprintf(fid,'    %c  %2d  %-12.8f  %2d  %-12.8f  %2d  %-12.8f\n',...
                  mol.labels(i),mol.iR(i),mol.R(i),mol.ialfa(i),mol.alfa(i),...
                  mol.ibeta(i),mol.beta(i));
        end
        fprintf(fid,'  end\n');

      else
        error('undefined format');
      end
      flgeom=0;
    end

    if strncmp(tline,'start',5)
      [buf,workprefix]=strread(tline,'%s %s');
      tline = strcat(tline,mol.desc);
    elseif strncmp(tline,'geometry',8)
      flgeom=1;
    end
    
    if ~isempty(strfind(tline,'driver')), driver_found=1; end
    if ~isempty(strfind(tline,'end')), driver_found=0; end
    if ~isempty(strfind(tline,'xyz')) && driver_found
      tline=[tline ' ' workname];
    end

    fprintf(fid,'%s\n',tline);
end

fclose(fid);
fclose(fidt);
