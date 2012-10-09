function savemolgs(odir,mol,format,order,gausstemplfile,route_add)

%saves molecule mol to output directory odir in Gaussian format
%format: 0 - create cartesian XYZ file
%        1 - create Z-matrix ZMT file
%order: ROW vector with atoms indexes that specifies atoms order in Z-matrix (it can
%be not complete than atoms that aren't specified are not sorted)
%gausstemplfile: name of template file
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-02-06
% Created        R O Zhurakivsky 2006-?-?

%2009-07  added route_add parameter

atomsind

if nargin<2
    error('output directory and/or molecular structure are not specified');
end
if nargin<3
    format=3; %XYZ
end
if nargin<4
    order=[];
%   [xxx,order]=sortrows(mol.pind); 
end

if isfield(mol,'pind')
  or=createbondchain(mol,order);
elseif format==3 %for unknown molecules
  order=reshape(order,1,numel(order));
  or=[order,setdiff(1:mol.atomnum,order)];
else
  error('savemol: Cannot save unknown molecule by ZMT');
end


if nargin<5
   gausstemplfile = [CD.templatesdir filesep 'gstempl_01.gjf'];
elseif exist(gausstemplfile,'file')~=2
   error('Template file not found');
end
if nargin<6
    route_add='';
end

%try

    fidt = fopen(gausstemplfile,'r');
    fname = strcat(odir,filesep,mol.desc,'.gjf');
    fid = fopen(fname,'w');

    dlm=strfind(mol.desc,'.');
    if isempty(dlm);
      workname=mol.desc;
    else
      workname=mol.desc(1:(dlm-1));
    end
    flgeom=0;
    flrouteend=0;
    fltitleend=0;
    addwfnfilename=0;
    addchkfilename=0;
    fl_additional_tpl=0;

    linenum=0;

    while 1
        linenum=linenum+1;

        tline = fgetl(fidt);
        fl_routesection=0;

        if ~ischar(tline)

            flgeom=1;

        else
            tlinetemp = lower(deblank(tline));
            if ~isempty(tlinetemp) && tlinetemp(1)=='%' 
                if ~isempty(strfind(tlinetemp,'chk')) 
                    addchkfilename=1;
                end
            end 
            if ~isempty(tlinetemp) && tlinetemp(1)=='#' 
                fl_routesection=1;
                
                if ~isempty(strfind(tlinetemp,'output')) && ~isempty(strfind(tlinetemp,'wfn'))
                    addwfnfilename=1;
                end
                if ~isempty(strfind(tlinetemp,'geom')) && ~isempty(strfind(tlinetemp,'connect')) %connectivity must be in additional section
                    fl_additional_tpl=1;
                end
                if ~isempty(strfind(tlinetemp,'freq')) && ~isempty(regexpi(tlinetemp,'ReadIsotopes'))
                    fl_additional_tpl=1;
                elseif ~isempty(regexpi(tlinetemp,'/gen'))
        		    fl_additional_tpl=1;
                elseif ~isempty(regexpi(tlinetemp,'(,|=|\()Read(,|\s|\))'))  %SCRF=(CPCM,Read)
        		    fl_additional_tpl=1;
                end
            end 
            if isempty(tlinetemp)
               if ~flrouteend
                   flrouteend=1;
               else
                   fltitleend=1;
               end
            end
        end
            
        fl_masses=1;
      	if ~(isfield(mol,'masses') && any(mol.masses))
           fl_masses=0;
        end;

        if flgeom
            if format==3 %create Gaussian input file from template

                for i=or
        	        if fl_masses && mol.masses(i)
                      fprintf(fid,'%s(Iso=%d)  %-12.8f  %-12.8f  %-12.8f\n',mol.labels{i},mol.masses(i),mol.x(i),mol.y(i),mol.z(i));
        	        else
                      fprintf(fid,'%s          %-12.8f  %-12.8f  %-12.8f\n',mol.labels{i},mol.x(i),mol.y(i),mol.z(i));
        	        end
                end

            elseif format==4  %create Gaussian input file from template with Z-matrix    

                mol = createzmt(mol,or);


                i = or(1);
      	      if fl_masses && mol.masses(i)
                    fprintf(fid,'%s(Iso=%d)\n',mol.labels{i},mol.masses(i));
      	      else
                    fprintf(fid,'%s\n',mol.labels{i});
      	      end

                i = or(2);
      	      if fl_masses && mol.masses(i)
                    fprintf(fid,'%s(Iso=%d)  %2d  %-12.8f\n',mol.labels{i},mol.masses(i),1,mol.R(i));
      	      else
                    fprintf(fid,'%s  %2d  %-12.8f\n',mol.labels{i},1,mol.R(i));
      	      end

                i = or(3);
      	      if fl_masses && mol.masses(i)
                    fprintf(fid,'%s(Iso=%d)  %2d  %-12.8f  %2d  %-12.8f\n',mol.labels{i},mol.masses(i),2,mol.R(i),1,mol.alfa(i));
    	          else
                    fprintf(fid,'%s  %2d  %-12.8f  %2d  %-12.8f\n',mol.labels{i},2,mol.R(i),1,mol.alfa(i));
                end

                for i=or(4:end)
        	        if fl_masses && mol.masses(i)
                      fprintf(fid,'%s(Iso=%d)  %2d  %-12.8f  %2d  %-12.8f  %2d  %-12.8f\n',...
                          mol.labels{i},mol.masses(i),mol.iR(i),mol.R(i),mol.ialfa(i),mol.alfa(i),...
                          mol.ibeta(i),mol.beta(i));
    	            else
                      fprintf(fid,'%s  %2d  %-12.8f  %2d  %-12.8f  %2d  %-12.8f\n',...
                          mol.labels{i},mol.iR(i),mol.R(i),mol.ialfa(i),mol.alfa(i),...
                          mol.ibeta(i),mol.beta(i));
    	            end
                end

            elseif format==5  %create Gaussian input file from template with Z-matrix and variables list   

                mol = createzmt(mol,or);

                fl_masses=1;
        	      if ~(isfield(mol,'masses') && any(mol.masses))
                    fl_masses=0;
                end;

                ivar=1; iR=1; ialfa=1; ibeta=1; 

                i = or(1);
      	      if fl_masses && mol.masses(i)
                    fprintf(fid,'%s(Iso=%d)\n',mol.labels{i},mol.masses(i));
      	      else
                    fprintf(fid,'%s\n',mol.labels{i});
      	      end

                i = or(2);
                labR=['R' int2str(iR)]; varlist(ivar,1:3)={1 labR mol.R(i)}; ivar=ivar+1; iR=iR+1; 
      	      if fl_masses && mol.masses(i)
                    fprintf(fid,'%s(Iso=%d)\t%2d]\t%s\n',mol.labels{i},mol.masses(i),1,labR);
      	      else
                    fprintf(fid,'%s\t\t%2d\t%s\n',mol.labels{i},1,labR);
      	      end

                i = or(3);
                labR=['R' int2str(iR)];    varlist(ivar,1:3)={1 labR mol.R(i)};    ivar=ivar+1; iR=iR+1; 
                labA=['A' int2str(ialfa)]; varlist(ivar,1:3)={2 labA mol.alfa(i)}; ivar=ivar+1; ialfa=ialfa+1; 
      	      if fl_masses && mol.masses(i)
                    fprintf(fid,'%s(Iso=%d)\t%2d\t%s\t%2d\t%s\n',mol.labels{i},mol.masses(i),2,labR,1,labA);
    	          else
                    fprintf(fid,'%s\t\t%2d\t%s\t%2d\t%s\n',mol.labels{i},2,labR,1,labA);
                end

                for i=or(4:end)
                  labR=['R' int2str(iR)];    varlist(ivar,1:3)={1 labR mol.R(i)};    ivar=ivar+1; iR=iR+1; 
                  labA=['A' int2str(ialfa)]; varlist(ivar,1:3)={2 labA mol.alfa(i)}; ivar=ivar+1; ialfa=ialfa+1; 
                  labD=['D' int2str(ibeta)]; varlist(ivar,1:3)={3 labD mol.beta(i)}; ivar=ivar+1; ibeta=ibeta+1; 
        	        if fl_masses && mol.masses(i)
                      fprintf(fid,'%s(Iso=%d)\t%2d\t%s\t%2d\t%s\t%2d\t%s\n',...
                          mol.labels{i},mol.masses(i),mol.iR(i),labR,mol.ialfa(i),labA,mol.ibeta(i),labD);
    	            else
                      fprintf(fid,'%s\t\t%2d\t%s\t%2d\t%s\t%2d\t%s\n',...
                          mol.labels{i},mol.iR(i),labR,mol.ialfa(i),labA,mol.ibeta(i),labD);
    	            end
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

            else
                error('undefined format');
            end


            fprintf(fid,'\n');

            if fl_additional_tpl

    	       dlm=strfind(gausstemplfile,'.');
    	       if isempty(dlm);
    	           gausstemplfile_add=[gausstemplfile '_add'];
    	       else
    	           gausstemplfile_add=[gausstemplfile(1:(dlm(end)-1)) '_add' gausstemplfile(dlm(end):end)];
    	       end
    	       fidt2 = fopen(gausstemplfile_add,'r');
               if fidt2 < 0 
                    error(['file with additional parameters not found: ' gausstemplfile_add]);
               end

    	       while 1

    	           tline2 = fgetl(fidt2);
                   if ~ischar(tline2), break, end
    	           fprintf(fid,'%s\n',tline2);
    	       end

    	       fclose(fidt2);
            end

            if addwfnfilename
                tline=[workname '.wfn'];
                fprintf(fid,'%s\n',tline);
            end

            break
        end %if flgeom


        if fltitleend
            fltitleend=0;
            fprintf(fid,'%s\n',mol.desc);
        end

        if addchkfilename
            tline=[tline workname '.chk'];
    	    addchkfilename=0;
        end
%%% my
        if fl_routesection
            fprintf(fid,'%s %s\n',tline,route_add);
        else
            fprintf(fid,'%s\n',tline);
        end

    end

    fclose(fid);
    fclose(fidt);
%catch
%    error('savemol:exception occured');
%    fclose(fid);
%    fclose(fidt);
%end

