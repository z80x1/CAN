%punch2aim.m extract WFN data (AIMPAC) from Gamess punch file
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-12-16
% Created        R O Zhurakivsky 2008-12-16

format compact

%----------------------------------------------------------------
infile='D:\zhr\_diplom\bin\pcg71c\work\h2o_wfn\punch'
%----------------------------------------------------------------

    ffname=infile;

    dlm=strfind(ffname,'\');

    fnamefull = ffname(dlm(end)+1:end);

    if ~exist(ffname,'file')
      warning(['File ' ffname ' is not found']);
      return
    end

    ffoutname=[ffname '.wfn'];

    disp(['Loaded file: ' ffname])

try
    fid=fopen(ffname,'r'); 	
    fido=fopen(ffoutname,'w'); 


    fl_aimpac=0;
    while 1
        
        tline = fgetl(fid);
        if ~ischar(tline)
    	    break
        end
        if isempty(tline)
            continue
        end

        if ~isempty(strfind(tline,'TOP OF INPUT FILE FOR BADER'))
          fl_aimpac=1;
          continue;
        end
        if ~isempty(strfind(tline,'END OF INPUT FILE FOR BADER'))
          fl_aimpac=0;
        end

        if fl_aimpac
    	   fprintf(fido,'%s\n',tline);
        end

    end

    fclose(fid);
    fclose(fido);
catch
    fclose(fid);
    fclose(fido);
end

disp(['Complete!'])

