function outfiles = getrecfilelist(indir,ext)
%returns recursive list if files with ext xtention inside directory indir
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-12-14 
% Created        R O Zhurakivsky 2008-12-14
    
    files = dir(strcat(indir));
    outfiles = [];
    for i=1:size(files,1)
      dlm=strfind(files(i).name,'.');
      
      if files(i).isdir 
        if ~strcmp(files(i).name,'.') && ~strcmp(files(i).name,'..')
          outfiles=[outfiles getrecfilelist([indir filesep files(i).name],ext)];
        end
      elseif ~isempty(dlm) && strcmp(files(i).name(dlm(end)+1:end),ext)
        if isempty(outfiles)
          outfiles=files(i);
        else
          outfiles(end+1).name=files(i).name;
          outfiles(end).date=files(i).date;
          outfiles(end).bytes=files(i).bytes;
          outfiles(end).isdir=files(i).isdir;
        end
        outfiles(end).path=indir;
      end

    end


end


