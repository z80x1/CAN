%???
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-11-14
% Created        R O Zhurakivsky 2006-?-?

format compact
atomsind

%----------------------------------

indir=[CD.xyzdir filesep 'rg40']
initfile='rg40_Aaba_tchi_000.tpl'
%initval=0;
step=1; %degrees
stepnum=20; %total number of jobs


var2change=[{'tchi'}];
%var2change=[{'tchi'} {'tchi2'}];
%----------------------------------

fullinitfile=[indir filesep initfile];
buf={};
try
    fidt = fopen(fullinitfile,'r');
    if fidt==-1
       error(['Can''t open file ' fullinitfile]);
    end

    while 1
        tline = fgetl(fidt);
        if ~ischar(tline), break, end
        buf(end+1)={tline};
    end

%    disp(buf');
    fclose(fidt);

catch
    disp(['error: ' lasterror.message]);
    if fidt~=-1
        fclose(fidt);
    end
    return
end

iline2change=[];
ivar2change=[];
ival2change=[];
for lineind = 1:numel(buf)
  for ivar=1:numel(var2change)
       token=regexp(buf{lineind}, [var2change{ivar} '\s+(-?\d+\.\d*)'],'tokens');
      if numel(token)
         iline2change(end+1)=lineind;
         ivar2change(end+1)=ivar;
         ival2change(end+1)=sscanf(token{1}{1},'%f');
      end
  end
end


[pathstr,name0,ext,versn]=fileparts(fullinitfile);

anglechanges=step*(ceil(-stepnum/2+1):floor(stepnum/2));

for curchange=anglechanges

    name=[name0 '_' var2change{1} num2str(curchange,'%+0.3d')];
    outfile=fullfile(pathstr,[name '.gjf' versn]);
    
    buf2=buf;
    for ivar=1:numel(iline2change)

       curval=curchange+ival2change(ivar);
       buf2{iline2change(ivar)}=regexprep(buf2{iline2change(ivar)},'-?\d+\.\d*',num2str(curval,'%0.5f'));
    end

%    try
        fido = fopen(outfile,'w');
        for lind=1:numel(buf2)
            if strfind(buf2{lind},'%chk=')
                fprintf(fido,'%s%s.chk\n',buf2{lind},name);
            else
                fprintf(fido,'%s\n',buf2{lind});
            end
        end
        fclose(fido);
%    catch
%        error('savemol:exception occured for output file');
%        fclose(fido);
%    end
end


