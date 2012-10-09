%groupexcel.m Concatenate xls files in directory to one file
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2009-01-06 
% Created        R O Zhurakivsky 2008-12-?

format compact
realtytype2include=[];
roomnumber2include=[];
district2exclude=[];
floors_total_min=0;
streets2include=[];

%----------------------------------------------------------------
indir='D:\zhr\ID\Brokband\12 December'
operationtype = 'rent' 
fl_onlyfull   = 1  %operate on only _full files
%realtytype2include=[{'Квартира'}];
realtytype2include=[{'Квартира'} {'Дом'} {'Гостинка'} {'Комната'}];
roomnumber2include=[1 2 3];
district2exclude=[{'Барышевский'} {'Бориспольский'} {'Броварской'} {'Васильковский'} {'Володарский'} {'Вышгородский'} {'Обуховcкий'} {'Обуховский'} {'Белоцерковский'} {'Бородянский'} {'Згуровский'} {'Переяслав-Хмельницкий'} {'Фастовский'} {'Мироновский'} {'Богуславский'} {'Яготинский'}];
%floors_total_min=0;
%streets2include=[{'Чистяк'}];
%----------------------------------------------------------------

sfiles=[];
sfiles=getrecfilelist(indir,'xls');  %sfiles now contains full filenames

numfiles = numel(sfiles);
if ~numfiles
  error('No input files found');
end


buf={};
data={};
enum=0;
etxt='';
fprocessed = 0;
for i=1:numfiles
    lineincluded=0;

    dlm=strfind(sfiles(i).name,'.');

    fnamefull = sfiles(i).name(1:(dlm(end)-1));

    ffname = fullfile(sfiles(i).path,sfiles(i).name); 

    if ~strncmpi(fnamefull,operationtype,numel(operationtype))
        continue;
    end
    if fl_onlyfull
      fullstr='_full';
      if ~strcmpi(fnamefull(end-numel(fullstr)+1:end),fullstr)
        continue;
      end
    end

    disp(['Loading file: ' ffname])

    if ~exist(ffname,'file')
      warning(['File ' ffname ' is not found']);
      continue;
    end

    [enum,etxt,buf]=xlsread(ffname,1);

    if fprocessed==0
       startrow=1;
       ffnamebase=ffname; %file on which output file will be constructed
    else
       startrow=2;
    end

    for j=startrow:size(buf,1)
      
      try
      if ~any(cell2mat(buf(j,:)))
          continue
      end
      catch
      end
       
      fl_include=[];
      if ~isempty(realtytype2include)
          fl_include(end+1)=0; 
          for tt=1:numel(realtytype2include)
              if strcmpi(buf(j,2),realtytype2include{tt})
                fl_include(end)=1; 
                break;
              end
          end    
      end
      if ~isempty(roomnumber2include) && any(roomnumber2include)
          fl_include(end+1)=0; 
          if ~all(isnan(buf{j,3})) && any(str2num(buf{j,3})) && any(str2num(buf{j,3})==roomnumber2include)
              fl_include(end)=1; 
          end
      end
      if ~isempty(district2exclude)
          district=lower(buf(j,6));
          fl_include(end+1)=1; 
          for tt=1:numel(district2exclude)
              d2exclude=lower(district2exclude{tt});
              if strncmp(district,d2exclude,numel(d2exclude))
                fl_include(end)=0; 
                break;
              end
          end    
      end
      if floors_total_min
          fl_include(end+1)=0; 
          if ~all(isnan(buf{j,13})) && any(str2num(buf{j,13})>=floors_total_min)
              fl_include(end)=1; 
          end
      end
      if ~isempty(streets2include)
          street=lower(buf(j,7));
          fl_include(end+1)=0; 
          for tt=1:numel(streets2include)
              street2include=lower(streets2include{tt});
              if strncmp(street,street2include,numel(street2include))
                fl_include(end)=1; 
                break;
              end
          end    
      end



      if all(fl_include)
          data(end+1,:)=buf(j,:);
          lineincluded=lineincluded+1;
      elseif ~fprocessed && j==startrow
          data(end+1,:)=buf(j,:);
      end
    end

    fprocessed=fprocessed+1;
    disp(['lines included: ' int2str(lineincluded) ', total lines: ' num2str(size(data,1)) ]);
end

if ~isempty(data)
      
      dlm=strfind(indir,'\');
      if dlm(end)==numel(indir) 
          i0=dlm(end-1)+1;
          i1=dlm(end)-1;
      else
          i0=dlm(end)+1;
          i1=numel(indir);
      end
      lastdirname=indir(i0:i1);
      
      ffnameout=[indir filesep '..' filesep operationtype '_' lastdirname '_complete' '.xls'];
      copyfile(ffnamebase,ffnameout);
      xlswrite(ffnameout,data,1,'A1');
end

disp(['Complete!'])
    %clear data buf etxt enum

