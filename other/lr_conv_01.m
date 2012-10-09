%convert some types of Lowrance type MP file objects
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-02-27 
% Created        R O Zhurakivsky 2007-?-?

clear 
format compact
time0 = cputime;
tic



%workdir='D:\work\GPS\060711\';
%filename='Ukr_beta_080706_4_lr.mp';
%workdir='D:\work\GPS\Crimea\';
%filename='vc_map-en-0.2_relief_lr.mp';
workdir='D:\work\navigation\work\maps\Ukraine_070224\'
filename='Ukraine_10022007_full_lev0+garmin.mp'

ffname=fullfile(workdir,filename)   %#ok

dlm=strfind(ffname,'.');
ffoutname = [ffname(1:dlm(end)-1) '_out' ffname(dlm(end):end)]  %#ok

if ~exist(ffname,'file')==2
  error('file not found')
end

%Points
garmintypes{1}  =num2cell(hex2dec('0'));
%Polylines
garmintypes{2}  =num2cell(hex2dec([...
 '07';...   %
 '0a';...   %Unpaved road
 '0b';...   %
 '0c';...
 '14';...
 '16';...   %Walkway/trail
 '1b';...   %Water or rail ferry
 '1f';...   %River
 '26';...   %Intermittent stream/ditch
 '08';...   %Highway ramp, low-speed
 '18';...   %Stream
 '28';...   %Oil or water pipeline
 '29';...   %Powerline
 '1c';...    %State/province boundary
 '1d';...   %Country/Parish boundary
 '00'...    %Road
 ]));
%Polygones
garmintypes{3}  =num2cell(hex2dec([...
 '02';...
 '28';...
 '4e';...
 '4f';...   %Scrub
 '51';...
 '53';...   
 '04';...
 '07';...
 '0b';...
 '0c';...
 '0d';...
 '13';...   %Building/Man-made area
 '17';...   %City park
 '1a';...   %Cemetery
 '4c';...   %Intermittent water
 '32';...    %Sea
 '50'...    %Forest
 ]));

%Points
lowrancetypes{1}=num2cell(hex2dec('600'));
lowrancetypes{2}=num2cell(hex2dec([...
 '17';...
 '10';...
 '05';...
 '05';...
 '54';...
 '42';...
 '05';...
 '32';...
 '30';...
 '17';...
 '32';...
 '50';...
 '52';...
 '40';...    %Trail-hiking
 '54';...   %Contour2
 '42'...    %Trail-bike
 ]));
lowrancetypes{3}=num2cell(hex2dec([...
 '01';...
 '10';...
 '25';...
 '20';...   %Natioanal forest
 '10';...
 '20';...
 '01';...
 '01';...
 '01';...
 '01';...
 '25';...
 '01';...
 '25';...
 '25';...   %National park
 '10';...   %Lake/River
 '15';...   %Ocean
 '20'...    %National forest
 ]));

%try
    fid=fopen(ffname,'r'); 	
%    ifile = importdata(ffname);
%    ifile = textscan(fid,'%s',inf,'bufsize',16384);
%    ifile = fread(fid,inf,'*uchar');
%    linenums = size(ifile,1);


    fido=fopen(ffoutname,'w'); 
    
%    linenum=0;
    objecttype=0;
    while 1
%        linenum=linenum+1;
%	if linenum>linenums
%	    break
%	end
%	tline=ifile{linenum};

        tline = fgetl(fid);
        if ~ischar(tline)
    	    break
        end

        if objecttype  %we are within object description section 
    	    otype = hex2dec(strread(tline,'Type=0x%s'));
            if ~isempty(otype)
                ind = find(cell2mat(garmintypes{objecttype})==otype);
                if ind
                  tlinenew=sprintf('Type=0x%x',lowrancetypes{objecttype}{ind});
            	  tline = tlinenew;
                end
                
            else
                warning('lr_conv:TYPEnotfound','TYPE keyword not found');
            end

            fprintf(fido,'%s\n',tline);
    	    objecttype=0;
            continue
        end

        if ~isempty(strfind(tline,'RGN')) && isempty(strfind(tline,'END'))
            if ~isempty(strfind(tline,'RGN10')) %point
                objecttype=1;
            elseif ~isempty(strfind(tline,'RGN20')) %cityname
                objecttype=1;
            elseif ~isempty(strfind(tline,'RGN40')) %polyline
                objecttype=2;
            elseif ~isempty(strfind(tline,'RGN80')) %polygon
                objecttype=3;
            end
        end
    	fprintf(fido,'%s\n',tline);
        continue

    end
    fclose(fid);
%    clear ifile
    fclose(fido);

%catch
%    fclose(fid);
%    fclose(fido);
%end


cputime-time0  %#ok
toc
