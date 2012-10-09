%convert some types of Lowrance type MP file objects
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-05-03 
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
filename='Ukraine_10022007_full_lr.mp'

ffname=fullfile(workdir,filename)   %#ok

dlm=strfind(ffname,'.');
ffoutname = [ffname(1:dlm(end)-1) '_out' ffname(dlm(end):end)]  %#ok

if ~exist(ffname,'file')==2
  error('file not found')
end

%in lcmmapedit 0.62 following garmin types are not converted:

%points:    0!!
%lowrance                         garmin  lowrance
%0x0300 Large city (2-5M)       -> 0x1 -> 
%0x0400 Large city (1-2M)       -> 0x1 -> 0x0200
%0x0500 Medium city (0.5-1M)    -> 0x2 -> 
%0x0600 Medium city (200-500k)  -> 0x2 -> 
%0x0700 Medium city (100-200k)  -> 0x2 -> 0x0400
%0x0800 Small city (50-100k)    -> 0x3 -> 
%0x0900 Small city (20-50k)     -> 0x3 -> 0x0600
%0x0a00 Small city/town (10-20k)-> 0x4 -> 
%0x0b00 Small city/town (5-10k) -> 0x4 -> 0x0600
%0x0c00 Settlement (2-5k)       -> 0x5 -> 
%0x0d00 Settlement (1-2k)       -> 0x5 -> 0x0600

%
%polylines: 00,07,09,0a!!,0b,0c,16,18!!,1b,1f!!,26
%            ^  ^    ^    ^  ^  ^  ^    ^  ^    ^
%polygones: 02,1a,28,4e,4f,50,51,53
%            ^  ^  ^  ^ ^   ^  ^ ^

%Points
garmintypes{1}  =num2cell(hex2dec([...
 '00';...
 '01';...
 '02';...
 '03';...
 '04';...
 '05';...
 ]));
lowrancetypes{1}=num2cell(hex2dec([...
 '600';...
 '200';...
 '300';...
 '400';...
 '500';...
 '600';...
 ]));

%Polylines
garmintypes{2}  =num2cell(hex2dec([...
 '00';...    %01 Road
 '07';...   %02
 '08';...   %11 Highway ramp, low-speed
 '09';...   %17 Highway ramp/high-speed
 '0a';...   %03 Unpaved road
 '0b';...   %04 Major highway connector
 '0c';...   %05 Roundabout
 '16';...   %07 Walkway/trail
 '18';...   %12 Stream
 '1b';...   %08 Water or rail ferry
 '1d';...   %16 Country/Parish boundary
 '1f';...   %09 River
 '26';...   %10 Intermittent stream/ditch
 '28';...   %13 Oil or water pipeline
 '29';...   %14 Powerline
 ]));
% '14';...   %06 Railroad done!
% '1c';...   %15 State/province boundary
lowrancetypes{2}=num2cell(hex2dec([...
 '42';...    %   Trail-bike
 '17';...   %
 '17';...   %
 '05';...
 '10';...   %
 '05';...   %
 '05';...   %
 '42';...   %
 '32';...   %
 '05';...   %
 '54';...   %   Contour2
 '32';...   %
 '30';...   %
 '50';...   %
 '52';...   %
 ]));
% '54';...   %
% '40';...   %   Trail-hiking

%Polygones
% '4b';... %Background 070306
garmintypes{3}  =num2cell(hex2dec([...
 '02';...
 '03';... %Rural housing area 070306
 '04';...
 '07';...
 '0b';...
 '0c';...
 '0d';...
 '13';...   %Building/Man-made area
 '17';...   %City park
 '1a';...   %Cemetery
 '28';...
 '32';...    %Sea
 '4c';...   %Intermittent water
 '4e';...
 '4f';...   %Scrub
 '50';...    %Forest
 '51';...
 '52';... %Tundra 070306
 '53';...   
 ]));
lowrancetypes{3}=num2cell(hex2dec([...
 '01';...
 '01';... %!
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
 '10';...   %Lake/River
 '25';...
 '20';...   %National forest
 '20';...    %National forest
 '10';...
 '25';...   %national park
 '20';...
 ]));

try
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
        if isempty(tline)
            continue
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

        if tline(1:2)=='[P'
            if findstr(tline,'[POI]') %point
                objecttype=1;
            elseif findstr(tline,'[POLYLINE]') %polyline
                objecttype=2;
            elseif findstr(tline,'[POLYGON]') %polygon
                objecttype=3;
            end
        end
    	fprintf(fido,'%s\n',char(lr_translit(tline)));
        continue

    end
    fclose(fid);
%    clear ifile
    fclose(fido);

catch
    fclose(fid);
    fclose(fido);
end


cputime-time0  %#ok
toc
