%convert some types of Lowrance type MP file objects
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-10-01 
% Created        R O Zhurakivsky 2007-?-?

%2008-0825 added 4b, 4d polygone types

%2008-0131
%Polyline:
%0x3e       Directed text      - 2remove  ???
%Polygons:
%0x00  - 2remove ???


clear 
format compact
time0 = cputime;
tic

%---- enter path to your file here -----------------------------
workdir='D:\navigation\install\_maps\belarus\'
filename='00981003_lr.mp'
%---------------------------------------------------------------

ffname=fullfile(workdir,filename)   %#ok

dlm=strfind(ffname,'.');
ffoutname = [ffname(1:dlm(end)-1) '_out' ffname(dlm(end):end)]  %#ok

if ~exist(ffname,'file')==2
  error('file not found')
end

%Current convertion map for Garmin -> Lowrance objects types convertion
%Polygons:
%0x01 Large urban area (>200k)       -> 0x01 Urban area              OK
%0x02 Small urban area (<200k)       -> {0x01 Urban area}            ?
%0x03 Rural housing area             -> {0x01 Urban area}            ?
%0x04 Military base                  -> {0x01 Urban area}            ???
%0x05 Parking lot                    -> {0x01 Urban area}            ?
%0x06 Parking garage                 -> {0x01 Urban area}            ?
%0x07 Airport                        -> {0x01 Urban area}            ???
%0x08 Shopping center                -> {0x01 Urban area}            ???
%0x0a University/collage             -> {0x01 Urban area}            ???
%0x0b Hospital                       -> {0x01 Urban area}            ???
%0x0c Industrial complex             -> {0x01 Urban area}            ???
%0x0d Reservation                    -> {0x25 National park}         ???
%0x0e Airport runway                 -> {0x01 Urban area}            ???
%0x13 Building/Man-made area         -> {0x01 Urban area}            ???
%0x17 City park                      -> {0x25 National park}         ???
%0x19 Sport complex                  -> {0x01 Urban area}            ??? %added 2008-0131
%0x1a Cemetery                       -> {0x25 National park}         ???
%0x28 Sea/ocean                      -> {0x10 Lake/River}            ??
%0x32 Sea                            -> {0x10 Lake/River}            ??
%0x3d Large lake (77-250km2)         -> 0x10 Lake/River
%0x3f Medium lake (11-25km2)         -> 0x10 Lake/River
%0x40 Small lake (0.25-11km2)        -> 0x10 Lake/River
%0x41 Small lake (<0.25km2)          -> 0x10 Lake/River
%0x45 Blue-unknown                   -> 0x10 Lake/River
%0x46 Major river (>1km)             -> 0x10 Lake/River
%0x48 Medium river (40-200m)         -> 0x10 Lake/River
%0x4a Map selection area             -> 2delete???
%0x4b Backgroung                     - delete this object - no convertion needed!
%0x4c Intermittent water             -> {0x10 Lake/River}            ???
%0x4e Orchard/plantation             -> {0x25 National park}         ???
%0x4f Scrub                          -> {0x20 National forest}       ???
%0x50 Forest                         -> {0x20 National forest}       ???
%0x51 Wetland/swamp                  -> {0x10 Lake/River}            ???
%0x53 Sand/tidal/mud flat            -> {0x20 National forest}       ???


%Polyline:
%0x00 Road                           -> {0x42 Trail - Bike}          ?
%0x01 Major highway                  -> 0x01 Highway - Interstate
%0x02 Principal highway              -> 0x03 Highway - US
%0x03 Other highway road             -> 0x05 Highway - State
%0x04 Arterial road                  -> 0x10 Road rural
%0x05 Collector road                 -> 0x15 Street - Major city     OK
%0x06 Residential street             -> 0x17 Street - Minor city     OK
%0x07 Alleyway/private driveway      -> {0x17 Street - Minor city}   ?
%0x08 Highway ramp, low-speed        -> {0x17 Street - Minor city}   ?
%0x09 Highway ramp/high-speed        -> {0x17 Street - Minor city}   ?
%0x0a Unpaved road                   -> {0x10 Street - Minor city}   ?
%0x0b Major highway connector        -> {0x05 Street - Minor city}   ?
%0x0c Roundabout                     -> {0x05 Street - Minor city}   ?
%0x0e (tunnel path)                  -> {0x10 Street - Minor city}   ?
%0x14 Railroad                       -> 0x70 Railroad                OK!
%0x16 Walkway/trail                  -> {0x42 Trail - Bike}          ?
%0x18 Stream                         -> {0x32 Stream}                OK!
%0x1a Water or rail ferry            -> {0x05 Highway - State}       ???  %added 2008-0131
%0x1b Water or rail ferry            -> {0x05 Highway - State}       ???
%0x1c State/province boundary        -> 0x7a State/province boundary OK!
%0x1d Country/Parish boundary        -> {0x54 Elevation}             ???
%0x1e International boundary         -> {}
%0x1f River                          -> {0x32 Stream}                ?
%0x22 Major land contour (1/1)       -> 0x66 Elevation - 3           ???
%0x26 Intermittent stream/ditch      -> {0x30 Stream - Intermittent} OK!
%0x27 Airport runway centerline      -> 0x01 Highway - Interstate
%0x28 Oil or water pipeline          -> {0x50 Contour - 0}           ???
%0x29 Powerline                      -> {0x52 Contour - 1}           ???
%


%Points:
%0x0000 ??? all unconverted point became zero type   -> {0x0600}
%0x0200 Large city (5-10M)           
%0x0300 Large city (2-5M)            -> 0x0100 Large city (over 10M)
%0x0400 Large city (1-2M)            -> 0x0100 Large city (over 10M)
%0x0500 Medium city (0.5-1M)         -> 0x0200 Large city (5-10M)
%0x0600 Medium city (200-500k)       -> 0x0200 Large city (5-10M)
%0x0700 Medium city (100-200k)       -> 0x0200 Large city (5-10M)
%0x0800 Small city (50-100k)         -> 0x0300 Medium city (0.5-1M)
%0x0900 Small city (20-50k)          -> 0x0300 Medium city (0.5-1M)
%0x0a00 Small city/town (10-20k)     -> 0x0400 Medium city (100-200k)
%0x0b00 Small city/town (5-10k)      -> 0x0400 Medium city (100-200k)
%0x0c00 Settlement (2-5k)            -> 0x0500 Small city (50-100k)
%0x0d00 Settlement (1-2k)            -> 0x0500 Small city (50-100k)
%0x1100 Settlement (less 100)        -> 0x0600 Small city (20-50k)
%0x2a03 Dining (barbecue)
%0x2a0e Dining (cafe/diner)
%0x2b00 Lodging
%0x2b01 Hotel/motel
%0x2b02 Bed&breakfast inn
%0x2b03 Camping/RV park
%0x2b04 Resort
%0x2c00 Attraction
%0x2c02 Museum/History
%0x2c04 Landmark
%0x2c07 Zoo/Aquarium
%0x2d00 Entertainment
%0x2d06 skiing centre/resort
%0x2e00 Shopping
%0x2e04 Shopping center
%0X2E0A Speciality retail
%0x2f01 Autofuel
%0x2f03 Autorepair
%0x2f04 Air transportation
%0x2f08 Ground transportation
%0x2f09 Marina/boat repair
%0x2f0c Rest area / tourist info
%0x2f0d Automobile club
%0x2f11 Business center
%0x3001 Police station
%0x3002 Hospital                     -> 0x3801 Hospirtal             OK!
%0x3006 Border crossing
%0x4100 Fishing spot
%0x4400 Gas station
%0x5000 Drinking water
%0x5200 Scenic area
%0x5300 Skiing
%0x5500 Dam
%0x5c00 Diving area
%0x6401 Bridge
%0x6402 Building
%0x6406 Crossing
%0x6508 Falls
%0x6511 Spring
%0x6610 Plain
%0x6611 Range
%0x6613 Ridge
%0x6614 Rock
%0x6616 Summit


%Points
garmintypes{1}  =num2cell(hex2dec([...
 '00';...   %!!!all unknown types
 ]));
lowrancetypes{1}=num2cell(hex2dec([...
 '600';...  %Small city (20-50k)
 ]));

%Polylines
garmintypes{2}  =num2cell(hex2dec([...
 '00';...   %Road
 '07';...   %Alleyway/private driveway
 '08';...   %Highway ramp, low-speed
 '09';...   %Highway ramp/high-speed
 '0a';...   %Unpaved road
 '0b';...   %Major highway connector
 '0c';...   %Roundabout
 '0e';...   %(tunnel path)
 '16';...   %Walkway/trail
 '18';...   %Stream
 '1a';...   %Water or rail ferry
 '1b';...   %Water or rail ferry
 '1d';...   %Country/Parish boundary
 '1f';...   %River
 '26';...   %Intermittent stream/ditch
 '27';...   %Airport runway centerline
 '28';...   %Oil or water pipeline
 '29';...   %Powerline
 ]));
lowrancetypes{2}=num2cell(hex2dec([...
 '42';...   %Trail-bike
 '17';...   %Street - Minor city
 '17';...   %Street - Minor city
 '05';...   %Highway - State
 '10';...   %Road rural
 '05';...   %Highway - State
 '05';...   %Highway - State
 '10';...   %Road rural
 '42';...   %Trail-bike
 '32';...   %Stream
 '05';...   %Highway - State
 '05';...   %Highway - State
 '54';...   %Contour2
 '32';...   %Stream
 '30';...   %Stream - Intermittent
 '01';...   %Highway - Interstate
 '50';...   %Contour - 0
 '52';...   %Contour - 1
 ]));

%Polygones
garmintypes{3}  =num2cell(hex2dec([...
 '02';...   %Small urban area (<200k)
 '03';...   %Rural housing area 070306
 '04';...   %Military base
 '05';...   %Parking lot   
 '06';...   %Parking garage
 '07';...   %Airport
 '08';...   %Shopping center   
 '0a';...   %University/collage
 '0b';...   %Hospital
 '0c';...   %Industrial complex
 '0d';...   %Reservation
 '0e';...   %Airport runway
 '13';...   %Building/Man-made area
 '17';...   %City park
 '19';...   %Sport complex 
 '1a';...   %Cemetery
 '28';...   %Sea/ocean
 '32';...   %Sea
 '45';...   %Blue-unknown
 '4b';...   %Map coverege area
 '4c';...   %Intermittent water
 '4d';...   %Glacier
 '4e';...   %Orchard/plantation
 '4f';...   %Scrub
 '50';...   %Forest
 '51';...   %Wetland/swamp
 '52';...   %Tundra 070306
 '53';...   %Sand/tidal/mud flat
 ]));
lowrancetypes{3}=num2cell(hex2dec([...
 '01';...   %Urban area
 '01';...   %Urban area
 '01';...   %Urban area
 '01';...   %Urban area
 '01';...   %Urban area
 '01';...   %Urban area
 '01';...   %Urban area
 '01';...   %Urban area
 '01';...   %Urban area
 '01';...   %Urban area
 '25';...   %National park
 '01';...   %Urban area
 '01';...   %Urban area
 '25';...   %National park
 '01';...   %Urban area
 '25';...   %National park
 '10';...   %Lake/River
 '15';...   %Ocean
 '10';...   %Lake/River
 '4b';...   %??? !!delete_this_type_in_Lowrance_maps!!
 '10';...   %Lake/River
 '12';...   %Swamp ?
 '25';...   %National park
 '20';...   %National forest
 '20';...   %National forest
 '12';...   %Swamp
 '25';...   %National park
 '20';...   %National forest
 ]));

%try
    fid=fopen(ffname,'r');  

    fido=fopen(ffoutname,'w'); 

    objecttype=0;
    lineind = 0;
    while 1
        lineind = lineind+1;
        
        tline = fgetl(fid);
        if ~ischar(tline)
            break
        end
        if isempty(tline) || numel(tline)<2
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
    fclose(fido);

%catch
%    disp(['Error: ' lasterror.message '. Last line ' lineind]);
%    fclose(fid);
%    fclose(fido);
%end

cputime-time0  %#ok
toc
