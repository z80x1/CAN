function [ms0,desc] = checknewconf(ms0,workdb,desc)

%CHECK IF CONFORMATION IS NEW ONE
%ms0 - molecule structure to analize
%desc is list of identifiers of existing conformations
%workdb is list of molecular structures
%returns new list of identifiers
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-01-21 
% Created        R O Zhurakivsky 2005-09-?

%2006-0510 removed 3 degrees limit on P change between confs to be assumed as similar ones
%2006-0511 P limit restored and is set to 30 degrees
numfiles = numel(workdb);
[iprev,inext,sdesc_prev,sdesc_next] = checksimdesc(ms0.prop.sdesc,desc);

if iprev
    if iprev<=numfiles
        dP0 = abs(ms0.prop.Pdeg - workdb(iprev).prop.Pdeg);
        dP=mod(dP0,360);
    else
        warning(['DB doesn''t contain info about #' int2str(iprev) ' conformation']);
        dP=0;
    end
elseif inext
    if inext<=numfiles
        dP0 = abs(ms0.prop.Pdeg - workdb(inext).prop.Pdeg);
        dP=mod(dP0,360);
    else
        warning(['DB doesn''t contain info about #' int2str(inext) ' conformation']);
        dP=0;
    end
end

if isfield(ms0.prop,'hta0') && ms0.prop.hta0<150  %exclude beta-ribose conformations (150 degree is critical value)
      ms0.new='B';
elseif strcmpar(desc,ms0.prop.sdesc) 
  disp([ms0.prop.sdesc ' exists'])
  ms0.new='N';
elseif iprev && dP<5 %check if differences in P value with similar conformation in neighbour P-angle families are more than 30 degrees
  disp([ms0.prop.sdesc ' is similar to ' sdesc_prev ' (dP = ' num2str(dP0,2) 'deg)'])
  ms0.new='P';
elseif inext && dP<5
  disp([ms0.prop.sdesc ' is similar to ' sdesc_next ' (dP = ' num2str(dP0,2) 'deg)'])
  ms0.new='X';
else
  ms0.new='Y';
  if isempty(desc)
    desc=ms0.prop.sdesc;
  else
    desc(end+1,:)=ms0.prop.sdesc;
  end
  disp([ms0.new ', ' ms0.prop.sdesc])
end
