%make all possible conformations by rotating around ordinary bonds and
%create NWChem or Gaussian input file
%using conformation's database
%
% Version 1.01    
% Last modified  R O Zhurakivsky 2009-03-15
% Created        R O Zhurakivsky 2005-09-?

%2009-0315 while rotating round chi seaching of bN9 changed be earlier
%than searching of bN1
%2009-0421 savemolgs: format=4 changed to 5

clear
format compact
global pind
global flPlot
pindsdef
atomsind

%--------------------------------------
moltype=390         %#ok
workname='r39010'   %#ok
%theory='mp2V2'      %#ok
theory='dftV3'      %#ok
flPlot=1            %#ok

minHHdist=1.2;  %minimal H to H atoms distance to try such conformation
%--------------------------------------

gtemplname=[CD.templatesdir filesep workname '_templ.gjf'];
odir=[CD.xyzdir filesep workname];
if exist(odir,'dir')~=7
   mkdir(odir);
end

diaryfname0=[odir filesep 'logfile'];
diaryfname=diaryfname0;
for i=2:inf
  if ~(exist(diaryfname,'file')==2)
    break
  end
  diaryfname = [diaryfname0 int2str(i)];
end
diary(diaryfname)


if ~strcmp(theory,'dft')
  theorystr = ['_' theory];
else
  theorystr = '';
end
workdbname=[CD.dbdir filesep 'r' int2str(moltype) '_g' theorystr '.mat']    %#ok
load(workdbname,'workdb');
recnum=numel(workdb);

decssize = numel(workdb(1).prop.sdesc);
desc=char(recnum,1:decssize);
sdescnew=cell(recnum,decssize);
for i=1:recnum
    desc(i,1:decssize) = workdb(i).prop.sdesc(1:decssize); %obtained conformations
    sdescnew(i) = {workdb(i).prop.sdescnew};
end

trieddesc = {workdb(:).desc}; %conformations are tried to obtain

for i=1:recnum
%    delete(gca);
    ms0=workdb(i);

    descstr=desc(i,:);
%    if ~( isempty(ms0.new) | ms0.new=='Y')
    if ms0.new~='Y'
        if ms0.new=='B'
            [descstr ' beta ribose conformation: processing discarded'] %#ok
        end
        continue
    end

    ipC1 = ms0.ind(find(find(strcmp(pind.labels,'pC1'))==ms0.pind)); %#ok
%    ipC2 = ms0.ind(find(find(strcmp(pind.labels,'pC2'))==ms0.pind)); 
    ipC3 = ms0.ind(find(find(strcmp(pind.labels,'pC3'))==ms0.pind)); %#ok
    ipC4 = ms0.ind(find(find(strcmp(pind.labels,'pC4'))==ms0.pind)); %#ok
    ipC5 = ms0.ind(find(find(strcmp(pind.labels,'pC5'))==ms0.pind)); %#ok
    ipO2 = ms0.ind(find(find(strcmp(pind.labels,'pO2'))==ms0.pind)); %#ok
    ipO3 = ms0.ind(find(find(strcmp(pind.labels,'pO3'))==ms0.pind)); %#ok
    ipO4 = ms0.ind(find(find(strcmp(pind.labels,'pO4'))==ms0.pind)); %#ok
    ipO5 = ms0.ind(find(find(strcmp(pind.labels,'pO5'))==ms0.pind)); %#ok
    ipH12= ms0.ind(find(find(strcmp(pind.labels,'pH12'))==ms0.pind)); %#ok
    ipH21= ms0.ind(find(find(strcmp(pind.labels,'pH21'))==ms0.pind)); %#ok
    ipH22= ms0.ind(find(find(strcmp(pind.labels,'pH22'))==ms0.pind)); %#ok
    ipH31= ms0.ind(find(find(strcmp(pind.labels,'pH31'))==ms0.pind)); %#ok
    ipH32= ms0.ind(find(find(strcmp(pind.labels,'pH32'))==ms0.pind)); %#ok
    ipH53= ms0.ind(find(find(strcmp(pind.labels,'pH53'))==ms0.pind)); %#ok
    ibN1= ms0.ind(find(find(strcmp(pind.labels,'bN1'))==ms0.pind)); %#ok
    ibN9= ms0.ind(find(find(strcmp(pind.labels,'bN9'))==ms0.pind)); %#ok

    if moltype==7 || moltype==14 %only for non 2'deoxy molecyles
      ipO2 = ms0.ind(find(find(strcmp(pind.labels,'pO2'))==ms0.pind)); %#ok
    end
    ipC2 = ms0.ind(find(find(strcmp(pind.labels,'pC2'))==ms0.pind)); %#ok

    disp(['Processing : ' descstr])
    order=[];
    ms0=createzmt(ms0,order);
    ms1=zmt2xyz(ms0);  %for testing purpose

%    [xxx,order]=sortrows(ms1.pind);
    order=[]; %atom order is now defined by createbondchain
    curind2rotate=1;


%rotates around tgamma 
    plotmol(ms1);
    curind2rotate=curind2rotate+1;
    descstr=rotatesdesc(descstr,curind2rotate);
    [iprev,inext] = checksimdesc(descstr,desc);
    foundintrieddesc=strcmpcellar(trieddesc,descstr);
    opt=[strcmpar(desc,descstr),numel(foundintrieddesc),iprev,inext];
    if ~any(opt)
       ind = findallrotates(ms0,ipC5,ipC4,order);
       [ms1,errlev] = rotaroundbond(ms0,ind,120);
       ms1.desc = strcat(workname,'_',descstr,'.',ms1.prop.sdesc,'-1000');
       if errlev
         disp([lastwarn ': skipped']);
       elseif adist(ms1,ipH22,ipH32)<minHHdist
         disp([descstr ' ' num2str(adist(ms1,ipH22,ipH32))]);
       else 
         disp(descstr)
         plotmol(ms1);
         savemol(odir,ms1,0); 
%         savemolnw(odir,ms1,4); 
         savemolgs(odir,ms1,5,order,gtemplname); 
         desc(end+1,:)=descstr;%#ok
       end
    else
       if opt(1), disp1=desc(opt(1),:); else disp1='0'; end;
       if opt(2), disp2=desc(foundintrieddesc,:); else disp2='0'; end;
       if opt(3), disp3=desc(opt(3),:); else disp3='0'; end;
       if opt(4), disp4=desc(opt(4),:); else disp4='0'; end;
       disp(['skipped gamma+120 because of opt=[' disp1 ' ' disp2 ' ' disp3 ' ' disp4 ']'])
    end
    descstr=rotatesdesc(descstr,curind2rotate);
    [iprev,inext] = checksimdesc(descstr,desc);
    foundintrieddesc=strcmpcellar(trieddesc,descstr);
    opt=[strcmpar(desc,descstr),numel(foundintrieddesc),iprev,inext];
    if ~any(opt)
       ind = findallrotates(ms0,ipC5,ipC4,order);
       [ms1,errlev] = rotaroundbond(ms0,ind,240);
       ms1.desc = strcat(workname,'_',descstr,'.',ms1.prop.sdesc,'-2000');
       if errlev
         disp([lastwarn ': skipped']);
       elseif adist(ms1,ipH22,ipH32)<minHHdist
         disp([descstr ' ' num2str(adist(ms1,ipH22,ipH32))]);
       else 
         disp(descstr)
         plotmol(ms1);
         savemol(odir,ms1,0); 
%         savemolnw(odir,ms1,4); 
         savemolgs(odir,ms1,5,order,gtemplname); 
         desc(end+1,:)=descstr;%#ok
       end 
    else
       if opt(1), disp1=desc(opt(1),:); else disp1='0'; end;
       if opt(2), disp2=desc(foundintrieddesc,:); else disp2='0'; end;
       if opt(3), disp3=desc(opt(3),:); else disp3='0'; end;
       if opt(4), disp4=desc(opt(4),:); else disp4='0'; end;
       disp(['skipped gamma+240 because of opt=[' disp1 ' ' disp2 ' ' disp3 ' ' disp4 ']'])
    end
    descstr=rotatesdesc(descstr,curind2rotate);
    
%rotates around tbeta 
    curind2rotate=curind2rotate+1;
    descstr=rotatesdesc(descstr,curind2rotate);
    [iprev,inext] = checksimdesc(descstr,desc);
    foundintrieddesc=strcmpcellar(trieddesc,descstr);
    opt=[strcmpar(desc,descstr),numel(foundintrieddesc),iprev,inext];
    if ~any(opt)
       ind = findallrotates(ms0,ipO5,ipC5,order);
       [ms1,errlev] = rotaroundbond(ms0,ind,120);
       ms1.desc = strcat(workname,'_',descstr,'.',ms1.prop.sdesc,'-0100');
       if errlev
         disp([lastwarn ': skipped']);
       elseif adist(ms1,ipH22,ipH32)<minHHdist
         disp([descstr ' ' num2str(adist(ms1,ipH22,ipH32))]);
       else 
         disp(descstr)
         plotmol(ms1);
         savemol(odir,ms1,0); 
%         savemolnw(odir,ms1,4); 
         savemolgs(odir,ms1,5,order,gtemplname); 
         desc(end+1,:)=descstr;%#ok
       end
    else
       if opt(1), disp1=desc(opt(1),:); else disp1='0'; end;
       if opt(2), disp2=desc(foundintrieddesc,:); else disp2='0'; end;
       if opt(3), disp3=desc(opt(3),:); else disp3='0'; end;
       if opt(4), disp4=desc(opt(4),:); else disp4='0'; end;
       disp(['skipped beta+120 because of opt=[' disp1 ' ' disp2 ' ' disp3 ' ' disp4 ']'])
    end
    descstr=rotatesdesc(descstr,curind2rotate);
    [iprev,inext] = checksimdesc(descstr,desc);
    foundintrieddesc=strcmpcellar(trieddesc,descstr);
    opt=[strcmpar(desc,descstr),numel(foundintrieddesc),iprev,inext];
    if ~any(opt)
       ind = findallrotates(ms0,ipO5,ipC5,order);
       [ms1,errlev] = rotaroundbond(ms0,ind,240);
       ms1.desc = strcat(workname,'_',descstr,'.',ms1.prop.sdesc,'-0200');
       if errlev
         disp([lastwarn ': skipped']);
       elseif adist(ms1,ipH22,ipH32)<minHHdist
         disp([descstr ' ' num2str(adist(ms1,ipH22,ipH32))]);
       else 
         disp(descstr)
         plotmol(ms1);
         savemol(odir,ms1,0); 
%         savemolnw(odir,ms1,4); 
         savemolgs(odir,ms1,5,order,gtemplname); 
         desc(end+1,:)=descstr;%#ok
       end
    else
       if opt(1), disp1=desc(opt(1),:); else disp1='0'; end;
       if opt(2), disp2=desc(foundintrieddesc,:); else disp2='0'; end;
       if opt(3), disp3=desc(opt(3),:); else disp3='0'; end;
       if opt(4), disp4=desc(opt(4),:); else disp4='0'; end;
       disp(['skipped beta+240 because of opt=[' disp1 ' ' disp2 ' ' disp3 ' ' disp4 ']'])
    end
    descstr=rotatesdesc(descstr,curind2rotate);

%rotates around tepsilon 
    curind2rotate=curind2rotate+1;
    descstr=rotatesdesc(descstr,curind2rotate);
    [iprev,inext] = checksimdesc(descstr,desc);
    foundintrieddesc=strcmpcellar(trieddesc,descstr);
    opt=[strcmpar(desc,descstr),numel(foundintrieddesc),iprev,inext];
    if ~any(opt)
       ind = findallrotates(ms0,ipO3,ipC3,order);
       [ms1,errlev] = rotaroundbond(ms0,ind,120);
       ms1.desc = strcat(workname,'_',descstr,'.',ms1.prop.sdesc,'-0010');
       if errlev
         disp([lastwarn ': skipped']);
       elseif adist(ms1,ipH22,ipH32)<minHHdist
         disp([descstr ' ' num2str(adist(ms1,ipH22,ipH32))]);
       else 
         disp(descstr)
         plotmol(ms1);
         savemol(odir,ms1,0); 
%         savemolnw(odir,ms1,4); 
         savemolgs(odir,ms1,5,order,gtemplname); 
         desc(end+1,:)=descstr;%#ok
       end
    else
       if opt(1), disp1=desc(opt(1),:); else disp1='0'; end;
       if opt(2), disp2=desc(foundintrieddesc,:); else disp2='0'; end;
       if opt(3), disp3=desc(opt(3),:); else disp3='0'; end;
       if opt(4), disp4=desc(opt(4),:); else disp4='0'; end;
       disp(['skipped epsilon+120 because of opt=[' disp1 ' ' disp2 ' ' disp3 ' ' disp4 ']'])
    end
    descstr=rotatesdesc(descstr,curind2rotate);
    [iprev,inext] = checksimdesc(descstr,desc);
    foundintrieddesc=strcmpcellar(trieddesc,descstr);
    opt=[strcmpar(desc,descstr),numel(foundintrieddesc),iprev,inext];
    if ~any(opt)
%keyboard
       ind = findallrotates(ms0,ipO3,ipC3,order);
       [ms1,errlev] = rotaroundbond(ms0,ind,240);
       ms1.desc = strcat(workname,'_',descstr,'.',ms1.prop.sdesc,'-0020');
       if errlev
         disp([lastwarn ': skipped']);
       elseif adist(ms1,ipH22,ipH32)<minHHdist
         disp([descstr ' ' num2str(adist(ms1,ipH22,ipH32))]);
       else 
         disp(descstr)
         plotmol(ms1);
         savemol(odir,ms1,0); 
%         savemolnw(odir,ms1,4); 
         savemolgs(odir,ms1,5,order,gtemplname); 
         desc(end+1,:)=descstr;%#ok
       end
    else
       if opt(1), disp1=desc(opt(1),:); else disp1='0'; end;
       if opt(2), disp2=desc(foundintrieddesc,:); else disp2='0'; end;
       if opt(3), disp3=desc(opt(3),:); else disp3='0'; end;
       if opt(4), disp4=desc(opt(4),:); else disp4='0'; end;
       disp(['skipped epsilon+240 because of opt=[' disp1 ' ' disp2 ' ' disp3 ' ' disp4 ']'])
    end
    descstr=rotatesdesc(descstr,curind2rotate);

%rotates around teta
    if any(moltype==[7,14]) || (moltype>100 && mod(moltype,100)==40) %ribo molecules
    curind2rotate=curind2rotate+1;
        descstr=rotatesdesc(descstr,curind2rotate);
        [iprev,inext] = checksimdesc(descstr,desc);
        foundintrieddesc=strcmpcellar(trieddesc,descstr);
        opt=[strcmpar(desc,descstr),numel(foundintrieddesc),iprev,inext];
        if ~any(opt)
           ind = findallrotates(ms0,ipO2,ipC2,order);
           [ms1,errlev] = rotaroundbond(ms0,ind,120);
           ms1.desc = strcat(workname,'_',descstr,'.',ms1.prop.sdesc,'-0001');
           if errlev
             disp([lastwarn ': skipped']);
           elseif adist(ms1,ipH22,ipH32)<minHHdist
             disp([descstr ' ' num2str(adist(ms1,ipH22,ipH32))]);
           else 
             disp(descstr)
             plotmol(ms1);
             savemol(odir,ms1,0); 
    %         savemolnw(odir,ms1,4); 
             savemolgs(odir,ms1,5,order,gtemplname); 
             desc(end+1,:)=descstr;%#ok
           end
        else
        if opt(1), disp1=desc(opt(1),:); else disp1='0'; end;
           if opt(2), disp2=desc(foundintrieddesc,:); else disp2='0'; end;
           if opt(3), disp3=desc(opt(3),:); else disp3='0'; end;
           if opt(4), disp4=desc(opt(4),:); else disp4='0'; end;
           disp(['skipped eta+120 because of opt=[' disp1 ' ' disp2 ' ' disp3 ' ' disp4 ']'])
        end
        descstr=rotatesdesc(descstr,curind2rotate);
        [iprev,inext] = checksimdesc(descstr,desc);
        foundintrieddesc=strcmpcellar(trieddesc,descstr);
        opt=[strcmpar(desc,descstr),numel(foundintrieddesc),iprev,inext];
        if ~any(opt)
           ind = findallrotates(ms0,ipO2,ipC2,order);
           [ms1,errlev] = rotaroundbond(ms0,ind,240);
           ms1.desc = strcat(workname,'_',descstr,'.',ms1.prop.sdesc,'-0002');
           if errlev
             disp([lastwarn ': skipped']);
           elseif adist(ms1,ipH22,ipH32)<minHHdist
             disp([descstr ' ' num2str(adist(ms1,ipH22,ipH32))]);
           else 
             disp(descstr)
             plotmol(ms1);
             savemol(odir,ms1,0); 
    %         savemolnw(odir,ms1,4); 
             savemolgs(odir,ms1,5,order,gtemplname); 
             desc(end+1,:)=descstr; %#ok
           end
        else
           if opt(1), disp1=desc(opt(1),:); else disp1='0'; end;
           if opt(2), disp2=desc(foundintrieddesc,:); else disp2='0'; end;
           if opt(3), disp3=desc(opt(3),:); else disp3='0'; end;
           if opt(4), disp4=desc(opt(4),:); else disp4='0'; end;
           disp(['skipped eta+240 because of opt=[' disp1 ' ' disp2 ' ' disp3 ' ' disp4 ']'])
        end
        descstr=rotatesdesc(descstr,curind2rotate);
    elseif i==1
    warning('Molecule detected as deoxyribo. Rotates round C2''O2'' are not performed.'); %#ok
    end


    if any(moltype~=[7,8]) %nucleosides/nucleotides

        begchiind='BCTV';
        endchiind='CBVT';
        chi2rotate=[60 -60 30 -30];

        chiconfind=find(descstr(end)==begchiind);

        descstr(end)=endchiind(chiconfind);

        [iprev,inext] = checksimdesc(descstr,desc);
        foundintrieddesc=strcmpcellar(trieddesc,descstr);
        opt=[strcmpar(desc,descstr),numel(foundintrieddesc),iprev,inext];
        if ~any(opt)
            if ~isempty(strcmpcellar(pind.labels(ms0.pind),'bN9')) %N1 atom is present in all nucleosides, so search of N9 must be done at first
                ind = findallrotates(ms0,ipC1,ibN9,order);
            elseif ~isempty(strcmpcellar(pind.labels(ms0.pind),'bN1'))
                ind = findallrotates(ms0,ipC1,ibN1,order);
            else
                error('Atoms of glycosidic bond are not found.');
            end

            [ms1,errlev] = rotaroundbond(ms0,ind,chi2rotate(chiconfind));

            ms1.desc = strcat(workname,'_',descstr,'.',ms1.prop.sdesc,'-0high');
            if errlev
              disp([lastwarn ': skipped']);
            elseif adist(ms1,ipH22,ipH32)<minHHdist
              disp([descstr ' ' num2str(adist(ms1,ipH22,ipH32))]);
            else 
              disp(descstr)
              plotmol(ms1);
              savemol(odir,ms1,0); 
     %         savemolnw(odir,ms1,4); 
              savemolgs(odir,ms1,5,order,gtemplname); 
              desc(end+1,:)=descstr;%#ok
            end
        else
           if opt(1), disp1=desc(opt(1),:); else disp1='0'; end;
           if opt(2), disp2=desc(foundintrieddesc,:); else disp2='0'; end;
           if opt(3), disp3=desc(opt(3),:); else disp3='0'; end;
           if opt(4), disp4=desc(opt(4),:); else disp4='0'; end;
           disp(['skipped chi change because of opt=[' disp1 ' ' disp2 ' ' disp3 ' ' disp4 ']'])
        end

        descstr(end)=endchiind(chiconfind);
    elseif i==1
        warning('Molecule seems to be without tchi torsion. Rotates round glycosidic bond are not performed.'); %#ok
    end


%    'Press any key.'
%    pause
end

disp(['Total ' int2str(size(desc,1)-size(trieddesc,2)) ' jobs created.']);
diary off

