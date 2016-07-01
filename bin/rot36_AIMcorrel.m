%rot36_AIMcorrel: analyze correlations for all Hbonds
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2009-04-25
% Created        R O Zhurakivsky 2009-04-24

tic
format compact

global pind
clear h;
clear HBs HBsgeom

atomsind
pindsdef

%---------------------------------
fl_analyzegeomcriteria=0 %#ok
fl_selectHBtypes=0 %#ok %display data only for selected Hbond types
Atype='C' %#ok A atom type for AH...B bond
Btype='N' %#ok B atom type for AH...B bond
fl_plot_AHB_AB=1; % E_AHB_AB or E_rho_deltarho
fl_plotmnk=0;
%---------------------------------

pcolor=[{'r'} {'g'} {'b'} {'k'} {'m'} {'c'} ...
        {[.8627 .0784 .2353]} {[1 .2706 0]} {[.5412 .1686 .8863]}  ...
        {[0 .5 0]} {[.5 0 0]} {[0 0 .5]} {[.5 .5 0]} {[.5 0 .5]} {[0 .5 .5]} {[.5 .5 .5]} ...
        {[.5 .1 .9]} {[.8 .2 .2]} {[.8 .8 .2]}...
        {[.9 .4 .9]} {[.2 .4 .6]} {[.6 .4 .6]} {[.6 .2 .2]} {[.8 .2 .8]} ...
        {[.2 .8 .8]} {[.2 .8 .2]} {[.2 .2 .8]} {[.4 .9 .1]} {[.1 .3 .6]} {'y'} {[.75 .75 .75]} {[.2745 .5098 .7059]}];
psign='x+*osdv^<>ph';
plotattr=[{'ro'} {'go'} {'bo'} {'ko'} {'mo'} {'co'} {'r+'} {'g+'} {'b+'} {'k+'} {'m+'} {'c+'}];

HBall=0;
HBallgeom=0;

HBtypes=zeros(0,3);
h=[]; %plot handles array
legs=[]; %legends array
legs_pinds=[]; %legends array
legs_h=[];

fl_firstmol=1;
for m_ind=1:numel(db)

    workdb=db(m_ind).workdb;
    recnum=numel(workdb);
    
    HB=[];
    HB.pinds=[];
    indHB=0;

    HB4=[];
    HB4.pinds=[];
    indHB4=0;
    
    HBgeom=[];
    HBgeom.distAB=[];
    indHBgeom=0;
    
    for i_rec=1:recnum
        ms0 = workdb(i_rec);
        HBnum=size(ms0.AIM.pinds,1);
       
       for i_HB=1:HBnum

           pinds = ms0.AIM.pinds(i_HB,:);
           if numel(pinds)<3 || ~pinds(3) %ignore 'VdV contacts' 
               continue
           end
           if numel(pinds)==4 && pinds(4) % H...H bonds
               indHB4=indHB4+1;
               HB4.pinds(indHB4,1:numel(pinds))=pinds;
               atomH1ind = find(ms0.pind==pinds(2));
               atomH2ind = find(ms0.pind==pinds(3));
               
               HB4.atomH1.x(indHB4) = ms0.x(atomH1ind);
               HB4.atomH1.y(indHB4) = ms0.y(atomH1ind);
               HB4.atomH1.z(indHB4) = ms0.z(atomH1ind);
               HB4.atomH2.x(indHB4) = ms0.x(atomH2ind);
               HB4.atomH2.y(indHB4) = ms0.y(atomH2ind);
               HB4.atomH2.z(indHB4) = ms0.z(atomH2ind);

               HB4.distHH(indHB4) = adist(ms0,atomH1ind,atomH2ind);
               continue
           end
           atomAind = find(ms0.pind==pinds(1));
           atomHind = find(ms0.pind==pinds(2));
           atomBind = find(ms0.pind==pinds(3));

           if fl_selectHBtypes && (ms0.labels{atomAind}~=Atype || ms0.labels{atomBind}~=Btype)
               continue
           end

           indHB=indHB+1;
           HB.pinds(indHB,1:numel(pinds))=pinds;

           HBtype = find(all(HBtypes==repmat(pinds(1:3),size(HBtypes,1),1),2));
           if HBtype
               HB.type(indHB)=HBtype;
           else
               HBtypes(end+1,1:3)=pinds(1:3);
               HB.type(indHB)=size(HBtypes,1);
           end

           HB.atomA.labels(indHB) = ms0.labels(atomAind);

           HB.atomH.x(indHB) = ms0.x(atomHind);
           HB.atomH.y(indHB) = ms0.y(atomHind);
           HB.atomH.z(indHB) = ms0.z(atomHind);

           HB.atomB.x(indHB) = ms0.x(atomBind);
           HB.atomB.y(indHB) = ms0.y(atomBind);
           HB.atomB.z(indHB) = ms0.z(atomBind);
           HB.atomB.labels(indHB) = ms0.labels(atomBind);

           HB.distAB(indHB) = adist(ms0,atomAind,atomBind);
           HB.distHB(indHB) = adist(ms0,atomHind,atomBind);
           HB.angAHB(indHB) = valang(ms0,atomAind,atomHind,atomBind);
            
           HB.rho(indHB) = ms0.AIM.ro(i_HB);
           if isfield(ms0.AIM,'deltaroAIM2000')
               HB.deltarho(indHB) = -4*ms0.AIM.deltaroAIM2000(i_HB);
           elseif isfield(ms0.AIM,'DelSqRho')
               HB.deltarho(indHB) = ms0.AIM.DelSqRho(i_HB);
           end

           if isfield(ms0.AIM,'V')
               HB.energy(indHB) = abs(0.5*ms0.AIM.V(i_HB)*CC.encoef);
           else
               HB.energy(indHB) = 0;
           end

           HB.G(indHB) = ms0.AIM.G(i_HB);

       end
       
if fl_analyzegeomcriteria
        %AH...B H-bonds searching by geometry criteria 
        %maybe searching at first for H 
        %HbondHind, HbondO1ind, HbondO2ind - hardindexes
        HbondHind=[];
        HbondO1ind=[];
        HbondO2ind=[];
        for i=setdiff(1:ms0.atomnum,strcmpcellar(ms0.labels,'H')) %all nonhydrogen atoms

           if fl_selectHBtypes && (ms0.labels{i}~=Atype) %filter by A atom type for AH...B bond
               continue
           end

            
            iAlink = [ms0.btA(ms0.btB==i) ms0.btB(ms0.btA==i)];
            for iH = iAlink  % find hydrogens connected to this atom
                if ms0.labels{iH}=='H'
                  for j=setdiff(1:ms0.atomnum,strcmpcellar(ms0.labels,'H')) 

                    if fl_selectHBtypes && (ms0.labels{j}~=Btype) %filter by B atom type for AH...B bond
                         continue
                    end
                      
                    if j~=i
                      HbondABdist=adist(ms0,i,j);
                      HbondHBdist=adist(ms0,iH,j);
                      HbondAHBang =valang(ms0,i,iH,j);
                      if HbondABdist<CG.lcritABdist && HbondHBdist<CG.lcritHBdist && HbondAHBang>CG.lcritAHBang %verify H bond geometry criteria
                        indHBgeom=indHBgeom+1;
                          
                        HBgeom.distAB(indHBgeom)=HbondABdist;
                        HBgeom.angAHB(indHBgeom)=HbondAHBang;
                        HBgeom.distHB(indHBgeom)=HbondHBdist;
                        HBgeom.distABdesc(indHBgeom)={['r' pind.labels{ms0.pind(i)} pind.labels{ms0.pind(j)}]};
                        HBgeom.angAHBdesc(indHBgeom)={['a' pind.labels{ms0.pind(i)} pind.labels{ms0.pind(iH)} pind.labels{ms0.pind(j)}]};
                        HbondHind(indHBgeom)=iH;
                        HbondO1ind(indHBgeom)=i;
                        HbondO2ind(indHBgeom)=j;

                        pinds = [ms0.pind(i) ms0.pind(iH) ms0.pind(j)];
                        HBgeom.pinds(indHBgeom,1:numel(pinds))=pinds;
                        HBtype = find(all(HBtypes==repmat(pinds(1:3),size(HBtypes,1),1),2));
                        if HBtype
                           HBgeom.type(indHBgeom)=HBtype;
                        else
                           HBtypes(end+1,1:3)=pinds(1:3);
                           HBgeom.type(indHBgeom)=size(HBtypes,1);
                        end

                      end
                    end
                  end
                end
            end
        end
end %fl_analyzegeomcriteria        
        
    end %1:recnum
    
    if ~indHB
        warning(['r' int2str(db(m_ind).moltype) ': zero number of bonds']);
        db(m_ind).HB=[];   
        db(m_ind).HB4=[];
        db(m_ind).HBgeom=[];
        continue;
    end
        
    db(m_ind).HB=HB;
    db(m_ind).HB4=HB4;
    db(m_ind).HBgeom=HBgeom;
    
%-------------- plotting
    
%    HB.atomA.labels
%    HB.atomB.labels
    if fl_firstmol
        figure
    end

    fl_mol_in_legend=0;
    HBtypes_num=size(HBtypes,1);
for Htype_ind=1:HBtypes_num
    inds = find(HB.type==Htype_ind);
    if isempty(inds), continue, end;

    if fl_plot_AHB_AB
      h(end+1)=plot3(HB.distAB(inds),HB.angAHB(inds),HB.energy(inds), [psign(mod(m_ind,numel(psign))+1)], 'MarkerSize', 5, 'Color', pcolor{mod(Htype_ind,numel(pcolor))+1} );
    else
      h(end+1)=plot3(HB.rho(inds),HB.deltarho(inds),HB.energy(inds), [psign(mod(m_ind,numel(psign))+1)], 'MarkerSize', 5, 'Color', pcolor{mod(Htype_ind,numel(pcolor))+1});
    end

    if isempty(legs_pinds) || isempty(find(all(repmat(HBtypes(Htype_ind,:),size(legs_pinds,1),1)==legs_pinds,2))) ...
        || (Htype_ind==HBtypes_num && ~fl_mol_in_legend) %add to legend only not used Hbond types and at least 1 legend per molecule
        legs_h(end+1)=h(end);
        legs{end+1}=['r' int2str(db(m_ind).moltype) ' : ' [pind.labels{ HBtypes(Htype_ind,:)}] ];
        legs_pinds(end+1,1:3)=HBtypes(Htype_ind,:);
        fl_mol_in_legend=1;
    end
    
end    

%    h(m_ind)=plot3(HB.distAB,HB.angAHB,HB.energy, plotattr{mod(m_ind,numel(plotattr))+1} );
%    legs{end+1}=['r' int2str(db(m_ind).moltype)  ];
%        plot3(HB.rho,HB.deltarho,HB.energy,plotattr{m_ind});

    HBall = HBall + indHB;
    HBallgeom = HBallgeom + indHBgeom;

    if fl_firstmol
        hold on
        grid on
        if fl_plot_AHB_AB
          xlabel('distance A...B, A');
          ylabel('angle AH...B, degree');
        else
          xlabel('\rho, a.e.');
          ylabel('\Delta \rho, a.e.');
        end
        zlabel('bond energy, kcal/mol');
        fl_firstmol=0;
    end

    disp(['r' int2str(db(m_ind).moltype) ': Total ' int2str(HBall) ' of H-bonds, current molecule - ' int2str(indHB) ', geom criteria: ' int2str(HBallgeom) ' - ' int2str(indHBgeom)]);
end %1:numel(db)

bondtypestr='';
if fl_selectHBtypes
    bondtypestr=[Atype 'H...' Btype ' '];
end
title([' Total ' int2str(HBall) ' of ' bondtypestr 'H-bonds']);

%inds=find(ishandle(h));
%hl=legend(h(inds), legs(inds),'Location','NEO');
inds=find(ishandle(legs_h));
hl=legend(legs_h(inds), legs(inds),'Location','NEO');
set(hl,'FontSize',6 );


HBs=[db.HB];
%HB4s=[db.HB4];
HBsgeom=[db.HBgeom];

max([HBs.distAB]) % 9mol CH...O 3.7244      %26mol all bonds 3.9538
min([HBs.angAHB]) % 9mol CH...O 100.7639    %26mol all bonds 98.5033
max([HBs.distHB]) % 9mol CH...O 2.9514      %26mol all bonds 3.1716


%output table with persentage of AIM bonds to geometry bond for each bond type
if fl_analyzegeomcriteria
    HBs_type=[HBs.type];
    trueHBsgeom_ind = [HBsgeom.distAB]<max([HBs.distAB]) & [HBsgeom.angAHB]>min([HBs.angAHB]) & [HBsgeom.distHB]<max([HBs.distHB]);
    HBsgeom_type=[HBsgeom.type];
    HBsgeom_type=HBsgeom_type(trueHBsgeom_ind);
    [res,I]=sortrows(HBtypes);
    for Htype_ind=I'
     
        HBAIMinds = find(HBs_type==Htype_ind);
        HBgeominds = find(HBsgeom_type==Htype_ind);

        disp(sprintf('%s: \tAIM: %d , \tgeom: %d , \t%%: %4.3f', [pind.labels{ HBtypes(Htype_ind,:)}],numel(HBAIMinds),numel(HBgeominds),numel(HBAIMinds)/numel(HBgeominds)))
    end    
end


if fl_plotmnk
  [a,b,err]=mnk([HBs.rho], [HBs.energy])
  [rho,ix]=sort([HBs.rho]);

  xxx=[0.002 0.024];
  plot3(xxx,[0 0],a+b*xxx,'r')

if 0
  energy=[HBs.energy];
  energy2=energy(ix);
  figure
  plot(rho,abs(a+b*rho-energy2)./energy2,'r')
  grid on
  xlabel('\rho');
  ylabel('abs ( a + b * \rho - E ) / E');
end

end


hold off
toc


if 0
  energy=[HBs.energy];
  sum(energy<0.6)/numel(energy); %#ok
end
