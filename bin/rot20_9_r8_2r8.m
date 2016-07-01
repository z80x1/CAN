%add azapirimidine ring to selected sugar conformations and make anti and syn conformations
%
% Version 1.0    
% Last modified  gator
% Created        gator 2006-?-?

format compact
clear
global pind
global flPlot



atomsind
pindsdef
flPlot=0;

workdbname = [CD.dbdir filesep 'r8_g_dftV2_or.mat']
workname='r91001';
odir=[CD.xyzdir filesep workname];
moltype=910;
gtemplname = [workname '_templ.gjf']
fullgtemplname = [CD.templatesdir filesep gtemplname]

load(workdbname,'workdb')

dO5P = 1.5946;

dOP = (1.47566+1.48314)/2;

recnum=numel(workdb);
for i=1:recnum           %!!!!!!!!!!!!  
  if workdb(i).new~='Y'
      continue
  end

    ms0=workdb(i);
    
    ['Loaded conformation: ' ms0.prop.sdesc]
%plotmol(ms0,'r',0);

%    maxind=0;
%    mc0.atomnum = length(mc0.labels);
%    mc0.ind = ((maxind+1):(maxind+mc0.atomnum))';
%    maxind=maxind+mc0.atomnum;

    indA  = ms0.ind(find(find(strcmp(pind.labels,'pO5')) == ms0.pind)); 
    indB = ms0.ind(find(find(strcmp(pind.labels,'pH53'))== ms0.pind)); 
    A=[ms0.x(indA) ms0.y(indA) ms0.z(indA)];
    B=[ms0.x(indB) ms0.y(indB) ms0.z(indB)];
    AB=B-A;

    ms0.pind(indB) = find(strcmp(pind.labels,'pP5'));
    ms0.labels(indB)={'P'};
    
    B=A+AB/norm(AB)*dO5P;
    ms0.x(indB)=B(1);
    ms0.y(indB)=B(2);
    ms0.z(indB)=B(3);

    %тетраедр AKLM , вершина А - задана (О53), центр O - P5, вектор EA - паралельний
    % вектору BK - C5'O5'
    %знайти координати K.L.M
    indF = ms0.ind(find(find(strcmp(pind.labels,'pC5'))==ms0.pind)); 
    F=[ms0.x(indF) ms0.y(indF) ms0.z(indF)];
    
    O=B;
    AF=A-F;  % AF=A-F было когда стартовали
    K=O+AF/norm(AF)*dO5P;
    
   

    L = fsolve(@(p) tetraider(p,A,B,K),[B(1);B(2);B(3)]);
    L = L';
    
    
    M = fsolve(@(p) tetraider2(p,L,A,K),[B(1);B(2);B(3)]);

    M = M';
  
    indO4 = ms0.ind(find(find(strcmp(pind.labels,'pO4'))== ms0.pind)); 
    Y=[ms0.x(indO4) ms0.y(indO4) ms0.z(indO4)];
     
            if norm(Y-M) < norm(Y-L)
               r = M;
               M = L;
               L = r; 
               disp('ZAMENA!') 
            end

     BM=M-B;
     M=B+BM/norm(BM)*dOP;
     BK=K-B;
     K=B+BK/norm(BK)*dOP;
     BL=L-B;      
     L=B+BL/norm(BL)*dO5P;       
            
            
    ms0.atomnum=ms0.atomnum+3;
    ms0.x=[ms0.x; K(1); L(1); M(1)];
    ms0.y=[ms0.y; K(2); L(2); M(2)];
    ms0.z=[ms0.z; K(3); L(3); M(3)];
    ms0.labels=[ms0.labels; {'O'}; {'O'}; {'O'}];
    ms0.ind=[ms0.ind; max(ms0.ind)+1; max(ms0.ind)+2; max(ms0.ind)+3];
    ms0.pind(end+1)=find(strcmp(pind.labels,'pOP1')); %K
    ms0.pind(end+1)=find(strcmp(pind.labels,'pOP2')); %L
    ms0.pind(end+1)=find(strcmp(pind.labels,'pOP3')); %M
    
%      plotmol(ms0,'r',1);
    
    indO32  = ms0.ind(find(find(strcmp(pind.labels,'pOP2')) == ms0.pind)); %L
    
    for j=1:recnum
        if workdb(j).new~='Y'
            continue
        end
        ms2=workdb(j);
       
        OM = M-O; % OM=M-O было когда стартовали
       
        indO3 = ms2.ind(find(find(strcmp(pind.labels,'pO3'))==ms2.pind)); 
        indH32 = ms2.ind(find(find(strcmp(pind.labels,'pH32'))==ms2.pind)); 
%        indh32= ms0.ind(find(find(strcmp(pind.labels,'pH32')) == ms0.pind)); 
    
        vect_mov = createvect(ms2,indO3,ms0,indO32);  %pPO3
        ms2.x = ms2.x + vect_mov(1);
        ms2.y = ms2.y + vect_mov(2);
        ms2.z = ms2.z + vect_mov(3);
    
        %vector H32'-O3'
        vect_ms2 = createvect(ms2,indH32,ms2,indO3);
        %rotates second residue so that H32'-O3' to be along P-OP3
        ms2 = rotvect3(OM, vect_ms2, ms2, indO3);
 
        ms2 = delatom( ms2, indH32 );  %H32

%        ind=find(ms0.ind==indO32); %true index
%        ms0_order= [ms0_order(1:ind-1); ms0_order(ind+1:end)];
        ms00 = delatom( ms0, indO32 ); %O3'
        
        [X,ms00_order]=sort(ms00.pind);
        [X,ms2_order]=sort(ms2.pind);
        ms2_order=ms2_order+max(ms00_order);

        mm=mergemol(ms00,ms2);
        mm=createbondtable(mm);    
    
    
%plotmol(mm,'r',0);
    
        
        indO41 = find(find(strcmp(pind.labels,'pO4'))== ms00.pind);
        indO42 = find(find(strcmp(pind.labels,'pO4'))== ms2.pind)+ms00.atomnum;
        Q=[ms0.x(indO41) ms0.y(indO41) ms0.z(indO41)];
        W=[mm.x(indO42) mm.y(indO42) mm.z(indO42)];
        regulating = norm(Q-W);
        
     
        indH = find(find(strcmp(pind.labels,'pH32'))== ms00.pind);
        indH2 = find(find(strcmp(pind.labels,'pH53'))== ms2.pind)+ms00.atomnum;
        Q2=[ms0.x(indH) ms0.y(indH) ms0.z(indH)];
        W2=[mm.x(indH2) mm.y(indH2) mm.z(indH2)];
        regulating2 = norm(Q2-W2);
        
        if regulating <= 3.4
            mm.desc = strcat('u',workname,'_',ms00.prop.sdesc,'_',ms2.prop.sdesc);
            route_add='opt=(GDIIS, maxcycle=50)';
           
        elseif regulating >= 3.54640 & regulating2 >= 2.86625 
            mm.desc = strcat('g',workname,'_',ms00.prop.sdesc,'_',ms2.prop.sdesc);           
            route_add='opt=(GDIIS)';
           
        else
            mm.desc = strcat(workname,'_',ms00.prop.sdesc,'_',ms2.prop.sdesc); 
            route_add='opt=(GDIIS, maxcycle=80)';
           
        end
           
   
        %        mcyd = identmol(mm,moltype); 

        if exist(odir)~=7
            mkdir(odir);
        end
        
        

        %savemol(odir,mm,0);
%        disp(['Written: ' mm.desc])
        order=[ms00_order; ms2_order];
       
        
        
     
        savemolgs(odir,mm,3,order,fullgtemplname,route_add); %Gaussian 
    
%      break
    end

%pause
%clf
%     break

end

