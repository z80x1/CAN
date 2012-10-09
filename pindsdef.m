%
% Version 1.0    
% Last modified  R O Zhurakivsky 2009-01-09 
% Created        R O Zhurakivsky 2005-09-?

format compact

%pretty atoms indexes :)
global pind

%Attention!!! New lables must be added only to the end of array!
pind.labels = [{'pC1'};{'pC2'};{'pC3'};{'pC4'};{'pC5'};...
{'pO4'};{'pO2'};{'pO3'};{'pO5'};{'pH11'};{'pH12'};{'pH21'};{'pH22'};...
{'pH31'};{'pH32'};{'pH4'};{'pH51'};{'pH52'};{'pH53'};... %#14-19
{'bN1'};{'bC2'};{'bN3'};{'bC4'};{'bC5'};{'bN6'};{'bN'};{'bO2'};{'bH5'};{'bH41'};{'bH42'};... % 6azaCyd  #21-30
{'bC6'};{'bH6'};... % Cyd
{'bH3'};{'bO4'};... % Urd, Thd & Cyd imino
{'bCC5'};{'bH51'};{'bH52'};{'bH53'};... % Thd #35-38
{'bN7'};{'bC8'};{'bN9'};{'bH2'};{'bCN6'};{'bH61'};{'bH62'};{'bH8'};... % Ade #39-46
{'bH1'};{'bCN2'};{'bH21'};{'bH22'};{'bCO6'};{'bH7'};{'bCO8'};... % Guo #47-53
{'bN5'}; %5aza structures atom
{'bBr5'}; %5bromine structures atom
{'bF6'}; %6flourine structures atom
{'bN8'}; %8aza structures atom
{'pP5'};{'pOP1'};{'pOP2'};{'pOP3'};{'pHP1'};{'pHP2'};{'xN'};{'xH1'};{'xH2'};{'xH3'};{'pP3'};... %nucleotides atoms  #58-68
{'bH4'};... %Urd, Thd lactim
{'pS4'};... %tionucleosides
{'pH42'};... %moltypes 2,3
{'bNC7'};{'bH71'};{'bH72'};{'bH73'};... %moltype 353 #72-75
{'bCl5'};{'bF5'};...
{'pNa'};...  %moltype 950 #78
{'pN31'};{'pN32'};{'pN33'};...  %moltype 413 AZT
];

pind.ind = 1:numel(pind.labels);

%pind.ind(end+1)=numel(pind.ind)+1;
pind.ind(end+1)=999;
pind.labels(end+1)={'dummy'};
