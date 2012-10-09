function out = lr_translit(in)
%transliterate cyrillic Win1251 string 
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-03-05 
% Created        R O Zhurakivsky 2007-?-?

in=unicode2native(in,'');

%                   [À], [Á], [Â], [Ã], [Ä], [Å], [Æ], [Ç], [È], [É], [Ê], [Ë], [Ì], [Í], [Î], [Ï], [Ğ], [Ñ], [Ò], [Ó], [Ô], [Õ], [Ö], [×], [Ø], [Ù], [Ú], [Û], [Ü], [İ], [Ş], [ß]
%intable1 = hex2dec(['C0';'C1';'C2';'C3';'C4';'C5';'C6';'C7';'C8';'C9';'CA';'CB';'CC';'CD';'CE';'CF';'D0';'D1';'D2';'D3';'D4';'D5';'D6';'D7';'D8';'D9';'DA';'DB';'DC';'DD';'DE';'DF']);
%intable2 = hex2dec(['E0';'E1';'E2';'E3';'E4';'E5';'E6';'E7';'E8';'E9';'EA';'EB';'EC';'ED';'EE';'EF';'F0';'F1';'F2';'F3';'F4';'F5';'F6';'F7';'F8';'F9';'FA';'FB';'FC';'FD';'FE';'FF']);
outtable=        ([{'A'}; {'B'}; {'V'}; {'G'}; {'D'}; {'E'}; {'Zh'};{'Z'}; {'Y'}; {'J'}; {'K'}; {'L'}; {'M'}; {'N'}; {'O'}; {'P'};...
                   {'R'}; {'S'}; {'T'}; {'U'}; {'F'}; {'H'}; {'C'}; {'Ch'};{'Sh'};{'Sch'};{''}; {'Y'}; {''''};{'E'}; {'Yu'};{'Ya'};...
                   {'a'}; {'b'}; {'v'}; {'g'}; {'d'}; {'e'}; {'zh'};{'z'}; {'y'}; {'j'}; {'k'}; {'l'}; {'m'}; {'n'}; {'o'}; {'p'};...
                   {'r'}; {'s'}; {'t'}; {'u'}; {'f'}; {'h'}; {'c'}; {'ch'};{'sh'};{'sch'};{''}; {'y'}; {''''};{'e'}; {'yu'};{'ya'}]);

%                   [¥], [¨], [ª], [¯], [²], [³], [´], [¸], [º], [¿]
%intable3 = hex2dec(['A5';'A8';'AA';'AF';'B2';'B3';'B4';'B8';'BA';'BF']);
%outtable3=        (['G'; 'Yo';'Ye';'Yi';'I'; 'i'; 'g'; 'yo';'ye';'yi']);


len=length(in);
out=zeros(1,2*len);
j=1;
for i=1:len
if in(i)>164
  if in(i)>191
     ind=in(i)-191;
     jend=j+length(outtable{ind});
     out(j:jend-1)=outtable{ind};
     j=jend;
  elseif in(i)==165 %[¥]
     out(j)='G';
     j=j+1;
  elseif in(i)==168 %[¨]
     out(j:j+1)='Yo';
     j=j+2;
  elseif in(i)==171 %[ª]
     out(j:j+1)='Ye';
     j=j+2;
  elseif in(i)==175 %[¯]
     out(j:j+1)='Yi';
     j=j+2;
  elseif in(i)==178 %[²]
     out(j)='I';
     j=j+1;
  elseif in(i)==179 %[³]
     out(j)='i';
     j=j+1;
  elseif in(i)==180 %[´]
     out(j)='g';
     j=j+1;
  elseif in(i)==184 %[¸]
     out(j:j+1)='yo';
     j=j+2;
  elseif in(i)==186 %[º]
     out(j:j+1)='ye';
     j=j+2;
  elseif in(i)==191 %[¿]
     out(j:j+1)='yi';
     j=j+2;
  else
     out(j)=in(i);
     j=j+1;
  end
else
  out(j)=in(i);
  j=j+1;
end
end
out=out(1:j-1);
