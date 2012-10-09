function display_skip_reason(desc, str, opt, foundintrieddesc)
       if opt(1), disp1=desc(opt(1),:); else disp1='0'; end;
       if opt(2), disp2=desc(foundintrieddesc,:); else disp2='0'; end;
       for i=1:size(disp2,1)
           if i==1 
               str2 = disp2(i,:);
           else
               str2 = [str2 '/' disp2(i,:)];
           end
       end
       if opt(3), disp3=desc(opt(3),:); else disp3='0'; end;
       if opt(4), disp4=desc(opt(4),:); else disp4='0'; end;
       disp(['skipped ' str ' because of opt=[' disp1 ' ' str2 ' ' disp3 ' ' disp4 ']'])
end

