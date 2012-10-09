classdef TorsionStat
   % Documentation example
   % A value class that implements a data type for polynonials
   % See Implementing a Class for Polynomials in the
   % MATLAB documentation for more information.
   
   properties
      data
   end
   
   % Class methods
   methods

		function obj = TorsionStat(data)
			obj.data = data;
        end 

        function str = char(obj)
          s{1} = ['min=' num2str(min(obj.data),4)];
          s{2} = [', max=' num2str(max(obj.data),4)];
          s{3} = [', mean=' num2str(mean(obj.data),4)]; 
          s{4} = [', std=' num2str(std(obj.data),4)]; 
          str = [s{:}];
        end % char

		function disp(obj)
		     % DISP Display object in MATLAB syntax
		     c = char(obj);
		     if iscell(c)
		        disp(['     ' c{:}])
		     else
		        disp(c)
		     end
		end % disp
        
		function d = stat(obj)
            d = {min(obj.data) max(obj.data) mean(obj.data) std(obj.data)}';
		end % stat
        
%       function b = subsref(a,s)
%          % SUBSREF Implementing the following syntax: 
%          % obj([1 ...])
%          % obj.coef
%          % obj.plot
%          % out = obj.method(args)
%          % out = obj.method
%          switch s(1).type
%             case '()'
%                ind = s.subs{:};
%                b = a.polyval(ind);
%             case '.'
%                switch s(1).subs
%                   case 'coef'
%                      b = a.data;
%                   case 'stat'
%                      b = a.stat;
%                   otherwise
%                      if length(s)>1
%                         b = a.(s(1).subs)(s(2).subs{:});
%                      else
%                         b = a.(s.subs);
%                      end
%                end
%             otherwise
%                error('Specify value for x as obj(x)')
%          end
%       end % subsref      

	end % methods 
end % classdef




