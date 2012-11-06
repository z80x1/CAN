classdef PlotsLimits < handle

    properties
      keys
      data
%      filename
    end
   
    % Class methods
    methods 

%     	function obj = PlotsLimits(f)
% 			obj.filename = f;
%         end 

        function res = addlimit(obj, key, data)
            obj.keys{end+1} = key;
            obj.data(end+1,:) = data;
            res = data;
        end 
        
        function res = updatelimit(obj, key, data)
            found = 0;
            for I = 1 : numel(obj.keys)
                if strcmp(obj.keys{I}, key)
                    found = 1;
                    obj.data(I,:) = data;
                end
            end
            if found
                res = data;
            else
                res = [];
            end
        end

        function res = getlimit(obj, key)
            res = [];
            for I = 1 : numel(obj.keys)
                if strcmp(obj.keys{I}, key)
                    res = obj.data(I,:);
                end
            end
        end 
        
        function res = correctlimit(obj, key, d)
%            res = [];
            fl_found = 0;
            for I = 1 : numel(obj.keys)
                if strcmp(obj.keys{I}, key)
                    fl_found = 1;
                    tmp = [d(1:2) min(obj.data(I,3), d(3)) max(obj.data(I,4), d(4))];
                    
                    res = updatelimit(obj, key, tmp);
                end
            end
            if ~fl_found
                res = addlimit(obj, key, d);
            end
        end 
        
%         function save(obj)
%             save(obj.filename,'obj');
%         end 

    	function disp(obj)
		     % DISP Display object in MATLAB syntax
	        disp(obj.data)
		end % disp

   end
end