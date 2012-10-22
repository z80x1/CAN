function rotatefig(src,evnt)
%rotates current axes
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2012-10-15
% Created        R O Zhurakivsky 2012-10-15

[az,el] = view;
if evnt.Key == 'j'
	view(az-1,el);
elseif evnt.Key == 'k'
	view(az+1,el);
elseif evnt.Key == 'i'
	view(az,el-1);
elseif evnt.Key == 'm'
	view(az,el+1);
end


