function [f,R] = worldparams(narg, varg)
% WORLDPARAMS internal function returning flattening and equatorial radius
% of planets. 
%
%   See also ECEF2LLA, GEOC2GEOD, GEOD2GEOC, GEOCRADIUS, LLA2ECEF, WGS84MODEL. 

%   Copyright 2000-2016 The MathWorks, Inc.

if (narg == 1)
    % use WGS84 model
    [f,R] = wgs84model;
end

if (narg == 2) 
    if ~ischar(varg{1}) && ~isstring(varg{1})
        error(message('aero:worldparams:notChar'));
    end
    % get world model type - only WGS84 supported for now
    type = char(varg{1});
    switch lower( type )
        case 'wgs84'
            [f,R] = wgs84model;
        otherwise
            error(message('aero:worldparams:invalidModel', type));
    end
end

if (narg == 3) 
    if (~isnumeric(varg{1}) || ~isscalar(varg{1}))
        error(message('aero:worldparams:notNumeric2'));
    end
    if (~isnumeric(varg{2}) || ~isscalar(varg{2}))
        error(message('aero:worldparams:notNumeric3'));
    end
    f = varg{1};
    R = varg{2};
end
