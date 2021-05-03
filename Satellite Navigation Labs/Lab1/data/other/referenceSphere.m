classdef referenceSphere < spheroid
%referenceSphere Reference sphere
%
%   A referenceSphere object represents a sphere with a specific name and
%   radius that can be used in map projections and other geodetic
%   operations. It has the following properties.
%
%       * A name (Name)
%
%       * A unit of length (LengthUnit) indicating the unit for the radius
%
%       * A scalar radius (Radius) 
%
%   S = referenceSphere returns a reference sphere object representing a
%   unit sphere.
%
%   S = referenceSphere(NAME) returns a reference sphere object
%   corresponding to NAME, which specifies an approximately spherical body.
%   NAME is one of the following: 'unit sphere', 'earth', 'sun', 'moon',
%   'mercury', 'venus', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune',
%   'pluto'.  NAME is case-insensitive.  The radius of the reference sphere
%   is given in meters.
%
%   S = referenceSphere(NAME, LENGTHUNIT) returns a reference sphere with
%   radius given in the specified length unit. LENGTHUNIT can be any length
%   unit supported by validateLengthUnit.
%
%   referenceSphere properties:
%      Name - Name of reference sphere
%      LengthUnit - Unit of length for radius
%      Radius - Radius of sphere
%
%   referenceSphere properties (read-only):
%      SemimajorAxis - Equatorial radius of sphere, a = Radius
%      SemiminorAxis - Distance from center of sphere to pole, b = Radius
%      InverseFlattening - Reciprocal of flattening, 1/f = Inf
%      Eccentricity - First eccentricity of sphere, ecc = 0
%      Flattening - Flattening of sphere, f = 0
%      ThirdFlattening - Third flattening of sphere, n = 0
%      MeanRadius - Mean radius of sphere
%      SurfaceArea - Surface area of sphere
%      Volume - Volume of sphere
%
%   The first 7 read-only properties are provided to ensure that reference
%   sphere objects can be used interchangeably with reference ellipsoid
%   objects in most contexts. These properties are omitted when an object
%   is displayed on the command line. The last 2 properties, SurfaceArea
%   and Volume, are also omitted, because their values are only needed in
%   certain cases.
%
%   Example
%   -------
%   % Construct a reference sphere that models the Earth as a sphere with
%   % a radius of 6,371,000 meters, then switch to use kilometers instead.
%   s = referenceSphere('Earth')
%   s.LengthUnit = 'kilometer'
%
%   % Surface area of the sphere in square kilometers
%   s.SurfaceArea
%
%   % Volume of the sphere in cubic kilometers
%   s.Volume
%
%   See also referenceEllipsoid, validateLengthUnit

% Copyright 2011-2017 The MathWorks, Inc.
    
    properties
        %Name Name of reference sphere
        %
        %   A value naming or describing the reference sphere.
        %
        %   Default value: 'Unit Sphere'
        Name = 'Unit Sphere';
    end
    
    properties (Dependent = true)
        %LengthUnit Unit of length for radius
        %
        %   An empty character vector, or any unit-of-length accepted by
        %   the validateLengthUnit function.
        %
        %   Default value: ''
        LengthUnit
        
        %Radius Radius of reference sphere
        %
        %   Positive, finite scalar.
        %
        %   Default value: 1
        Radius
    end
    
    properties (Dependent = true)
        %SemimajorAxis Equatorial radius of sphere
        %
        %   Positive, finite scalar. Its value is equal to Radius
        %   and cannot be set.
        SemimajorAxis
        
        %SemiminorAxis Distance from center of sphere to pole
        %
        %   Positive, finite scalar. Its value is equal to Radius
        %   and cannot be set.
        SemiminorAxis
    end
    
    properties (Constant = true)
        %InverseFlattening Reciprocal of flattening
        %
        %   This property is provided for consistency with the oblate
        %   spheroid class. Its value is always Inf.
        InverseFlattening = Inf;
        
        %Eccentricity First eccentricity of sphere
        %
        %   This property is provided for consistency with the oblate
        %   spheroid class. Its value is always 0.
        Eccentricity = 0;
        
        %Flattening Flattening of sphere
        %
        %   This property is provided for consistency with the oblate
        %   spheroid class. Its value is always 0.
        Flattening = 0;
        
        
        %ThirdFlattening Third flattening of sphere
        %
        %   This property is provided for consistency with the oblate
        %   spheroid class. Its value is always 0.
        ThirdFlattening = 0;
    end
    
    properties (Dependent = true)
        %MeanRadius Mean radius of sphere
        %
        %   This property is provided for consistency with the oblate
        %   spheroid class. Its value is always equal to Radius.
        MeanRadius
        
        %SurfaceArea Surface area of sphere
        %
        %   Surface area of the sphere in units of area consistent with the
        %   LengthUnit property value. For example, if LengthUnit is
        %   'kilometer' then SurfaceArea is in square kilometers.
        SurfaceArea
        
        %Volume Volume of sphere
        %
        %   Volume of the sphere in units of volume consistent with the
        %   LengthUnit property value. For example, if LengthUnit is
        %   'kilometer' then Volume is in cubic kilometers.
        Volume
    end
    
    properties (Access = private, Hidden = true)
        a = 1;             % Stores value of Radius
        pLengthUnit = '';  % Stores value of LengthUnit
    end
    
    properties (Constant, Hidden, Access = protected)
        % Control the display of referenceSphere objects.

        DerivedProperties = {'SemimajorAxis','SemiminorAxis', ...
            'InverseFlattening','Eccentricity','Flattening', ...
            'ThirdFlattening','MeanRadius','SurfaceArea','Volume'};
        
        DisplayFormat = '';   % Always use the current setting.
    end
    
    %--------------------------- Constructor ------------------------------

    methods
        
        function self = referenceSphere(name, lengthUnit)
            
            if nargin > 0
                
                % Ensure that NAME is a string or character vector.
                validateattributes(name,{'char','string'}, ...
                    {'nonempty','scalartext'},'referenceSphere','NAME',1)
                
                % Make 'unitsphere' equivalent to 'Unit Sphere'.
                if strncmpi(name,'unitsphere',length(name))
                    name = 'Unit Sphere';
                end
                
                % Names and radii in meters
                spheres = { ...
                    'Unit Sphere',    1; ...
                    'Earth',   6371000; ...
                    'Sun',     6.9446e+08; ...
                    'Moon',    1738000; ...
                    'Mercury', 2439000; ...
                    'Venus',   6051000; ...
                    'Mars',    3.39e+06; ...
                    'Jupiter', 6.9882e+07; ...
                    'Saturn',  5.8235e+07; ...
                    'Uranus',  2.5362e+07; ...
                    'Neptune', 2.4622e+07; ...
                    'Pluto',   1151000 ...
                    };
                
                validNames = spheres(:,1)';
                
                name = validatestring(name, validNames, ...
                    'referenceSphere', 'NAME', 1);
 
                self.Name = name;
                self.Radius = spheres{strcmp(name, validNames), 2};
                
                % Set LengthUnit to 'meter', except for the unit sphere.
                if ~strcmp(self.Name, 'Unit Sphere')
                    self.pLengthUnit = 'meter';
                end
                
                if nargin > 1
                    self.LengthUnit = validateLengthUnit(lengthUnit);
                end
                
            end
            
        end
        
    end
    
    %--------------------------- Get methods ------------------------------
    
    methods
        
        function lengthUnit = get.LengthUnit(self)
            lengthUnit = self.pLengthUnit;
        end
        
        function radius = get.Radius(self)
            radius = self.a;
        end
        
        function a = get.SemimajorAxis(self)
            a = self.a;
        end
        
        function b = get.SemiminorAxis(self)
            b = self.a;
        end
        
        function radius = get.MeanRadius(self)
            radius = self.a;
        end
        
        function surfarea = get.SurfaceArea(self)
            surfarea = 4 * pi * self.a^2;
        end
        
        function vol = get.Volume(self)
            vol = (4*pi/3) * self.a^3;
        end
    end
    
    %---------------------------- Set methods -----------------------------
    
    methods
        
        function self = set.Name(self, name)
            % Accept any valid name, including empty.
            if ~isempty(name)
                validateattributes(name, {'char','string'}, {'scalartext'}, '', 'Name')
            else
                validateattributes(name, {'char','string'}, {}, '', 'Name')
            end
            self.Name = name;
        end
        
        
        function self = set.LengthUnit(self, unit)
            % If the length unit is not yet designated, validate the input
            % and assign it. Otherwise rescale the radius, then make
            % the assignment.
            if isempty(unit)
                % It would be unusual to change LengthUnit to empty, but
                % there's no reason to disallow it. Ensure an empty char,
                % vs. empty of some class other than char.
                self.pLengthUnit = '';
            else
                unit = validateLengthUnit(unit);
                if ~isempty(self.pLengthUnit)
                    ratio = unitsratio(unit, self.pLengthUnit);
                    self.a = ratio * self.a;
                end
                self.pLengthUnit = unit;
            end
        end
        
        
        function self = set.Radius(self, radius)
            validateattributes(radius, ...
                {'double'}, {'real','positive','finite','scalar'}, ...
                '', 'Radius');
            self.a = radius;
        end
        
    end
    
    %------------------- 3-D Coordinate Transformations -------------------
    
    methods (Hidden)
        
        function [x, y, z] = geodetic2ecef(spheroid, phi, lambda, h, angleUnit)
            
            if ~isobject(spheroid) && isobject(h)
                % Older syntax with PHI and LAMBDA in radians:
                %    [X, Y, Z] = geodetic2ecef(PHI, LAMBDA, H, SPHEROID)
                % Also do the right thing with this undocumented syntax:
                %    [X, Y, Z] = geodetic2ecef(PHI, LAMBDA, H, SPHEROID, angleUnit)
                [phi, lambda, h, spheroid] = deal(spheroid, phi, lambda, h);
                inDegrees = ~(nargin < 5 || ~map.geodesy.isDegree(angleUnit));
            else
                inDegrees = nargin < 5 || map.geodesy.isDegree(angleUnit);
            end
            
            if ~(isscalar(spheroid) && isobject(spheroid))
                validateattributes(spheroid, ...
                    {'oblateSpheroid','referenceEllipsoid','referenceSphere'}, ...
                    {'scalar'}, '', 'spheroid')
            end
        
            r = spheroid.a + h;
            if inDegrees
                rho = r .* cosd(phi);
                x = rho .* cosd(lambda);
                y = rho .* sind(lambda);
                z = r .* sind(phi);
            else
                [x, y, z] = sph2cart(lambda, phi, r);
            end
            
        end
        
        
        function [phi, lambda, h] = ecef2geodetic(spheroid, x, y, z, angleUnit)
            
            if ~isobject(spheroid) && isobject(z)
                % Older syntax with PHI and LAMBDA in radians:
                %    [PHI, LAMBDA, H] = ecef2geodetic(X, Y, Z, SPHEROID)
                % Also do the right thing with this undocumented syntax:
                %    [PHI, LAMBDA, H] = ecef2geodetic(X, Y, Z, SPHEROID, angleUnit)
                
                [x, y, z, spheroid] = deal(spheroid, x, y, z);
                inDegrees = ~(nargin < 5 || ~map.geodesy.isDegree(angleUnit));
            else
                inDegrees = nargin < 5 || map.geodesy.isDegree(angleUnit);
            end
            
            if ~(isscalar(spheroid) && isobject(spheroid)) 
                validateattributes(spheroid, ...
                    {'oblateSpheroid','referenceEllipsoid','referenceSphere'}, ...
                    {'scalar'}, '', 'spheroid')
            end
            
            if inDegrees
                rho = hypot(x,y);
                phi = atan2d(z,rho);
                lambda = atan2d(y,x);
                r = hypot(z,rho);
            else
                [lambda, phi, r] = cart2sph(x, y, z);
            end
            h = r - spheroid.a;           
        end
        
        
        function [deltaX, deltaY, deltaZ] = ecefOffset(spheroid, ...
                phi1, lambda1, h1, phi2, lambda2, h2, angleUnit)
            
            if ~(isscalar(spheroid) && isobject(spheroid))
                validateattributes(spheroid, ...
                    {'oblateSpheroid','referenceEllipsoid','referenceSphere'}, ...
                    {'scalar'}, '', 'spheroid')
            end
            
            % See implementation notes for oblateSpheroid/ecefOffset. The
            % following uses the same strategy, but the expressions are
            % simpler because the eccentricity is 0.

            if nargin < 8 || map.geodesy.isDegree(angleUnit)
                sinfun = @sind;
                cosfun = @cosd;
            else
                sinfun = @sin;
                cosfun = @cos;
            end
            
            s1 = sinfun(phi1);
            c1 = cosfun(phi1);
            
            s2 = sinfun(phi2);
            c2 = cosfun(phi2);
            
            p1 = c1 .* cosfun(lambda1);
            p2 = c2 .* cosfun(lambda2);
            
            q1 = c1 .* sinfun(lambda1);
            q2 = c2 .* sinfun(lambda2);
            
            deltaX = spheroid.a * (p2 - p1) + (h2 .* p2 - h1 .* p1);
            deltaY = spheroid.a * (q2 - q1) + (h2 .* q2 - h1 .* q1);
            deltaZ = spheroid.a * (s2 - s1) + (h2 .* s2 - h1 .* s1);
        end        
    end
end
