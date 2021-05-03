clear all
close all
clc

%% a)Convert latitude and longitude of the station from degrees to radians

%latitude and longitude
latitude = [45 3 48.114];
longitude = [7 39 40.605];

%convert minutes and seconds in degree
LatitudeInDegrees = dms2degrees(latitude);
LongitudeInDegrees = dms2degrees(longitude);

%convert degrees in radiants
phi = deg2rad(LatitudeInDegrees); %latitude in radiant (phi)
lambda = deg2rad(LongitudeInDegrees); %longitude in radiant (lambda)

%% b)Calculate e^2 and W

%flattening
f = 1/298.257223;

%For the WGS84 ellipsoid to model Earth, the defining values are
% a (equatorial radius): 6 378 137.0 m
% 1/f (inverse flattening): 298.257 223 563
% from which one derives
% b (polar radius): 6 356 752.3142 m

%eccentricity 
e2 = 2*f-f^2;

%W 
W = sqrt((1-e2*sin(phi)^2));

%% c) Calculate the station position (in ECEF coordinates) Xi, Yi, Zi 
%semi-major axis of the Earth [m]
a = 6378137;

%coordinates
Xi = (a*cos(phi)*cos(lambda))/W;
Yi = (a*cos(phi)*sin(lambda))/W;
Zi = (a*(1-e2)*sin(phi))/W;

%% d) Calculate local coordinates e, n, u

%Satellite coordinates (ECEF)
Xsat = 15487292.829;
Ysat = 6543538.932;
Zsat = 20727274.429;

%Convert ECEF coordinates to local ENU coordinates
matrix = [-sin(lambda), cos(lambda), 0; -sin(phi)*cos(lambda), -sin(phi)*sin(lambda), cos(phi); cos(phi)*cos(lambda), cos(phi)*sin(lambda), sin(phi); ];
LocalCoordinates = matrix * [Xsat-Xi; Ysat-Yi; Zsat-Zi];
e = LocalCoordinates(1);
n = LocalCoordinates(2);
u = LocalCoordinates(3);

%% Estimate Azimuth and Elevation

azimuth = atan(u/sqrt(n^2+e^2));
az=rad2deg(azimuth);

elevation = atan(e/n);
elev = rad2deg(elevation);

fprintf('The azimuth is %f degrees\n',az);
fprintf('The elevation is %f degrees\n',elev);


