%%
clear all
close all
clc

%% Exercise 1: Ellipsoid's parameters
%Name Ellipsoid, semi-mahor axis a, flattening (alpha) 
ellipsoid = [{'Delambre', 6376985,1/308.6};
{'Everest', 6377276,1/300.8};
{'Bessel', 6377397,1/299.1528128}; 
{'Fisher', 6378338,1/288.5};
{'Clarke', 6378249,1/293.5};
{'Helmert', 6378140,1/298.3};
{'Hayford', 6378388,1/297};
{'Krassovsky', 6378245,1/298.3}; 
{'WGS84', 6378137,1/298.257223563}];


for i= 1:size(ellipsoid,1)
    ellipsoid{i,4} = ellipsoid{i,2}*(1 - ellipsoid{i,3}); %c = a(1-alpha) semi-minor axis
    ellipsoid{i,5} = (ellipsoid{i,2}^2 - ellipsoid{i,4}^2)/ellipsoid{i,2}^2; %first eccentricity
    ellipsoid{i,6} = (ellipsoid{i,2}^2 - ellipsoid{i,4}^2)/ellipsoid{i,4}^2; %second eccentricity
end

% writecell(ellipsoid,'Ellipsoid''s parameters.dat') %obtain a .txt file with all parameters

%% Exercise 2: Coordinates transformation from Geographic to Geocentric

%matrix containing the two points: latitude and longitude are trasformed
%from degree-minutes-seconds to degrees with function dms2degrees()
coordinates_geo = [dms2degrees([44 45 01.03930]),dms2degrees([7 39 40.605]),322.4909;
    dms2degrees([44 47 10.90505]),dms2degrees([7 30 26.53939]),305.7367];


%Using the functions x(), y(), z() defined in other files, the cartesian
%coordinates for both the points and by considering both the ellipsoids:
%WGS84 and Hayford. The parameters are taken from the previous exercise
wgs84_coordinates_cart = [];
hay_coordinates_cart = [];
for i = 1:size(coordinates_geo,1)
    wgs84_coordinates_cart(i,1) = x(ellipsoid, coordinates_geo(i,1),coordinates_geo(i,2), coordinates_geo(i,3),'WGS84');
    wgs84_coordinates_cart(i,2) = y(ellipsoid, coordinates_geo(i,1),coordinates_geo(i,2), coordinates_geo(i,3),'WGS84');
    wgs84_coordinates_cart(i,3) = z(ellipsoid, coordinates_geo(i,1),coordinates_geo(i,2), coordinates_geo(i,3),'WGS84');
    hay_coordinates_cart(i,1) = x(ellipsoid, coordinates_geo(i,1),coordinates_geo(i,2), coordinates_geo(i,3),'Hayford');
    hay_coordinates_cart(i,2) = y(ellipsoid, coordinates_geo(i,1),coordinates_geo(i,2), coordinates_geo(i,3),'Hayford');
    hay_coordinates_cart(i,3) = z(ellipsoid, coordinates_geo(i,1),coordinates_geo(i,2), coordinates_geo(i,3),'Hayford');    
end

% writematrix(wgs84_coordinates_cart,'ECEF Coordinates WGS84(Ex2).txt');
% writematrix(hay_coordinates_cart,'ECEF Coordinates Hayford(Ex2).txt');


%Try to change «h1» , adding 2000 m.
h2000coordinates_cart = [];

h2000coordinates_cart(1,1) = x(ellipsoid, coordinates_geo(1,1),coordinates_geo(1,2), coordinates_geo(1,3)+2000,'WGS84');
h2000coordinates_cart(1,2) = y(ellipsoid, coordinates_geo(1,1),coordinates_geo(1,2), coordinates_geo(1,3)+2000,'WGS84');
h2000coordinates_cart(1,3) = z(ellipsoid, coordinates_geo(1,1),coordinates_geo(1,2), coordinates_geo(1,3)+2000,'WGS84');



%% Exercise 3: Coordinates transformation from Geocentric to Geographic

points = [4499525.4271, 585034.1293, 4467910.3596;
          4495694.2695, 592457.8605, 4470744.7781;
          4503484.7172, 578160.7507, 4465024.3002;
          4498329.3715, 562840.7651, 4472537.6125];
      
epsilon = 1e-8;

wgs84coord_geograph = [];
haycoord_geograph = [];
for i = 1: size(points,1)
    wgs84coord_geograph(i,2) = atan(points(i,2)/points(i,1)); %longitude (lambda)
    [wgs84coord_geograph(i,1), wgs84coord_geograph(i,3)] = xyz2llh(ellipsoid, points(i,:), wgs84coord_geograph(i,2), epsilon, 'WGS84');
    haycoord_geograph(i,2) = atan(points(i,2)/points(i,1)); %longitude (lambda)
    [haycoord_geograph(i,1), haycoord_geograph(i,3)] = xyz2llh(ellipsoid, points(i,:), haycoord_geograph(i,2), epsilon, 'Hayford');
end

%From rad to deg
for i = 1: size(points,1)
    wgs84coord_geograph(i,1) = rad2deg(wgs84coord_geograph(i,1));
    wgs84coord_geograph(i,2) = rad2deg(wgs84coord_geograph(i,2));
    haycoord_geograph(i,1) = rad2deg(haycoord_geograph(i,1));
    haycoord_geograph(i,2) = rad2deg(haycoord_geograph(i,2));
end

% Save two files for the report
% writematrix(wgs84coord_geograph,'Geographic CoordinatesWGS84(Ex3).txt');
% writematrix(haycoord_geograph,'Geographic CoordinatesHayford(Ex3).txt');
% writematrix(points, 'points(ex3).txt');

%% Exercise 5: Helmert parameters
pointsHelmert = load('points_helmert.TXT');

%In order to avoid singular matrix, it is better to divide the coordinates (X, Y, Z) by 10^6.
pointsHelDiv = pointsHelmert/1e+6;

A = zeros(36,7);
l0 = zeros(36,1);

for i = 1:12
    A(1+3*i,1) = 1;
    A(2+3*i,2) = 1;
    A(3+3*i,3) = 1;
    A(1+3*i,4) = pointsHelDiv(i,1);
    A(2+3*i,4) = pointsHelDiv(i,2);
    A(3+3*i,4) = pointsHelDiv(i,3);
    A(2+3*i,5) = -pointsHelDiv(i,3);
    A(3+3*i,5) = pointsHelDiv(i,2);
    A(1+3*i,6) = -pointsHelDiv(i,3);
    A(3+3*i,6) = pointsHelDiv(i,1);
    A(1+3*i,7) = pointsHelDiv(i,2);
    A(2+3*i,7) = -pointsHelDiv(i,1);
    l0(1+3*i,1) = pointsHelDiv(i,4) - pointsHelDiv(i,1);
    l0(2+3*i,1) = pointsHelDiv(i,5) - pointsHelDiv(i,2);
    l0(3+3*i,1) = pointsHelDiv(i,6) - pointsHelDiv(i,3);
end

x = inv(A' * A) * (A' * l0);
res = A * x - l0;
x = x*1e+6;
res = res*1e+6;

r = max(res);





















