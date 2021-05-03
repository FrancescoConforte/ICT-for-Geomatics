clear all
close all
clc

pos_sat = load('pos_sat.txt');

%% Calculate the coordinates Xi, Yi, Zi of the receiver;
%latitude and longitude
latitude = [45 3 48]; %phi
longitude = [7 39 41]; %lambda

%convert minutes and seconds in degree
phi = dms2degrees(latitude); %latitude in degrees
lambda = dms2degrees(longitude); %longitude in degrees

h = 0; %ellipsoidal height in meters

wgs84 = wgs84Ellipsoid('meter');
[x,y,z] = geodetic2ecef(wgs84,phi,lambda,h); %mapping toolbox


%% b) Calculate ρ and the D matrix;

rho = zeros(size(pos_sat,1),1);

for i = 1:size(pos_sat,1)
    rho(i) = sqrt((pos_sat(i,2)-x)^2 + (pos_sat(i,3)-y)^2 + (pos_sat(i,4)-z)^2);
end

D = [];

for i = 1:size(pos_sat,1)
    D(i,1) = (pos_sat(i,2)-x)/rho(i);
    D(i,2) = (pos_sat(i,3)-y)/rho(i);
    D(i,3) = (pos_sat(i,4)-z)/rho(i);
    D(i,4) = -1;
end

%% c) Calculate R and Qxx matrix;

%convert degrees in radiants
phi = deg2rad(phi); %latitude in radiant (phi)
lambda = deg2rad(lambda); %longitude in radiant (lambda)

%compute R
R = [-sin(lambda), cos(lambda), 0; -sin(phi)*cos(lambda), -sin(phi)*sin(lambda), cos(phi); cos(phi)*cos(lambda), cos(phi)*sin(lambda), sin(phi); ];
   
%compute Qxx
Qxx = inv((transpose(D)*D));

%% d) Calculate Quu matrix and estimate the GDOP, PDOP and HDOP

%calculate Quu matrix
Qxx_star = Qxx(1:3,1:3); 
Quu_star = R*Qxx_star*transpose(R);
Quu = Quu_star(1:3,1:3);
Quu = [Quu Qxx(1:3,end)];
Quu = [Quu;Qxx(end,:)];

%compute GDOP, PDOP, HDOP
GDOP = sqrt(trace(Quu));
PDOP = sqrt(trace(Quu)-Quu(4,4));
HDOP = sqrt(Quu(1,1)+Quu(2,2));

fprintf('The GDOP is %f \n',GDOP);
fprintf('The PDOP is %f \n',PDOP);
fprintf('The HDOP is %f \n',HDOP);