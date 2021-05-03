clear all
close all
clc

%% Exercise 1: Modulus of linear deformation

phi = deg2rad([30;37;45;60]); %latitude
m_l = [];
lambda = [];
for i = 1: size(phi,1)
    count = 1;
    for j = 0:0.5:3
        m_l(i,count) = 0.9996*(1 + (deg2rad(j)^2/2)*cos(phi(i))^2);
        lambda(count) = j; %longitude in deg
        count = count + 1;
    end 
    plot(lambda(1,:),m_l(i,:),'-o')
    hold on
end
grid on
ylabel('Modulus of linear deformation')
xlabel('Longitude [degree]')
legend('Latitude 30°','Latitude 37°','Latitude 45°','Latitude 60°','Location','best')
hold off

%% Exercise 2: From Geographic to Cartographic coordiantes with Hirvonen transformation

a = 6378137;
alpha = 1 / 298.257223563;
c = a * (1 - alpha);
e2 = (a^2 - c^2) / a^2;
second_eccentricity = (a^2 - c^2) / c^2;

%Convert from dms to degrees
lat1 = dms2degrees([45 3 45.717]);
lon1 = dms2degrees([7 47 26.292]);
lat2 = dms2degrees([38 32 34.649]);
lon2 = dms2degrees([16 50 06.493]);

%Points converted in radians
points = [deg2rad(lat1), deg2rad(lon1); deg2rad(lat2), deg2rad(lon2)]; % latitude and longitude of the two points 

centr_mer = [deg2rad(9); deg2rad(15)]; %central meridian for the two points in radians
mc = 0.9996;
false_east = 500000;
A1 = 1 - (e2/4) - (3*e2^2/64) - (5*e2^3/256);
A2 = (3*e2/8) + (3*e2^2/32) + (45*e2^3/1024);
A4 = (15*e2^2/256) + (45*e2^3/1024);
A6 = (35*e2^3/3072);
R_p = a^2/c;


x = zeros(2,1);
y = zeros(2,1);
East = zeros(2,1);
North = zeros(2,1);
for i = 1:2
    lambda_prime = points(i,2) - centr_mer(i);
    v1 = sqrt(1+(second_eccentricity*cos(points(i,1))^2));
    csi = atan(tan(points(i,1))/cos(v1*lambda_prime));
    v = sqrt(1+(second_eccentricity*(cos(csi))^2));
    x(i) = R_p*asinh((cos(csi)*tan(lambda_prime))/v);
    y(i) = a*((A1*csi) - (A2*sin(2*csi)) + (A4*sin(4*csi)) - (A6*sin(6*csi)));
    East(i) = (x(i)*mc) + false_east;
    North(i) = y(i)*mc;
end

%% Exercise 3: from cartographic (E,N) to geographic (phi,lambda)

E = 470139.66;
N = 5031468.37;

x2 = (E - false_east)/mc;
y2 = N/mc;

theta = y2/(a*A1);
E1 = (a-c)/(a+c);
B2 = (3*E1)/2 - (27*E1^3)/32;
B4 = (21*E1^2)/16 - (55*E1^4)/32;
B6 = (151*E1^3)/96;
B8 = (1097*E1^4)/512;
csi2 = theta + (B2*sin(2*theta)) + (B4*sin(4*theta)) + (B6*sin(6*theta)) + (B8*sin(8*theta));
v = sqrt(1+(second_eccentricity*(cos(csi2))^2));
lambda_prime2 = atan((v*sinh(x2/R_p))/cos(csi2)); 

phi = atan(tan(csi2)*cos(v*lambda_prime2));
lambda = lambda_prime2 + centr_mer(1);

phi = rad2deg(phi);
lambda = rad2deg(lambda);

latitude = degrees2dms(phi);
longitude = degrees2dms(lambda);


