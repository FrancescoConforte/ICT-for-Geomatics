function value_z = z(ellipsoid, phi, lambda, h, n)
    if strcmp(n, 'WGS84') %WGS84
       a = ellipsoid{9,2}; %semi-major axis
       e2 = ellipsoid{9,5}; %first eccentricity
       W = sqrt((1- e2 *sin(deg2rad(phi))^2));
    elseif strcmp(n, 'Hayford')
        a = ellipsoid{7,2}; %semi-major axis
        e2 = ellipsoid{7,5}; %first eccentricity
        W = sqrt((1- e2 *sin(deg2rad(phi))^2));
    end
    value_z = ((a/W)*(1 - e2) + h)*sin(deg2rad(phi));
end