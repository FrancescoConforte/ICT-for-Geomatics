%% Function to find longitude and height from ECEF coordinates
function [phi, h] = xyz2llh(ellipsoid, points, lambda, threshold, n)
    if strcmp(n, 'WGS84')
        h0=0;
        r = sqrt(points(1)^2+points(2)^2);       
        phi0 = atan(points(3)/r);
        while 1
            e2 = ellipsoid{9,5};
            N = ellipsoid{9,2}/sqrt((1- e2 *sin(phi0)^2));
            h = points(1)/(cos(phi0)*cos(lambda))-N;
            phi = atan(points(3)/(r*(1-((e2*N)/(N+h)))));
            if ((phi-phi0) < threshold) && ((h-h0) < threshold)
                break;
            else
                phi0 = phi;
                h0 = h;
            end
        end
    elseif strcmp(n, 'Hayford')
        h0=0;
        r = sqrt(points(1)^2+points(2)^2);       
        phi0 = atan(points(3)/r);
        while 1
            e2 = ellipsoid{7,5};
            N = ellipsoid{7,2}/sqrt((1- e2 *sin(phi0)^2));
            h = points(1)/(cos(phi0)*cos(lambda))-N;
            phi = atan(points(3)/(r*(1-((e2*N)/(N+h)))));
            if ((phi-phi0) < threshold) && ((h-h0) < threshold)
                break;
            else
                phi0 = phi;
                h0 = h;
            end
        end
    end
end