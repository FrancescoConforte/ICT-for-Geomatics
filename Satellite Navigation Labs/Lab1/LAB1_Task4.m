clear all
close all
clc

%% TASK A-4: Estimate the sigma_j,UERE for all satellites using the dataset from the realisticUERE folder.
%% Implement the Weigthed Least Mean Square (WLMS) repeating Tasks 2 and 3.

%load('dataset_1_20180329T160947.mat')
%load('dataset_2_20180329T160900.mat')
%load('dataset_3_20180329T161023.mat')
%load('dataset_4_20180329T161103.mat')
load('dataset_5_20180329T161418.mat')
%load('dataset_6_20180329T161139.mat')

%GPS
N = size(RHO.GPS,2); %3600
J_gps = size(RHO.GPS,1); %32

R_gps = [];
i = 1;
for j = 1:J_gps
    if isnan(RHO.GPS(j,:))
        continue;
    else
        %remove the variation along the time
        rho_removed(i,:) = diff(RHO.GPS(j,:),2); 
%         plot(rho','.');
%         hold on
        R_gps(i,i) = var(rho_removed(i,:));
        %Estimate the sigma_j_UERE for all visible satellites
        sigma_j_UERE_gps(i,1) = j;
        sigma_j_UERE_gps(i,2) = std(rho_removed(i,:));
        i = i + 1;
    end
end
% plot(rho_removed(10,:),'-o');
% hold on
% plot(rho_removed(5,:),'-o');
% grid on
% xlabel('time [s]');
% ylabel('meters');
% title('d^{2}\rho/dt^{2}');
% legend('PRN10','PRN5');
sigma_UERE_GPS = sqrt(trace(R_gps)); %standard deviation of the position error
W_gps = inv(R_gps);

%GLONASS
N = size(RHO.GLO,2); %3600
J_glo = size(RHO.GLO,1); 

R_glo = [];
i = 1;
rho_removed = [];
for j = 1:J_glo
    if isnan(RHO.GLO(j,:))
        continue;
    else
        %remove the variation along the time
        rho_removed(i,:) = diff(RHO.GLO(j,:),2); 
%         plot(rho','.');
%         hold on
        R_glo(i,i) = var(rho_removed(i,:));
        %Estimate the sigma_j_UERE for all satellites
        sigma_j_UERE_glo(i,1) = j;
        sigma_j_UERE_glo(i,2) = std(rho_removed(i,:));
        i = i + 1;
    end
end
sigma_UERE_GLO = sqrt(trace(R_glo)); %standard deviation of the position error
W_glo = inv(R_glo);

%GALILEO
N = size(RHO.GAL,2); %3600
J_gal = size(RHO.GAL,1); 

R_gal = [];
i = 1;
rho_removed = [];
for j = 1:J_gal
    if isnan(RHO.GAL(j,:))
        continue;
    else
        %remove the variation along the time
        rho_removed(i,:) = diff(RHO.GAL(j,:),2); 
%         plot(rho_removed','.');
%         hold on
        R_gal(i,i) = var(rho_removed(i,:));
        %Estimate the sigma_j_UERE for all satellites
        sigma_j_UERE_gal(i,1) = j;
        sigma_j_UERE_gal(i,2) = std(rho_removed(i,:));
        i = i + 1;
    end
end
sigma_UERE_GAL = sqrt(trace(R_gal)); %standard deviation of the position error
W_gal = inv(R_gal);

%BEIDOU
N = size(RHO.BEI,2); %3600
J_bei = size(RHO.BEI,1); 

R_bei = [];
i = 1;
rho_removed = [];
for j = 1:J_bei
    if isnan(RHO.BEI(j,:))
        continue;
    else
        %remove the variation along the time
        rho_removed(i,:) = diff(RHO.BEI(j,:),2); 
%         plot(rho_removed','.');
%         hold on
        R_bei(i,i) = var(rho_removed(i,:));
        %Estimate the sigma_j_UERE for all satellites
        sigma_j_UERE_bei(i,1) = j;
        sigma_j_UERE_bei(i,2) = std(rho_removed(i,:));
        i = i + 1;
    end
end
sigma_UERE_BEI = sqrt(trace(R_bei)); %standard deviation of the position error
W_bei = inv(R_bei);

%GPS with WLMS

K = 9; 
for n = 1:N
    rho = [];
    satellitePosition = [];
    rho_hat = [];
    r_hat = [];
    H = [];
    i = 1;
    for j = 1:J_gps
        if isnan(RHO.GPS(j,n))
            continue;
        else
            rho(1,i) = RHO.GPS(j,n); %build the measurement vector rho_n = [rho1,n .... rhoJ,n] of visible satellites
            satellitePosition(i,:) = SAT_POS_ECEF.GPS(j).pos(n,:); %retrieve the corresponding satellites coordinates yj,n = [xj,n yj,n zj,n]
            i = i + 1;
        end
    end
    
    x_vector_hat = zeros(1,4); %inizialization of the estimated position/state
    for k = 1:K
        for j = 1:size(satellitePosition,1)
            %compute PVT solution 
            r_hat(j) = sqrt((satellitePosition(j,1)-x_vector_hat(1,1))^2 + (satellitePosition(j,2)-x_vector_hat(1,2))^2 + ...
                (satellitePosition(j,3)-x_vector_hat(1,3))^2); %compute the euclidean distance between the estimated user point and the satellite
            rho_hat(1,j) = r_hat(j) + x_vector_hat(1,4); %estimation of pseudoranges at the k-th iteration
            H(j,1) = (satellitePosition(j,1)-x_vector_hat(1,1))/r_hat(j);
            H(j,2) = (satellitePosition(j,2)-x_vector_hat(1,2))/r_hat(j);
            H(j,3) = (satellitePosition(j,3)-x_vector_hat(1,3))/r_hat(j);
            H(j,4) = 1;
        end
        
        deltaRho_hat = rho_hat - rho;
        H_bar = inv(transpose(H)*W_gps*H)*transpose(H)*W_gps;
        deltaX_hat = H_bar*deltaRho_hat';
        x_vector_hat = x_vector_hat + deltaX_hat'; %update of user location          
    end
    %convert from ECEF to LLA and create a .kml file for GoogleEarth
    lla(n,:) = ecef2lla(x_vector_hat(1:3));
    x_vector(n,:) = x_vector_hat;
end
position_GPS_WLMS = mean(lla);
% writeKML_GoogleEarth('Dataset2WLMS_gps',position(1),position(2),position(3));

%% TASK A-3: Compute the standard deviation of the position error over time n=1,2,..,N
sigma_x2 = var(x_vector(:,1));
sigma_y2 = var(x_vector(:,2));
sigma_z2 = var(x_vector(:,3));
sigma_but2 = var(x_vector(:,4));
sigma_x_WLMS_gps = sqrt(sigma_x2 + sigma_y2 + sigma_z2);

%GPS with LMS

K = 9; 
lla = [];

for n = 1:N
    rho = [];
    satellitePosition = [];
    rho_hat = [];
    r_hat = [];
    H = [];
    i = 1;
    for j = 1:J_gps
        if isnan(RHO.GPS(j,n))
            continue;
        else
            rho(1,i) = RHO.GPS(j,n); %build the measurement vector rho_n = [rho1,n .... rhoJ,n] of visible satellites
            satellitePosition(i,:) = SAT_POS_ECEF.GPS(j).pos(n,:); %retrieve the corresponding satellites coordinates yj,n = [xj,n yj,n zj,n]
            i = i + 1;
        end
    end
    
    x_vector_hat = zeros(1,4); %inizialization of the estimated position/state
    for k = 1:K
        for j = 1:size(satellitePosition,1)
            %compute PVT solution 
            r_hat(j) = sqrt((satellitePosition(j,1)-x_vector_hat(1,1))^2 + (satellitePosition(j,2)-x_vector_hat(1,2))^2 + ...
                (satellitePosition(j,3)-x_vector_hat(1,3))^2); %compute the euclidean distance between the estimated user point and the satellite
            rho_hat(1,j) = r_hat(j) + x_vector_hat(1,4); %estimation of pseudoranges at the k-th iteration
            H(j,1) = (satellitePosition(j,1)-x_vector_hat(1,1))/r_hat(j);
            H(j,2) = (satellitePosition(j,2)-x_vector_hat(1,2))/r_hat(j);
            H(j,3) = (satellitePosition(j,3)-x_vector_hat(1,3))/r_hat(j);
            H(j,4) = 1;
        end
        deltaRho_hat = rho_hat - rho;
        deltaX_hat = inv((transpose(H)*H))*transpose(H)*deltaRho_hat';
        x_vector_hat = x_vector_hat + deltaX_hat'; %update of user location          
    end
    %convert from ECEF to LLA and create a .kml file for GoogleEarth
    lla(n,:) = ecef2lla(x_vector_hat(1:3));
    x_vector(n,:) = x_vector_hat;
end
position_GPS_LMS = mean(lla);
% writeKML_GoogleEarth('Dataset2LMS_gps',position(1),position(2),position(3));

%% TASK A-3: Compute the standard deviation of the position error over time n=1,2,..,N
sigma_x2 = var(x_vector(:,1));
sigma_y2 = var(x_vector(:,2));
sigma_z2 = var(x_vector(:,3));
sigma_but2 = var(x_vector(:,4));
sigma_x_LMS_gps = sqrt(sigma_x2 + sigma_y2 + sigma_z2);


%GLONASS with WLMS

K = 9; 
for n = 1:N
    rho = [];
    satellitePosition = [];
    rho_hat = [];
    r_hat = [];
    H = [];
    i = 1;
    for j = 1:J_glo
        if isnan(RHO.GLO(j,n))
            continue;
        else
            rho(1,i) = RHO.GLO(j,n); %build the measurement vector rho_n = [rho1,n .... rhoJ,n] of visible satellites
            satellitePosition(i,:) = SAT_POS_ECEF.GLO(j).pos(n,:); %retrieve the corresponding satellites coordinates yj,n = [xj,n yj,n zj,n]
            i = i + 1;
        end
    end
    
    x_vector_hat = zeros(1,4); %inizialization of the estimated position/state
    for k = 1:K
        for j = 1:size(satellitePosition,1)
            %compute PVT solution 
            r_hat(j) = sqrt((satellitePosition(j,1)-x_vector_hat(1,1))^2 + (satellitePosition(j,2)-x_vector_hat(1,2))^2 + ...
                (satellitePosition(j,3)-x_vector_hat(1,3))^2); %compute the euclidean distance between the estimated user point and the satellite
            rho_hat(1,j) = r_hat(j) + x_vector_hat(1,4); %estimation of pseudoranges at the k-th iteration
            H(j,1) = (satellitePosition(j,1)-x_vector_hat(1,1))/r_hat(j);
            H(j,2) = (satellitePosition(j,2)-x_vector_hat(1,2))/r_hat(j);
            H(j,3) = (satellitePosition(j,3)-x_vector_hat(1,3))/r_hat(j);
            H(j,4) = 1;
        end
        
        deltaRho_hat = rho_hat - rho;
        H_bar = inv(transpose(H)*(W_glo*H))*transpose(H)*W_glo;
        deltaX_hat = H_bar*deltaRho_hat';
        x_vector_hat = x_vector_hat + deltaX_hat'; %update of user location          
    end
    %convert from ECEF to LLA and create a .kml file for GoogleEarth
    lla(n,:) = ecef2lla(x_vector_hat(1:3));
    x_vector(n,:) = x_vector_hat;
end
position_GLONASS_WLMS = mean(lla);
% writeKML_GoogleEarth('Dataset2WLMS_glonass',position(1),position(2),position(3));

%% TASK A-3: Compute the standard deviation of the position error over time n=1,2,..,N
sigma_x2 = var(x_vector(:,1));
sigma_y2 = var(x_vector(:,2));
sigma_z2 = var(x_vector(:,3));
sigma_but2 = var(x_vector(:,4));
sigma_x_WLMS_glonass = sqrt(sigma_x2 + sigma_y2 + sigma_z2);

%GLONASS with LMS

K = 9; 
lla = [];

for n = 1:N
    rho = [];
    satellitePosition = [];
    rho_hat = [];
    r_hat = [];
    H = [];
    i = 1;
    for j = 1:J_glo
        if isnan(RHO.GLO(j,n))
            continue;
        else
            rho(1,i) = RHO.GLO(j,n); %build the measurement vector rho_n = [rho1,n .... rhoJ,n] of visible satellites
            satellitePosition(i,:) = SAT_POS_ECEF.GLO(j).pos(n,:); %retrieve the corresponding satellites coordinates yj,n = [xj,n yj,n zj,n]
            i = i + 1;
        end
    end
    
    x_vector_hat = zeros(1,4); %inizialization of the estimated position/state
    for k = 1:K
        for j = 1:size(satellitePosition,1)
            %compute PVT solution 
            r_hat(j) = sqrt((satellitePosition(j,1)-x_vector_hat(1,1))^2 + (satellitePosition(j,2)-x_vector_hat(1,2))^2 + ...
                (satellitePosition(j,3)-x_vector_hat(1,3))^2); %compute the euclidean distance between the estimated user point and the satellite
            rho_hat(1,j) = r_hat(j) + x_vector_hat(1,4); %estimation of pseudoranges at the k-th iteration
            H(j,1) = (satellitePosition(j,1)-x_vector_hat(1,1))/r_hat(j);
            H(j,2) = (satellitePosition(j,2)-x_vector_hat(1,2))/r_hat(j);
            H(j,3) = (satellitePosition(j,3)-x_vector_hat(1,3))/r_hat(j);
            H(j,4) = 1;
        end
        deltaRho_hat = rho_hat - rho;
        deltaX_hat = inv((transpose(H)*H))*transpose(H)*deltaRho_hat';
        x_vector_hat = x_vector_hat + deltaX_hat'; %update of user location          
    end
    %convert from ECEF to LLA and create a .kml file for GoogleEarth
    lla(n,:) = ecef2lla(x_vector_hat(1:3));
    x_vector(n,:) = x_vector_hat;
end
position_GLONASS_LMS = mean(lla);
%writeKML_GoogleEarth('Dataset2LMS_glonass',position(1),position(2),position(3));

%% TASK A-3: Compute the standard deviation of the position error over time n=1,2,..,N
sigma_x2 = var(x_vector(:,1));
sigma_y2 = var(x_vector(:,2));
sigma_z2 = var(x_vector(:,3));
sigma_but2 = var(x_vector(:,4));
sigma_x_LMS_glonass = sqrt(sigma_x2 + sigma_y2 + sigma_z2);

%GALILEO with WLMS

K = 9; 
for n = 1:N
    rho = [];
    satellitePosition = [];
    rho_hat = [];
    r_hat = [];
    H = [];
    i = 1;
    for j = 1:J_gal
        if isnan(RHO.GAL(j,n))
            continue;
        else
            rho(1,i) = RHO.GAL(j,n); %build the measurement vector rho_n = [rho1,n .... rhoJ,n] of visible satellites
            satellitePosition(i,:) = SAT_POS_ECEF.GAL(j).pos(n,:); %retrieve the corresponding satellites coordinates yj,n = [xj,n yj,n zj,n]
            i = i + 1;
        end
    end
    
    x_vector_hat = zeros(1,4); %inizialization of the estimated position/state
    for k = 1:K
        for j = 1:size(satellitePosition,1)
            %compute PVT solution 
            r_hat(j) = sqrt((satellitePosition(j,1)-x_vector_hat(1,1))^2 + (satellitePosition(j,2)-x_vector_hat(1,2))^2 + ...
                (satellitePosition(j,3)-x_vector_hat(1,3))^2); %compute the euclidean distance between the estimated user point and the satellite
            rho_hat(1,j) = r_hat(j) + x_vector_hat(1,4); %estimation of pseudoranges at the k-th iteration
            H(j,1) = (satellitePosition(j,1)-x_vector_hat(1,1))/r_hat(j);
            H(j,2) = (satellitePosition(j,2)-x_vector_hat(1,2))/r_hat(j);
            H(j,3) = (satellitePosition(j,3)-x_vector_hat(1,3))/r_hat(j);
            H(j,4) = 1;
        end
        
        deltaRho_hat = rho_hat - rho;
        H_bar = inv(transpose(H)*W_gal*H)*transpose(H)*W_gal;
        deltaX_hat = H_bar*deltaRho_hat';
        x_vector_hat = x_vector_hat + deltaX_hat'; %update of user location          
    end
    %convert from ECEF to LLA and create a .kml file for GoogleEarth
    lla(n,:) = ecef2lla(x_vector_hat(1:3));
    x_vector(n,:) = x_vector_hat;
end
position_GALILEO_WLMS = mean(lla);
% writeKML_GoogleEarth('Dataset2WLMS_galileo',position(1),position(2),position(3));

%% TASK A-3: Compute the standard deviation of the position error over time n=1,2,..,N
sigma_x2 = var(x_vector(:,1));
sigma_y2 = var(x_vector(:,2));
sigma_z2 = var(x_vector(:,3));
sigma_but2 = var(x_vector(:,4));
sigma_x_WLMS_galileo = sqrt(sigma_x2 + sigma_y2 + sigma_z2);

%GALILEO with LMS

K = 9; 
lla = [];

for n = 1:N
    rho = [];
    satellitePosition = [];
    rho_hat = [];
    r_hat = [];
    H = [];
    i = 1;
    for j = 1:J_gal
        if isnan(RHO.GAL(j,n))
            continue;
        else
            rho(1,i) = RHO.GAL(j,n); %build the measurement vector rho_n = [rho1,n .... rhoJ,n] of visible satellites
            satellitePosition(i,:) = SAT_POS_ECEF.GAL(j).pos(n,:); %retrieve the corresponding satellites coordinates yj,n = [xj,n yj,n zj,n]
            i = i + 1;
        end
    end
    
    x_vector_hat = zeros(1,4); %inizialization of the estimated position/state
    for k = 1:K
        for j = 1:size(satellitePosition,1)
            %compute PVT solution 
            r_hat(j) = sqrt((satellitePosition(j,1)-x_vector_hat(1,1))^2 + (satellitePosition(j,2)-x_vector_hat(1,2))^2 + ...
                (satellitePosition(j,3)-x_vector_hat(1,3))^2); %compute the euclidean distance between the estimated user point and the satellite
            rho_hat(1,j) = r_hat(j) + x_vector_hat(1,4); %estimation of pseudoranges at the k-th iteration
            H(j,1) = (satellitePosition(j,1)-x_vector_hat(1,1))/r_hat(j);
            H(j,2) = (satellitePosition(j,2)-x_vector_hat(1,2))/r_hat(j);
            H(j,3) = (satellitePosition(j,3)-x_vector_hat(1,3))/r_hat(j);
            H(j,4) = 1;
        end
        deltaRho_hat = rho_hat - rho;
        deltaX_hat = inv((transpose(H)*H))*transpose(H)*deltaRho_hat';
        x_vector_hat = x_vector_hat + deltaX_hat'; %update of user location          
    end
    %convert from ECEF to LLA and create a .kml file for GoogleEarth
    lla(n,:) = ecef2lla(x_vector_hat(1:3));
    x_vector(n,:) = x_vector_hat;
end
position_GALILEO_LMS = mean(lla);
% writeKML_GoogleEarth('Dataset2LMS_galileo',position(1),position(2),position(3));

%% TASK A-3: Compute the standard deviation of the position error over time n=1,2,..,N
sigma_x2 = var(x_vector(:,1));
sigma_y2 = var(x_vector(:,2));
sigma_z2 = var(x_vector(:,3));
sigma_but2 = var(x_vector(:,4));
sigma_x_LMS_galileo = sqrt(sigma_x2 + sigma_y2 + sigma_z2);
% 
% 
%BEIDOU with WLMS
K = 9; 
for n = 1:N
    rho = [];
    satellitePosition = [];
    rho_hat = [];
    r_hat = [];
    H = [];
    i = 1;
    for j = 1:J_bei
        if isnan(RHO.BEI(j,n))
            continue;
        else
            rho(1,i) = RHO.BEI(j,n); %build the measurement vector rho_n = [rho1,n .... rhoJ,n] of visible satellites
            satellitePosition(i,:) = SAT_POS_ECEF.BEI(j).pos(n,:); %retrieve the corresponding satellites coordinates yj,n = [xj,n yj,n zj,n]
            i = i + 1;
        end
    end
    
    x_vector_hat = zeros(1,4); %inizialization of the estimated position/state
    for k = 1:K
        for j = 1:size(satellitePosition,1)
            %compute PVT solution 
            r_hat(j) = sqrt((satellitePosition(j,1)-x_vector_hat(1,1))^2 + (satellitePosition(j,2)-x_vector_hat(1,2))^2 + ...
                (satellitePosition(j,3)-x_vector_hat(1,3))^2); %compute the euclidean distance between the estimated user point and the satellite
            rho_hat(1,j) = r_hat(j) + x_vector_hat(1,4); %estimation of pseudoranges at the k-th iteration
            H(j,1) = (satellitePosition(j,1)-x_vector_hat(1,1))/r_hat(j);
            H(j,2) = (satellitePosition(j,2)-x_vector_hat(1,2))/r_hat(j);
            H(j,3) = (satellitePosition(j,3)-x_vector_hat(1,3))/r_hat(j);
            H(j,4) = 1;
        end
        
        deltaRho_hat = rho_hat - rho;
        H_bar = inv(transpose(H)*(W_bei*H))*transpose(H)*W_bei;
        deltaX_hat = H_bar*deltaRho_hat';
        x_vector_hat = x_vector_hat + deltaX_hat'; %update of user location          
    end
    %convert from ECEF to LLA and create a .kml file for GoogleEarth
    lla(n,:) = ecef2lla(x_vector_hat(1:3));
    x_vector(n,:) = x_vector_hat;
end
position_BEIDOU_WLMS = mean(lla);
% writeKML_GoogleEarth('Dataset2WLMS_beidou',position(1),position(2),position(3));

%% TASK A-3: Compute the standard deviation of the position error over time n=1,2,..,N
sigma_x2 = var(x_vector(:,1));
sigma_y2 = var(x_vector(:,2));
sigma_z2 = var(x_vector(:,3));
sigma_but2 = var(x_vector(:,4));
sigma_x_WLMS_beidou = sqrt(sigma_x2 + sigma_y2 + sigma_z2);

%BEIDOU with LMS
K = 9; 
lla = [];

for n = 1:N
    rho = [];
    satellitePosition = [];
    rho_hat = [];
    r_hat = [];
    H = [];
    i = 1;
    for j = 1:J_bei
        if isnan(RHO.BEI(j,n))
            continue;
        else
            rho(1,i) = RHO.BEI(j,n); %build the measurement vector rho_n = [rho1,n .... rhoJ,n] of visible satellites
            satellitePosition(i,:) = SAT_POS_ECEF.BEI(j).pos(n,:); %retrieve the corresponding satellites coordinates yj,n = [xj,n yj,n zj,n]
            i = i + 1;
        end
    end
    
    x_vector_hat = zeros(1,4); %inizialization of the estimated position/state
    for k = 1:K
        for j = 1:size(satellitePosition,1)
            %compute PVT solution 
            r_hat(j) = sqrt((satellitePosition(j,1)-x_vector_hat(1,1))^2 + (satellitePosition(j,2)-x_vector_hat(1,2))^2 + ...
                (satellitePosition(j,3)-x_vector_hat(1,3))^2); %compute the euclidean distance between the estimated user point and the satellite
            rho_hat(1,j) = r_hat(j) + x_vector_hat(1,4); %estimation of pseudoranges at the k-th iteration
            H(j,1) = (satellitePosition(j,1)-x_vector_hat(1,1))/r_hat(j);
            H(j,2) = (satellitePosition(j,2)-x_vector_hat(1,2))/r_hat(j);
            H(j,3) = (satellitePosition(j,3)-x_vector_hat(1,3))/r_hat(j);
            H(j,4) = 1;
        end
        deltaRho_hat = rho_hat - rho;
        deltaX_hat = inv((transpose(H)*H))*transpose(H)*deltaRho_hat';
        x_vector_hat = x_vector_hat + deltaX_hat'; %update of user location          
    end
    %convert from ECEF to LLA and create a .kml file for GoogleEarth
    lla(n,:) = ecef2lla(x_vector_hat(1:3));
    x_vector(n,:) = x_vector_hat;
end
position_BEIDOU_LMS = mean(lla);
% writeKML_GoogleEarth('Dataset2LMS_beidou',position(1),position(2),position(3));

%% TASK A-3: Compute the standard deviation of the position error over time n=1,2,..,N
sigma_x2 = var(x_vector(:,1));
sigma_y2 = var(x_vector(:,2));
sigma_z2 = var(x_vector(:,3));
sigma_but2 = var(x_vector(:,4));
sigma_x_LMS_beidou = sqrt(sigma_x2 + sigma_y2 + sigma_z2);
