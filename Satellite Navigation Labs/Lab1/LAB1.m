%%
close all
clear all
clc

%%
%load('dataset_1_20180328T122038.mat');
%load('dataset_2_20180328T122158.mat');
%load('dataset_3_20180328T121914.mat');
%load('dataset_4_20180328T121804.mat');
load('dataset_5_20180328T121529.mat');
%load('dataset_6_20180328T121701.mat');

N = size(RHO.GPS,2); %3600
J_gps = size(RHO.GPS,1);
J_gal = size(RHO.GAL,1);
J_bei = size(RHO.BEI,1);
J_glo = size(RHO.GLO,1);

%% TASK A-1: Visibility and pseudoranges GPS system



sat_gps = zeros(1,N);
vis_gps = figure; hold on
for i = 1:N
    for j = 1:J_gps
        if isnan(RHO.GPS(j,i))
            %gps(j,i) = 0;
            continue;
        else
            sat_gps(1,i) = sat_gps(1,i) + 1;
        end
    end
end
plot(sat_gps(1,:),'.b','Linewidth',1.5);
hold off
grid on
axis([0 3600 0 7])
yticks(0:1:7);
ylabel('Number of visible satellites')
xlabel('Instants (seconds)')
title('Visibility of GPS satellites')
figure(vis_gps)

pseud_gps = figure; 
color = jet(32);
for i =1:J_gps
    if isnan(RHO.GPS(i,:))
        continue;
    else
        txt = ['Satellite ',num2str(i)];
        p = plot(RHO.GPS(i,:),'-','Linewidth',1.5,'Color',color(i,:),'DisplayName', txt);
        hold on
    end
end
hold off
grid on
axis([0 3600 1.8e07 3e07])
ylabel('Pseudorange')
xlabel('Instants (seconds)')
title('Pseudoranges of GPS satellites')
legend('NumColumns',2,'Location','best')
figure(pseud_gps)

%% Visibility and pseudoranges BEIDOU system 
 

sat_bei = zeros(1,N);
vis_bei = figure; hold on
for i = 1:N
    for j = 1:J_bei
        if isnan(RHO.BEI(j,i))
            %bei(j,i) = 0;
            continue;
        else
            sat_bei(1,i) = sat_bei(1,i) + 1;
        end
    end
end
plot(sat_bei(1,:),'.b','Linewidth',1.5);
hold off
grid on
axis([0 3600 0 J_bei])
yticks(0:1:J_bei);
ylabel('Number of visible satellites')
xlabel('Instants (seconds)')
title('Visibility of BEIDOU satellites')
figure(vis_bei)

pseud_bei = figure; 
rng(1);
color = jet(14);
color = color(randperm(size(color, 1)), :);
for i =1:J_bei
    if isnan(RHO.BEI(i,:))
        continue;
    elseif i == 5
        txt = ['Satellite ',num2str(i)];
        p = plot(RHO.BEI(i,:),'-','Linewidth',1.5,'Color','k','DisplayName', txt);
        hold on
    else
        txt = ['Satellite ',num2str(i)];
        p = plot(RHO.BEI(i,:),'-','Linewidth',1.5,'Color',color(i,:),'DisplayName', txt);
        hold on
    end
end
hold off
grid on
axis([0 3600 2.5e07 4.5e07])
ylabel('Pseudorange')
xlabel('Instants (seconds)')
title('Pseudoranges of BEIDOU satellites')
legend('NumColumns',2,'Location','best')
figure(pseud_bei)

%% Visibility GALILEO system



sat_gal = zeros(1,N);
vis_gal = figure; hold on
for i = 1:N
    for j = 1:J_gal
        if isnan(RHO.GAL(j,i))
            %gal(j,i) = 0;
            continue;
        else
            sat_gal(1,i) = sat_gal(1,i) + 1;
        end
    end
end
plot(sat_gal(1,:),'.b','Linewidth',1.5);
hold off
grid on
axis([0 3600 0 J_gal])
yticks(0:1:J_gal);
ylabel('Number of visible satellites')
xlabel('Instants (seconds)')
title('Visibility of GALILEO satellites')
figure(vis_gal)

pseud_gal = figure; 
rng(2);
color = jet(30);
color = color(randperm(size(color, 1)), :);
for i =1:J_gal
    if isnan(RHO.GAL(i,:))
        continue;
    else
        txt = ['Satellite ',num2str(i)];
        p = plot(RHO.GAL(i,:),'-','Linewidth',1.5,'Color',color(i,:),'DisplayName', txt);
        hold on
    end
end
hold off
grid on
axis([0 3600 2e07 3.2e07])
ylabel('Pseudorange')
xlabel('Instants (seconds)')
title('Pseudoranges of GALILEO satellites')
legend('NumColumns',2,'Location','best')
figure(pseud_gal)

%% Visibility and pseudoranges GLONASS system


sat_glo = zeros(1,N);
vis_glo = figure; hold on
for i = 1:N
    for j = 1:J_glo
        if isnan(RHO.GLO(j,i))
            %glo(j,i) = 0;
            continue;
        else
            sat_glo(1,i) = sat_glo(1,i) + 1;
        end
    end
end
plot(sat_glo(1,:),'.b','Linewidth',1.5);
hold off
grid on
axis([0 3600 0 20])
yticks(0:1:20);
ylabel('Number of visible satellites')
xlabel('Instants (seconds)')
title('Visibility of GLONASS satellites')
figure(vis_glo)

pseud_glo = figure; 
rng(2);
color = jet(54);
color = color(randperm(size(color, 1)), :);
for i =1:J_glo
    if isnan(RHO.GLO(i,:))
        continue;
    else
        txt = ['Satellite ',num2str(i)];
        p = plot(RHO.GLO(i,:),'-','Linewidth',1.5,'Color',color(i,:),'DisplayName', txt);
        hold on
    end
end
hold off
grid on
axis([0 3600 1.55e07 3e07])
ylabel('Pseudorange')
xlabel('Instants (seconds)')
title('Pseudoranges of GLONASS satellites')
legend('NumColumns',2,'Location','best')
figure(pseud_glo)

%% TASK A-2: Develop a LMS positioning algorithm and estimate the user state at each time instant n &&


%GPS
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
        deltaX_hat = inv((transpose(H)*H))*transpose(H)*deltaRho_hat';
        x_vector_hat = x_vector_hat + deltaX_hat'; %update of user location          
    end
    %convert from ECEF to LLA and create a .kml file for GoogleEarth
    lla(n,:) = ecef2lla(x_vector_hat(1:3));
    x_vector(n,:) = x_vector_hat;
    G_matrix = inv((transpose(H)*H));
    GDOP(n) = sqrt(trace(G_matrix));
end
position = mean(lla);
GDOP_gps = mean(GDOP);

% writeKML_GoogleEarth('Dataset2_GPS',position(1),position(2),position(3));

%% TASK A-3: Compute the standard deviation of the position error over time n=1,2,..,N
sigma_x2 = var(x_vector(:,1));
sigma_y2 = var(x_vector(:,2));
sigma_z2 = var(x_vector(:,3));
sigma_but2 = var(x_vector(:,4));
sigma_x_gps = sqrt(sigma_x2 + sigma_y2 + sigma_z2);

%GLONASS
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
        deltaX_hat = inv((transpose(H)*H))*transpose(H)*deltaRho_hat';
        x_vector_hat = x_vector_hat + deltaX_hat'; %update of user location          
    end
    %convert from ECEF to LLA and create a .kml file for GoogleEarth
    lla(n,:) = ecef2lla(x_vector_hat(1:3));
    x_vector(n,:) = x_vector_hat;
    GDOP(n) = sqrt(trace(inv((transpose(H)*H))));
end
position = mean(lla);
GDOP_glonass = mean(GDOP);
% writeKML_GoogleEarth('Dataset2_GLONASS',position(1),position(2),position(3));
% 
%% TASK A-3: Compute the standard deviation of the position error over time n=1,2,..,N
sigma_x2 = var(x_vector(:,1));
sigma_y2 = var(x_vector(:,2));
sigma_z2 = var(x_vector(:,3));
sigma_but2 = var(x_vector(:,4));
sigma_x_glonass = sqrt(sigma_x2 + sigma_y2 + sigma_z2);

GALILEO
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
            compute PVT solution 
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
    convert from ECEF to LLA and create a .kml file for GoogleEarth
    lla(n,:) = ecef2lla(x_vector_hat(1:3));
    x_vector(n,:) = x_vector_hat;
    GDOP(n) = sqrt(trace(inv((transpose(H)*H))));
end
position = mean(lla);
GDOP_galileo = mean(GDOP);
% writeKML_GoogleEarth('Dataset5',position(1),position(2),position(3));

% TASK A-3: Compute the standard deviation of the position error over time n=1,2,..,N
sigma_x2 = var(x_vector(:,1));
sigma_y2 = var(x_vector(:,2));
sigma_z2 = var(x_vector(:,3));
sigma_but2 = var(x_vector(:,4));
sigma_x_galileo = sqrt(sigma_x2 + sigma_y2 + sigma_z2);

%BEIDOU
K = 9; 
H =[];
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
    GDOP(n) = sqrt(trace(inv((transpose(H)*H))));
end
position = mean(lla);
GDOP_beidou = mean(GDOP);
% writeKML_GoogleEarth('Dataset5',position(1),position(2),position(3));

%% TASK A-3: Compute the standard deviation of the position error over time n=1,2,..,N
sigma_x2 = var(x_vector(:,1));
sigma_y2 = var(x_vector(:,2));
sigma_z2 = var(x_vector(:,3));
sigma_but2 = var(x_vector(:,4));
sigma_x_beidou = sqrt(sigma_x2 + sigma_y2 + sigma_z2);

