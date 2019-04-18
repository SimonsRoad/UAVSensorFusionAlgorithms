function [RMSE_EKF, CRLB_EKF, RMSE_UKF, CRLB_UKF, RMSE_PF, CRLB_PF] = PerformanceCompare( plotting, nIterations)
%PERFORMANCE_CHECKER This function will run the three filters nIterations
% times, it will calculate the RMSE and plot the results in svg format.
% This function is temporary so no effort in proper and efficient coding is done.
% Example of use [RMSE_EKF, CRLB_EKF, RMSE_UKF, CRLB_UKF, RMSE_PF, CRLB_PF] = measure_performance(false, 10)
% TODO: This function should be run after tuning individual Q and R for
% each filter.
% TODO2: Add possibility of modifying P_init


%% Initialisation of variables for each filter
% Global variables needed for place_jammer() and place_uav()
global x_bnd y_bnd d2r
x_bnd=12*10^3;                                                      %   x area boundary [m]
y_bnd=12*10^3;                                                      %   y area boundary [m]
d2r=pi/180;                                                         %   Value in rad = Value in deg * d2r
% Matrix for storing each state (2) at each increment of time (1801) for
% every iteraion (nIterations)
x_state_ekf = zeros(2, 1801, nIterations);
x_state_ukf = zeros(2, 1801, nIterations);
x_state_pf = zeros(2, 1801, nIterations);
x_state_akf = zeros(2, 1801, nIterations);
x_state_hinf = zeros(2, 1801, nIterations);
x_state_ahinf = zeros(2, 1801, nIterations);
x_state_upf = zeros(2, 1801, nIterations);
x_state_epf = zeros(2, 1801, nIterations);

% RMSE for each filter and for each t of all the iterations
RMSE_EKF = zeros(nIterations, 1801);
CRLB_EKF = zeros(nIterations, 1801);

RMSE_UKF = zeros(nIterations, 1801);
CRLB_UKF = zeros(nIterations, 1801);

RMSE_PF = zeros(nIterations, 1801);
CRLB_PF = zeros(nIterations, 1801);

RMSE_AKF = zeros(nIterations, 1801);
CRLB_AKF = zeros(nIterations, 1801);

RMSE_HINF = zeros(nIterations, 1801);
CRLB_HINF = zeros(nIterations, 1801);

RMSE_AHINF = zeros(nIterations, 1801);
CRLB_AHINF = zeros(nIterations, 1801);

RMSE_EPF = zeros(nIterations, 1801);
CRLB_EPF = zeros(nIterations, 1801);

RMSE_UPF = zeros(nIterations, 1801);
CRLB_UPF = zeros(nIterations, 1801);

% Computing Time for each filter and for each t of all the iterations
TIME_EKF = zeros(nIterations,1);
TIME_UKF = zeros(nIterations,1);
TIME_PF = zeros(nIterations,1);
TIME_AKF = zeros(nIterations,1);
TIME_HINF = zeros(nIterations,1);
TIME_AHINF = zeros(nIterations,1);
TIME_EPF = zeros(nIterations,1);
TIME_UPF = zeros(nIterations,1);


% Jammer and UAV variables
x_jammer = zeros(1,2,nIterations);
x_uav = zeros(1,2,nIterations);
psi_uav = zeros(1,1,nIterations);

%EKF
Q_EKF = 0.1;
R_EKF = 0.08;
%UKF
Q_UKF = 0.1;
R_UKF = 0.08;

%PF
Q_PF = 5;
R_PF = 0.08;

%AKF
Q_AKF = 0.1;
R_AKF = 0.08;

%HINF
Q_HINF = 0.1;
R_HINF = 0.08;

%AHINF
Q_AHINF = 0.1;
R_AHINF = 0.08;

%EPF
Q_EPF = 5;
R_EPF = 0.08;

%UPF
Q_UPF = 0.1;
R_UPF = 0.08;

% Number of Particles
N = 1000;

%All
Pinit =sqrt(4000);
Pupf =sqrt(400);

%%Samples creation which are used by all algorithms
for i=1:nIterations
  x_jammer(:,:,i) = place_jammer();
  [x_uav(:,:,i), psi_uav(:,:,i)] = place_uav();
end

for i=1:nIterations
    [x_state_ekf(:, :, i), x_t_vec, P_cov_EKF(:,:,:,i),TIME_EKF(i)] = Main_isotropic_EKF(plotting, Q_EKF, R_EKF, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i), Pinit);
    % For each t calculate the distance between the estimate and the real location
    for k=1 : size(x_state_ekf, 2)
        % Calculate the norm of all the states 
        RMSE_EKF(i, k) = norm((x_state_ekf(:,k, i)- x_t_vec'));
        CRLB_EKF(i, k) = norm(diag(sqrt(diag(P_cov_EKF(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end

% Calculate the average error along all iterations for each t in the process
RMSE_EKF = mean(RMSE_EKF);
ERROR_EKF = RMSE_EKF(1800)
CRLB_EKF = mean(CRLB_EKF);
TIME_EKF = mean(TIME_EKF)

for i=1:nIterations
    [x_state_ukf(:, :, i), x_t_vec_UKF, P_cov_UKF(:,:,:,i),TIME_UKF(i)] = Main_isotropic_UKF(plotting, Q_UKF, R_UKF, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i), Pinit);
    for k=1 : size(x_state_ukf, 2)
        RMSE_UKF(i, k) = norm((x_state_ukf(:,k, i)- x_t_vec_UKF'));
        CRLB_UKF(i, k) = norm(diag(sqrt(diag(P_cov_UKF(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end
RMSE_UKF = mean(RMSE_UKF);
ERROR_UKF = RMSE_UKF(1800)
CRLB_UKF = mean(CRLB_UKF);
TIME_UKF = mean(TIME_UKF)

for i=1:nIterations
    [x_state_pf(:, :, i), x_t_vec, P_cov(:,:,:,i),TIME_PF(i)] = Main_isotropic_PF(plotting, Q_PF, R_PF, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i), Pinit,N);
    for k=1 : size(x_state_pf, 2)
        RMSE_PF(i, k) = norm((x_state_pf(:,k, i)- x_t_vec'));
        CRLB_PF(i, k) = norm(diag(sqrt(diag(P_cov(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end
RMSE_PF = mean(RMSE_PF);
ERROR_PF = RMSE_PF(1800)
CRLB_PF = mean(CRLB_PF);
TIME_PF = mean(TIME_PF)

for i=1:nIterations
    [x_state_akf(:, :, i), x_t_vec, P_cov(:,:,:,i),TIME_AKF(i)] = Main_isotropic_AKF(plotting, Q_AKF, R_AKF, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i), Pinit);
    for k=1 : size(x_state_akf, 2)
        RMSE_AKF(i, k) = norm((x_state_akf(:,k, i)- x_t_vec'));
        CRLB_AKF(i, k) = norm(diag(sqrt(diag(P_cov(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end
RMSE_AKF = mean(RMSE_AKF);
ERROR_AKF = RMSE_AKF(1800)
CRLB_AKF = mean(CRLB_AKF);
TIME_AKF = mean(TIME_AKF)

for i=1:nIterations
    [x_state_hinf(:, :, i), x_t_vec, P_cov(:,:,:,i),TIME_HINF(i)] = Main_isotropic_HINF(plotting, Q_HINF, R_HINF, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i), Pinit);
    for k=1 : size(x_state_hinf, 2)
        RMSE_HINF(i, k) = norm((x_state_hinf(:,k, i)- x_t_vec'));
        CRLB_HINF(i, k) = norm(diag(sqrt(diag(P_cov(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end
RMSE_HINF = mean(RMSE_HINF);
ERROR_HINF = RMSE_HINF(1800)
CRLB_HINF = mean(CRLB_HINF);
TIME_HINF = mean(TIME_HINF)

for i=1:nIterations
    [x_state_ahinf(:, :, i), x_t_vec, P_cov(:,:,:,i),TIME_AHINF(i)] = Main_isotropic_AHINF(plotting, Q_AHINF, R_AHINF, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i), Pinit);
    for k=1 : size(x_state_ahinf, 2)
        RMSE_AHINF(i, k) = norm((x_state_ahinf(:,k, i)- x_t_vec'));
        CRLB_AHINF(i, k) = norm(diag(sqrt(diag(P_cov(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end
RMSE_AHINF = mean(RMSE_AHINF);
ERROR_AHINF = RMSE_AHINF(1800)
CRLB_AHINF = mean(CRLB_AHINF);
TIME_AHINF = mean(TIME_AHINF)

for i=1:nIterations
    [x_state_epf(:, :, i), x_t_vec, P_cov(:,:,:,i),TIME_EPF(i)] = Main_isotropic_EPF(plotting, Q_EPF, R_EPF, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i), Pinit,N);
    for k=1 : size(x_state_epf, 2)
        RMSE_EPF(i, k) = norm((x_state_epf(:,k, i)- x_t_vec'));
        CRLB_EPF(i, k) = norm(diag(sqrt(diag(P_cov(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end
RMSE_EPF = mean(RMSE_EPF);
ERROR_EPF = RMSE_EPF(1800)
CRLB_EPF = mean(CRLB_EPF);
TIME_EPF = mean(TIME_EPF)

for i=1:nIterations
    [x_state_upf(:, :, i), x_t_vec, P_cov(:,:,:,i),TIME_UPF(i)] = Main_isotropic_UPF(plotting, Q_UPF, R_UPF, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i), Pupf,N);
    for k=1 : size(x_state_upf, 2)
        RMSE_UPF(i, k) = norm((x_state_upf(:,k, i)- x_t_vec'));
        CRLB_UPF(i, k) = norm(diag(sqrt(diag(P_cov(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end
RMSE_UPF = mean(RMSE_UPF);
ERROR_UPF = RMSE_UPF(1800)
CRLB_UPF = mean(CRLB_UPF);
TIME_UPF = mean(TIME_UPF)


%% Individual plots for each filter -- Filter vs Cramer Rao Lower Bound
mkdir('Graphs')

% EKF
CRLB_EKF(1,1:25) = NaN;
figure('Visible','off')
plot(RMSE_EKF)
hold on
plot(CRLB_EKF, 'g')
legend('RMSE', 'CRLB')
title('EKF RMSE (m)')
xlabel('Time step')
ylabel('RMSE (m)')
print('./Graphs/EKF','-dsvg')

% EKF-Final
xlim([1000, 1800])
ylim('auto')
title('Final track EKF RMSE (m)')
print('./Graphs/EKF_final','-dsvg')

% UKF
CRLB_UKF(1,1:25) = NaN;
figure('Visible','off')
plot(RMSE_UKF)
hold on
plot(CRLB_UKF, 'g')
legend('RMSE', 'CRLB')
title('UKF RMSE (m)')
xlabel('Time step')
ylabel('RMSE (m)')
print('./Graphs/UKF','-dsvg')

% UKF-Final
xlim([1000, 1800])
ylim('auto')
title('Final track UKF RMSE (m)')
print('./Graphs/UKF_final','-dsvg')

% PF
CRLB_PF(1,1:25) = NaN;
figure('Visible','off')
plot(RMSE_PF)
hold on
plot(CRLB_PF, 'g')
legend('RMSE', 'CRLB')
title('PF RMSE (m)')
xlabel('Time step')
ylabel('RMSE (m)')
print('./Graphs/PF','-dsvg')

% PF-Final
xlim([1000, 1800])
ylim('auto')
title('Final track PF RMSE (m)')
print('./Graphs/PF_final','-dsvg')

% AKF
CRLB_AKF(1,1:25) = NaN;
figure('Visible','off')
plot(RMSE_AKF)
hold on
plot(CRLB_AKF, 'g')
legend('RMSE', 'CRLB')
title('AKF RMSE (m)')
xlabel('Time step')
ylabel('RMSE (m)')
print('./Graphs/AKF','-dsvg')

% AKF-Final
xlim([1000, 1800])
ylim('auto')
title('Final track AKF RMSE (m)')
print('./Graphs/AKF_final','-dsvg')

% HINF
CRLB_HINF(1,1:25) = NaN;
figure('Visible','off')
plot(RMSE_HINF)
hold on
plot(CRLB_HINF, 'g')
legend('RMSE', 'CRLB')
title('HINF RMSE (m)')
xlabel('Time step')
ylabel('RMSE (m)')
print('./Graphs/HINF','-dsvg')

% HINF-Final
xlim([1000, 1800])
ylim('auto')
title('Final track HINF RMSE (m)')
print('./Graphs/HINF_final','-dsvg')

% AHINF
CRLB_AHINF(1,1:25) = NaN;
figure('Visible','off')
plot(RMSE_AHINF)
hold on
plot(CRLB_AHINF, 'g')
legend('RMSE', 'CRLB')
title('AHINF RMSE (m)')
xlabel('Time step')
ylabel('RMSE (m)')
print('./Graphs/AHINF','-dsvg')

% AHINF-Final
xlim([1000, 1800])
ylim('auto')
title('Final track AHINF RMSE (m)')
print('./Graphs/AHINF_final','-dsvg')

% EPF
CRLB_EPF(1,1:25) = NaN;
figure('Visible','off')
plot(RMSE_EPF)
hold on
plot(CRLB_EPF, 'g')
legend('RMSE', 'CRLB')
title('EPF RMSE (m)')
xlabel('Time step')
ylabel('RMSE (m)')
print('./Graphs/EPF','-dsvg')

% EPF-Final
xlim([1000, 1800])
ylim('auto')
title('Final track PF RMSE (m)')
print('./Graphs/EPF_final','-dsvg')

% UPF
CRLB_UPF(1,1:25) = NaN;
figure('Visible','off')
plot(RMSE_UPF)
hold on
plot(CRLB_UPF, 'g')
legend('RMSE', 'CRLB')
title('UPF RMSE (m)')
xlabel('Time step')
ylabel('RMSE (m)')
print('./Graphs/UPF','-dsvg')

% UPF-Final
xlim([1000, 1800])
ylim('auto')
title('Final track PF RMSE (m)')
print('./Graphs/UPF_final','-dsvg')


%% Comparison plots between filters -- EKF vs UKF vs PF
clear plot;
figure('Visible','off');
plot(RMSE_EKF);
hold on;
plot(RMSE_UKF);
plot(RMSE_PF );
plot(RMSE_AKF );
plot(RMSE_HINF);
plot(RMSE_AHINF);
plot(RMSE_EPF);
plot(RMSE_UPF);
legend ('EKF','UKF','PF','AKF','HINF','AHINF','EPF','UPF')
xlabel('Time step')
ylabel('RMSE (m)')
print('./Graphs/All','-dsvg')
% Final part
xlim([1000, 1800])
ylim('auto')
title('Final track (m)')
print('./Graphs/All_final','-dsvg')
end

