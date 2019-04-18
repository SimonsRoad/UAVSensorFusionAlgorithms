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

% RMSE for each filter and for each t of all the iterations
RMSE_EKF = zeros(nIterations, 1801);
CRLB_EKF = zeros(nIterations, 1801);
RMSE_UKF = zeros(nIterations, 1801);
CRLB_UKF = zeros(nIterations, 1801);
RMSE_PF = zeros(nIterations, 1801);
CRLB_PF = zeros(nIterations, 1801);

% Computing Time for each filter and for each t of all the iterations
TIME_EKF = zeros(nIterations,1);
TIME_UKF = zeros(nIterations,1);
TIME_PF = zeros(nIterations,1);


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
N = 1000;

%All
Pinit =sqrt(4000);

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
CRLB_PF = mean(CRLB_PF);
TIME_PF = mean(TIME_PF)

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

%% Comparison plots between filters -- EKF vs UKF vs PF
clear plot;
figure('Visible','off');
plot(RMSE_EKF, 'r');
hold on;
plot(RMSE_UKF, 'b');
plot(RMSE_PF , 'g');
legend ('EKF','UKF','PF')
xlabel('Time step')
ylabel('RMSE (m)')
print('./Graphs/All','-dsvg')
% Final part
xlim([1000, 1800])
ylim('auto')
title('Final track (m)')
print('./Graphs/All_final','-dsvg')
end

