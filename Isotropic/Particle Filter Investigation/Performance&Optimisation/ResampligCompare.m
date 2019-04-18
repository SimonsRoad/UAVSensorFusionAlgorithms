function [RMSE_NM, CRLB_NM, RMSE_SYS, CRLB_SYS] = PerformanceCompare( plotting, nIterations)
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

x_state_nm = zeros(2, 1801, nIterations);
x_state_sys = zeros(2, 1801, nIterations);
x_state_lse = zeros(2, 1801, nIterations);
x_state_res = zeros(2, 1801, nIterations);
x_state_rsr = zeros(2, 1801, nIterations);
x_state_srt = zeros(2, 1801, nIterations);

% RMSE for each filter and for each t of all the iterations

RMSE_SYS = zeros(nIterations, 1801);
CRLB_SYS = zeros(nIterations, 1801);

RMSE_NM = zeros(nIterations, 1801);
CRLB_NM = zeros(nIterations, 1801);

RMSE_LSE = zeros(nIterations, 1801);
CRLB_LSE = zeros(nIterations, 1801);

RMSE_RES = zeros(nIterations, 1801);
CRLB_RES = zeros(nIterations, 1801);

RMSE_RSR = zeros(nIterations, 1801);
CRLB_RSR = zeros(nIterations, 1801);

RMSE_SRT = zeros(nIterations, 1801);
CRLB_SRT = zeros(nIterations, 1801);

% Computing Time for each filter and for each t of all the iterations

TIME_SYS = zeros(nIterations,1);
TIME_NM = zeros(nIterations,1);
TIME_LSE = zeros(nIterations,1);
TIME_RES = zeros(nIterations,1);
TIME_RSR = zeros(nIterations,1);
TIME_SRT = zeros(nIterations,1);


% Jammer and UAV variables
x_jammer = zeros(1,2,nIterations);
x_uav = zeros(1,2,nIterations);
psi_uav = zeros(1,1,nIterations);

%SYS
Q_SYS = 0.01;
R_SYS = 0.08;

%NM
Q_NM = 0.01;
R_NM = 0.08;

%NM
Q_LSE = 0.01;
R_LSE = 0.08;

%NM
Q_RES = 0.01;
R_RES = 0.08;

%NM
Q_RSR = 0.01;
R_RSR = 0.08;

%NM
Q_SRT = 0.01;
R_SRT = 0.08;

% Number of Particles
N = 500;

%All
Pinit =sqrt(4000);

%%Samples creation which are used by all algorithms
for i=1:nIterations
  x_jammer(:,:,i) = place_jammer();
  [x_uav(:,:,i), psi_uav(:,:,i)] = place_uav();
end

for i=1:nIterations
    [x_state_sys(:, :, i), x_t_vec, P_cov_SYS(:,:,:,i),TIME_SYS(i)] = Main_isotropic_PF_SYS(plotting, Q_SYS, R_SYS, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i), Pinit,N);
    % For each t calculate the distance between the estimate and the real location
    for k=1 : size(x_state_sys, 2)
        % Calculate the norm of all the states 
        RMSE_SYS(i, k) = norm((x_state_sys(:,k, i)- x_t_vec'));
        CRLB_SYS(i, k) = norm(diag(sqrt(diag(P_cov_SYS(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end

% Calculate the average error along all iterations for each t in the process
RMSE_SYS = mean(RMSE_SYS);
ERROR_SYS = RMSE_SYS(1800)
CRLB_SYS = mean(CRLB_SYS);
TIME_SYS = mean(TIME_SYS)

for i=1:nIterations
    [x_state_nm(:, :, i), x_t_vec_NM, P_cov_NM(:,:,:,i),TIME_NM(i)] = Main_isotropic_PF_NM(plotting, Q_NM, R_NM, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i), Pinit,N);
    for k=1 : size(x_state_nm, 2)
        RMSE_NM(i, k) = norm((x_state_nm(:,k, i)- x_t_vec_NM'));
        CRLB_NM(i, k) = norm(diag(sqrt(diag(P_cov_NM(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end
RMSE_NM = mean(RMSE_NM);
ERROR_NM = RMSE_NM(1800)
CRLB_NM = mean(CRLB_NM);
TIME_NM = mean(TIME_NM)

for i=1:nIterations
    [x_state_lse(:, :, i), x_t_vec_LSE, P_cov_LSE(:,:,:,i),TIME_LSE(i)] = Main_isotropic_PF_LSE(plotting, Q_LSE, R_LSE, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i), Pinit,N);
    for k=1 : size(x_state_lse, 2)
        RMSE_LSE(i, k) = norm((x_state_lse(:,k, i)- x_t_vec_LSE'));
        CRLB_LSE(i, k) = norm(diag(sqrt(diag(P_cov_LSE(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end
RMSE_LSE = mean(RMSE_LSE);
ERROR_LSE = RMSE_LSE(1800)
CRLB_LSE = mean(CRLB_LSE);
TIME_LSE = mean(TIME_LSE)

for i=1:nIterations
    [x_state_res(:, :, i), x_t_vec_RES, P_cov_RES(:,:,:,i),TIME_RES(i)] = Main_isotropic_PF_RES(plotting, Q_RES, R_RES, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i), Pinit,N);
    for k=1 : size(x_state_res, 2)
        RMSE_RES(i, k) = norm((x_state_res(:,k, i)- x_t_vec_RES'));
        CRLB_RES(i, k) = norm(diag(sqrt(diag(P_cov_RES(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end
RMSE_RES = mean(RMSE_RES);
ERROR_RES = RMSE_RES(1800)
CRLB_RES = mean(CRLB_RES);
TIME_RES = mean(TIME_RES)

for i=1:nIterations
    [x_state_rsr(:, :, i), x_t_vec_RSR, P_cov_RSR(:,:,:,i),TIME_RSR(i)] = Main_isotropic_PF_RSR(plotting, Q_RSR, R_RSR, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i), Pinit,N);
    for k=1 : size(x_state_rsr, 2)
        RMSE_RSR(i, k) = norm((x_state_rsr(:,k, i)- x_t_vec_RSR'));
        CRLB_RSR(i, k) = norm(diag(sqrt(diag(P_cov_RSR(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end
RMSE_RSR = mean(RMSE_RSR);
ERROR_RSR = RMSE_RSR(1800)
CRLB_RSR = mean(CRLB_RSR);
TIME_RSR = mean(TIME_RSR)

for i=1:nIterations
    [x_state_srt(:, :, i), x_t_vec_SRT, P_cov_SRT(:,:,:,i),TIME_SRT(i)] = Main_isotropic_PF_SRT(plotting, Q_SRT, R_SRT, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i), Pinit,N);
    for k=1 : size(x_state_srt, 2)
        RMSE_SRT(i, k) = norm((x_state_srt(:,k, i)- x_t_vec_SRT'));
        CRLB_SRT(i, k) = norm(diag(sqrt(diag(P_cov_SRT(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end
RMSE_SRT = mean(RMSE_SRT);
ERROR_SRT = RMSE_SRT(1800)
CRLB_SRT = mean(CRLB_SRT);
TIME_SRT = mean(TIME_SRT)


%% Individual plots for each filter -- Filter vs Cramer Rao Lower Bound
mkdir('Resampling')

% SYS
CRLB_SYS(1,1:25) = NaN;
figure('Visible','off')
plot(RMSE_SYS)
hold on
plot(CRLB_SYS, 'g')
legend('RMSE', 'CRLB')
title('SYS RMSE (m)')
xlabel('Time step')
ylabel('RMSE (m)')
print('./Resampling/SYS','-dsvg')

% SYS-Final
xlim([1000, 1800])
ylim('auto')
title('Final track PF RMSE (m)')
print('./Resampling/SYS_final','-dsvg')

% NM
CRLB_NM(1,1:25) = NaN;
figure('Visible','off')
plot(RMSE_NM)
hold on
plot(CRLB_NM, 'g')
legend('RMSE', 'CRLB')
title('NM RMSE (m)')
xlabel('Time step')
ylabel('RMSE (m)')
print('./Resampling/NM','-dsvg')

% NM-Final
xlim([1000, 1800])
ylim('auto')
title('Final track PF RMSE (m)')
print('./Resampling/NM_final','-dsvg')


%% Comparison plots between filters -- EKF vs UKF vs PF
clear plot;
figure('Visible','off');
plot(RMSE_NM);
hold on;
plot(RMSE_SYS);
plot(RMSE_LSE);
plot(RMSE_RES);
plot(RMSE_RSR);
plot(RMSE_SRT);
legend ('Multinomial','Systematic','Local-Selection','Residual','Residual-Systematic','Stratified')
xlabel('Time step')
ylabel('RMSE (m)')
print('./Resampling/All','-dsvg')
% Final part
xlim([1000, 1800])
ylim('auto')
title('Final track (m)')
print('./Resampling/All_final','-dsvg')
end

