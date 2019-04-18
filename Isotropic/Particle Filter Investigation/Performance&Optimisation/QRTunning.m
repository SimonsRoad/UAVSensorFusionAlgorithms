function [meanRMSEvalues, meanRMSE_finalValues, Q, R, P_init] = QR_Tunning( plotting, nIterations)
%PERFORMANCE_CHECKER This function will run the three filters nIterations
%times and will calculate the RMSE.
%  This function is temporary so no effort in proper and efficient coding is done.
%  INPUT
% plotting - Boolean to determine if we want to plot the proccess while simulating
% nIterations - Number of iterations to run
% filter_function - A function handler pointing to the filter rhat should
% be used.
% OUTPUT 
% meanRMSEvalues Is the mean of the RMSE of all the iterations. Used to measure general performance
% meanRMSE_finalValues is the mean of the RMSE of the last iterations. Used to measure the quality of the convergence.
% EXAMPLE OF USE
% [meanRMSEvalues, meanRMSE_finalValues] = tune_filter( false, 10, @Main_isotropic_EKF)
% TODO: Allow tuning of  P_init (Initial covariance matrix)

%   Workspace cleaning
%clc; close all; clear all;

%% Initialise variables
x_state = zeros(2, 1801, nIterations);
RMSE = zeros(nIterations, 1801);

% P is same in every use case
P_init = sqrt(4000);

% Global variables needed for place_jammer() and place_uav()
global x_bnd y_bnd d2r
x_bnd=12*10^3;                                                      %   x area boundary [m]
y_bnd=12*10^3;                                                      %   y area boundary [m]
d2r=pi/180;                                                         %   Value in rad = Value in deg * d2r
% Jammer and UAV variables
x_jammer = zeros(1,2,nIterations);
x_uav = zeros(1,2,nIterations);
psi_uav = zeros(1,1,nIterations);

for i=1:nIterations
  x_jammer(:,:,i) = place_jammer();
  [x_uav(:,:,i), psi_uav(:,:,i)] = place_uav();
end

%% Start iterations for Q Optimization
%==========================================================================

% Values over we want to iterate
% StartValue of R
Q=0;
% Min, Max and Step msize for the iteration
qMin=0;
qStep=0.1;
qMax=1;
for q = qMin : qStep : qMax
    Q = [Q, q];
end
% In this case is Q static, but needs to be adapted on the dynamic sze of R
R = 0.01;
Rq_s=0.01
for m=1:(size(Q,2)-1)  
    R = [R, 0.01];
end

meanRMSEvalues = zeros(size(Q,2), size(R,2));
meanRMSE_finalValues = zeros(size(Q,2), size(R,2));

% In this case is Q static, but needs to be adapted on the dynamic sze of R
for m=1:size(Q,2)
    for k=1:size(R,2)
        for i=1:nIterations
            [x_state(:, :, i), x_t_vec] = Main_isotropic_EKF(plotting, Q(m), R(k), x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i),P_init);
            for l=1 : size(x_state, 2)
                RMSE(i, l) = norm((x_state(:,l, i)- x_t_vec'));
            end
        end
        meanRMSEvalues(m,k) = mean(mean(RMSE));
        meanRMSE_finalValues(m,k) = mean(mean(RMSE(:, 1500:end)));
        RMSE = zeros(nIterations, 1801);
    end
end

% Mean of overall simualation
meanSimValues=mean(meanRMSEvalues')

%%Ploting of Q related results
mkdir('tuning_results')
h = figure;
%surf(Q,R, meanRMSEvalues')
plot(Q, meanSimValues)
set(gca, 'XScale', 'log', 'YScale', 'log', 'ZLim', [0 2000], 'FontSize', 10)
caxis([0 2000])
colormap jet
xlabel('Q')
ylabel('RMSE (m)')
title(sprintf('RMSE vs Q\n R = %.2f and P_init = %d \n Random Samples = %d',Rq_s ,P_init,nIterations))
shading interp

print('./tuning_results/Q','-dsvg')

%%End Q Optimization
%==========================================================================



%% Start iterations for R Optimization
%==========================================================================
% reset of reused values
clear Q R meanRMSEvalues meanRMSE_finalValues
% Values over we want to iterate
% StartValue of R
R=0;
% Min, Max and Step msize for the iteration
rMin=0;
rStep=1;
rMax=10;
% Ieration loop for R
for r = rMin : rStep : rMax
    R = [R, r];
end
% In this case is Q static, but needs to be adapted on the dynamic sze of R
Q = 0.01;
Qr_s=0.01
for m=1:(size(R,2)-1)  
    Q = [Q, 0.01];
end

% Iteration Loop for optimization
for m=1:size(Q,2)
    for k=1:size(R,2)
        for i=1:nIterations
            [x_state(:, :, i), x_t_vec] = Main_isotropic_EKF(plotting, Q(m), R(k), x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i),P_init);
            for l=1 : size(x_state, 2)
                RMSE(i, l) = norm((x_state(:,l, i)- x_t_vec'));
            end
        end
        meanRMSEvalues(m,k) = mean(mean(RMSE));
        meanRMSE_finalValues(m,k) = mean(mean(RMSE(:, 1500:end)));
        RMSE = zeros(nIterations, 1801);
    end
end

% Mean of overall simualation
meanSimValues=mean(meanRMSEvalues')

%%Ploting of Q related results
h = figure;
%surf(Q,R, meanRMSEvalues')
plot(R, meanSimValues)
set(gca, 'XScale', 'log', 'YScale', 'log', 'ZLim', [0 2000], 'FontSize', 10)
caxis([0 2000])
colormap jet
xlabel('R')
ylabel('RMSE (m)')
%zlabel('RMSE (m)')
title(sprintf('RMSE vs R\n R = %.2f and P_init = %d \n Random Samples = %d',Rq_s ,P_init,nIterations))
shading interp

print('./tuning_results/R','-dsvg')

%%End Q Optimization
%==========================================================================

end

