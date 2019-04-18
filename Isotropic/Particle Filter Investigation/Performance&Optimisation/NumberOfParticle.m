function [meanRMSEvalues, meanRMSE_finalValues, Q, R, P_init] = NumberOfParticle( plotting, nIterations)
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
close all;
%% Initialise variables
x_state = zeros(2, 1801, nIterations);
RMSE = zeros(nIterations, 1801);

% P is same in every use case
P_init = 400;

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

% Min, Max and Step msize for the iteration
nMin=100;
nStep=100;
nMax=10000;
N=nMin;
for n = (nMin+nStep) : nStep : nMax
    N = [N, n];
end

%% Statistical Properties
Q = 0.1;
R = 0.1;
P_ini = sqrt(4000);

meanRMSEvalues = zeros(size(Q,2), size(R,2));
meanRMSE_finalValues = zeros(size(Q,2), size(R,2));

% Computing Time for each filter and for each t of all the iterations
steps=((nMax-nMin)/nStep);
run_TIME_PF = zeros(nIterations,1);
avg_TIME_PF = zeros(steps,1);

% In this case is Q static, but needs to be adapted on the dynamic sze of R
for m=1:size(N,2)
    for i=1:nIterations
        [x_state_pf(:, :, i), x_t_vec, P_cov(:,:,:,i),run_TIME_PF(i)] = Main_isotropic_PF_NM(plotting, Q, R, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i),P_ini,N(m));
        for l=1 : size(x_state, 2)
            RMSE(m, i) = norm((x_state(:,l, i)- x_t_vec'));
        end
    end
    meanRMSEvalues(m) = mean(mean(RMSE))
    meanRMSE_finalValues(m) = mean(mean(RMSE(:, 1500:end)));
    RMSE = zeros(nIterations, 1801);
    avg_TIME_PF(m) = mean(run_TIME_PF);
end
% meanRMSEvalues = mean(meanRMSEvalues')
% % Mean of overall simualation
% meanSimValues=mean(meanRMSEvalues');
% plot(N,avg_TIME_PF)

%%Ploting of N related results
mkdir('NumberOfParticles')
h = figure;
size(N)
size(meanRMSEvalues)
plot(N, meanRMSEvalues);
set(gca, 'XScale', 'log', 'YScale', 'log', 'ZLim', [0 2000], 'FontSize', 10)
colormap jet
xlabel('N (Number of Particle)')
ylabel('RMSE (m)')
title(sprintf('RMSE vs N\n R = %.2f and Q = %d',R ,Q))
shading interp
print('./NumberOfParticles/Mean_N','-dsvg')
% savefig(h, strcat('./Results/', func2str(filter_function), int2str(nIterations)))

figure('Visible','off')
plot(N, avg_TIME_PF);
set(gca, 'XScale', 'log', 'YScale', 'log', 'ZLim', [0 2000], 'FontSize', 10)
colormap jet
xlabel('N (Number of Particle)')
ylabel('Time Average (s)')
title(sprintf('Computational Time vs N\n R = %.2f and Q = %d',R ,Q))
shading interp
print('./NumberOfParticles/Time_N','-dsvg')
% savefig(h, strcat('./Results/', func2str(filter_function), int2str(nIterations)))


%%End Q Optimization
%==========================================================================
end

