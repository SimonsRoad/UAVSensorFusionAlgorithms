function [x_state,P_cov,K_EKF_gain]=UKF_form(x_init,x_current,h_0,alpha,x_state_ini,P_cov_ini,F_KF,G_KF,Q_KF,R_KF,G_t_1,G_t_2)
%% Unscendet Kalman Filter
    persistent n m kappa
    persistent firstRun X_s P_s

    %**********************************************************************
    %Init
    if isempty(firstRun)
      X_s = x_state_ini;
      P_s = P_cov_ini;
      n = 2;
      m = 1;
      kappa = 0;
      firstRun = 1;
    end

    %computing the sigma points
    [Xi W] = SigmaPoints(X_s, P_s, kappa);

    % give the sigma points to the non-linear function and store the results
    % in seperate vectore
    fXi = zeros(n, 2*n+1);
    for k = 1:2*n+1
      fXi(:, k) = Xi(:,k);
    end

    % exectuing the unscented transformation for estimated values
    [xp Pp] = UT(fXi, W, Q_KF);

    % maps sigma points to the measurement space and store in seperate vectore
    hXi = zeros(m, 2*n+1);
    for k = 1:2*n+1

      hXi(:, k) = hx(fXi(:,k),x_init,x_current,h_0,G_t_1,G_t_2);
    end
    
    % exectuing the unscented transformation to measured values
    [zp Pz] = UT(hXi, W, R_KF);

    % Compute Cross Co-relation Matrix between state space and predicted space
    Pxz = zeros(n, m);
    for k = 1:2*n+1
      Pxz = Pxz + W(k)*(fXi(:,k) - xp)*(hXi(:,k) - zp)';
    end

    % Comput Kalam gain with Predicted Covariance Matrix
    K_EKF_gain = Pxz*inv(Pz);

    %% Update INS states (position only) %%
    X_s = xp + K_EKF_gain*(alpha - zp);
    x_state = X_s;

    %% Precidtinc merged Covariance %%
    P_s = Pp - K_EKF_gain*Pz*K_EKF_gain';
    P_cov = P_s;
end
%%
function h=hx(X_s,uav_init_pos, uav_actual_pos,h_0,G_t_1,G_t_2)

uav_init_pos = [uav_init_pos, h_0];
uav_actual_pos = [uav_actual_pos, h_0];

X_predicted = [X_s; h_0];

h=(G_t_2/G_t_1)*norm(X_predicted - uav_init_pos')^2 / norm(X_predicted - uav_actual_pos')^2;
end


%%
function [xPts wPts] = SigmaPoints(x, P, kappa)
% Number of sigma points and scaling terms
n    = size(x(:),1);
nPts = 2*n+1;            % we're using the symmetric SUT

%Design Parameters
alpha =0.000001;
beta =0.05;
% Recalculate kappa according to scaling parameters
lambda = alpha^2*(n+kappa)-n;

% Allocate space

wPts=zeros(1,nPts);
xPts=zeros(n,nPts);

% Calculate matrix square root of weighted covariance matrix
Psqrtm=(chol((n+lambda)*P))';  

% Array of the sigma points
xPts=[zeros(size(P,1),1) -Psqrtm Psqrtm];

% Add mean back in
xPts = xPts + repmat(x,1,nPts);  

% Array of the weights for each sigma point
wPts=[lambda 0.5*ones(1,nPts-1) 0]/(n+lambda);

% Now calculate the zero'th covariance term weight
wPts(nPts+1) = wPts(1) + (1-alpha^2) + beta;
end

function [xm xcov] = UT(Xi, W, noiseCov)  
%
%
    [n, kmax] = size(Xi);

    xm = 0;
    for k=1:kmax
      xm = xm + W(k)*Xi(:, k);
    end

    xcov = zeros(n, n);
    for k=1:kmax
      xcov = xcov + W(k)*(Xi(:, k) - xm)*(Xi(:, k) - xm)';
    end
    xcov = xcov + noiseCov;
end