function [x_state,P_cov,K_EKF_gain]=UKF_form(x_init,x_current,h_0,alpha,X_s,P_s,F_KF,G_KF,Q_KF,R_KF,G_t_1,G_t_2)

% Parameter:
%   x_vec_all(1,:)              inital uav position vector
%   x_vec_all(k,:)              current uav position vector
%   h_0                         Constant altitude of the UAV[m]
%   P_r_filt_ratio(k,1)         Alpha: Power ratio between initial and current : see 'alpha' in report
%   x_state_ini                 Initial state guess - Middle of the area is the first guess
%   P_cov_ini                   Initial state covariance guess - Change if needed
%   F_KF                        Dynamics matrix: unity because model is static
%   G_KF                        Noise matrix: unity for pure additive gaussian noise
%   Q_KF                        Process noise matrix: better to be small std for position and power
%   R_KF                        Specify noise on alpha: enable if wanted

    persistent n m kappa
    persistent firstRun

    %**********************************************************************
    %Init
    if isempty(firstRun)
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
% This function returns the scaled symmetric sigma point distribution.
%
%  [xPts, wPts, nPts] = scaledSymmetricSigmaPoints(x,P,alpha,beta,kappa)  
%
% Inputs:
%	 x	      mean
%	 P	      covariance
%        alpha        scaling parameter 1
%        beta         extra weight on zero'th point
%	 kappa	      scaling parameter 2 (usually set to default 0)
%
% Outputs:
%        xPts	 The sigma points
%        wPts	 The weights on the points
%	 nPts	 The number of points
%

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
P = nearestSPD(P);
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
%%
function Ahat = nearestSPD(A)
    % nearestSPD - the nearest (in Frobenius norm) Symmetric Positive Definite matrix to A
    % usage: Ahat = nearestSPD(A)
    
    % arguments: (input)
    %  A - square matrix, which will be converted to the nearest Symmetric
    %    Positive Definite Matrix.
    %
    % Arguments: (output)
    %  Ahat - The matrix chosen as the nearest SPD matrix to A.

    if nargin ~= 1
      error('Exactly one argument must be provided.')
    end

    % test for a square matrix A
    [r,c] = size(A);
    if r ~= c
      error('A must be a square matrix.')
    elseif (r == 1) && (A <= 0)
      % A was scalar and non-positive, so just return eps
      Ahat = eps;
      return
    end

    % symmetrize A into B
    B = (A + A')/2;

    % Compute the symmetric polar factor of B. Call it H.
    % Clearly H is itself SPD.
    [U,Sigma,V] = svd(B);
    H = V*Sigma*V';

    % get Ahat in the above formula
    Ahat = (B+H)/2;

    % ensure symmetry
    Ahat = (Ahat + Ahat')/2;

    % test that Ahat is in fact PD. if it is not so, then tweak it just a bit.
    p = 1;
    k = 0;
    while p ~= 0
      [R,p] = chol(Ahat);
      k = k + 1;
      if p ~= 0
        % Ahat failed the chol test. It must have been just a hair off,
        % due to floating point trash, so it is simplest now just to
        % tweak by adding a tiny multiple of an identity matrix.
        mineig = min(eig(Ahat));
        Ahat = Ahat + (-mineig*k.^2 + eps(mineig))*eye(size(A));
      end
    end
end
%%
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