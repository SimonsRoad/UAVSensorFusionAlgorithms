
function [x_state,P_cov,K_EKF_gain]=EPF_form(xy1,xy2,h_0,P_r_filt_ratio ,x_state_ini,P_cov_ini,F_KF,G_KF,Q_KF,R_KF,N)

    persistent firstRun
    persistent n m kappa
    persistent Po P_s
    
    % initialize the variables
    alpha=P_r_filt_ratio;   % receiver power ratio
    F=F_KF;                 % phai in equations
    G=G_KF;                 % Noise matrix: (unity for pure additive gaussian noise)
    Q=Q_KF;                 % Process noise matrix: better to be small std for position and power
    R=R_KF;                 % noise on alpha: enabled if wanted 
    x_N = Q_KF;             % Noise covariance in the system (i.e. process noise in the state update, here, we'll use a gaussian.)
    x_R = R;                % Noise covariance in the measurement (i.e. the Quail creates complex illusions in its trail!)
     

       %% initilize
    %initilize our initial, prior particle distribution as a gaussian around
    %the true initial value
    if isempty(firstRun)
        
        X_s = x_state_ini;
           
        %========================================
        %InitializeParticle Filter
        for i = 1 : N     
            % This is the region in which the target initialy is assumed
            % The target is randomly placed in this area (information from code extracted)
            P_init(1,:)=4000*rand(1,N) + 4000; 
            P_init(2,:)=4000*rand(1,N) + 4000;
            
            % Initial Covariance for all particles
            P_s_init(:,:,i) = P_cov_ini;
        end
        
        %Particle Array
        Po = P_init;  % Particle starts off
        P_s = P_s_init;
        %========================================
        %Design Paramters UKF
        n = 2;
        m = 1;
        kappa = 1;
        
        firstRun = 1;
    end
    
    %%========================================
    %% Particle Filter Loop
    % Current Target Measurment
    Z = alpha;
    
    % Baseline estimation
    for i = 1 : N

%         Po(:, i) = Po(:, i) + G*sqrt(Q) * [randn; randn];    
        
%         % Error covariance extrapolation
%         P_s(:,:,i) = F*(P_s(:,:,i))*F' + Q;
        
        %computing the sigma points
        [Xi W] = SigmaPoints(Po(:, i), P_s(:,:,i), kappa);
        
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
          hXi(:, k) = hx(fXi(:,k),xy1,xy2,h_0); 
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
        % Distance to observation
        diff =  Z - zp;    
        
        % to predict a new landmark position
        Po_pr(:, i) = Po(:, i) + K_EKF_gain*(diff);

        %% Precidtinc merged Covariance %%
        P_s(:,:,i) = Pp - K_EKF_gain*Pz*K_EKF_gain'; 
        
        w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(diff)^2 / 2 / R);   % Finding the weight
    end

    % Normalize to get the sampling weight
    wsum = sum(w);
    
    for i = 1 : N
        
        w(i) = w(i) / wsum;
    end
    %=========================
    
    % Resample - Residual Systematic Resampling Method
    
    M = length(w);
    i =1;
    u = rand/N;
    j = 0;
    while j < M
        j = j + 1;
        Ns = floor(N*(w(j)-u))+1;
          counter = 1;
          while counter <= Ns
            Po(:,i)=Po_pr(:,j); 
            i = i + 1; counter = counter + 1;
          end
        u = u + Ns/N-w(j);
    end 
    %=====
%     hold on
%     p2= plot(Po(1,:)/10^3,Po(2,:)/10^3,'r.');
%     pause(0.001)
%     delete(p2)
    % Mean of the System
    Pest = mean(Po');


    %% Equation 6: State update
    
    X_s = Pest;
    %X_s=MovAvgFilter(X_s');
    x_state = X_s;

    %% Equation 7: Error covariance update (not updated bc only internal use)
    P_cov=P_cov_ini;

end 

%% ===============================================
function h=hx(X_s,uav_init_pos, uav_actual_pos,h_0)

uav_init_pos = [uav_init_pos, h_0];
uav_actual_pos = [uav_actual_pos, h_0];

X_predicted = [X_s; h_0];

h=norm(X_predicted - uav_init_pos')^2 / norm(X_predicted - uav_actual_pos')^2;
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
alpha =0.001;
beta =2;
% Recalculate kappa according to scaling parameters
lambda = alpha^2*(n+kappa)-n;

% Allocate space

wPts=zeros(1,nPts);
xPts=zeros(n,nPts);

% Calculate matrix square root of weighted covariance matrix
Psqrtm=(chol(((n+lambda)*P)))';  

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


function avg = MovAvgFilter(x)

    persistent prevAvg n xbuf
    persistent firstRun


    if isempty(firstRun)
      n    = 100;
      
      xbuf = ones(size(x,1),n);
      
      for i=1:n
          xbuf(:,i)=x;
      end
      
      prevAvg = x;

      firstRun = 1;  
    end


    for m=1:(n-1)
      xbuf(:,m) = xbuf(:,(m+1));
    end
    
    xbuf(:,n) = x;

    avg = prevAvg + (x - xbuf(:,1)) / n;


    prevAvg = avg;
end