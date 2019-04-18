
function [x_state,P_cov,K_EKF_gain]=EPF_form(xy1,xy2,h_0,P_r_filt_ratio,x_state_ini,P_cov_ini,F_KF,G_KF,Q_KF,R_KF,G_t_1,G_t_2)

    persistent firstRun
    persistent X_s P_s x_hat N
    persistent  Q_adpt R_adpt alpha_adpt eps_adpt 
    persistent n m kappa
    persistent Po
    

    % initialize the variables
    alpha=P_r_filt_ratio;   % receiver power ratio
    F=F_KF;                 % phai in equations
    G=G_KF;                 % Noise matrix: (unity for pure additive gaussian noise)
    Q=Q_KF;                 % Process noise matrix: better to be small std for position and power
    R=R_KF;                 % noise on alpha: enabled if wanted 
    N = 100;                 % The number of particles the system generates. The larger this is, the better your approximation, but the more computation you need.
    
    %% initilize
    %initilize our initial, prior particle distribution as a gaussian around
    %the true initial value
    if isempty(firstRun)
        
        X_s = x_state_ini;
        Q_adpt = Q_KF;
        R_adpt=R_KF;
        alpha_adpt = 1;
        eps_adpt = 0;
           
        %========================================
        %InitializeParticle Filter
        for i = 1 : N     
            % This is the region in which the target initialy is assumed
            % The target is randomly placed in this area (information from code extracted)
            P_init(1,:)=6000*rand(1,N) + 3000; 
            P_init(2,:)=6000*rand(1,N) + 3000;
            
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
        
%         % Measurment Value of Particle
%         Zhat = hk(xy1,xy2,Po(:, i),x_state_ini,h_0)+ sqrt(R) * randn;    
        
        % Error covariance extrapolation
        P_s(:,:,i) = F*(P_s(:,:,i))*F' + G*sqrt(Q) * [randn; randn];
        
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
          hXi(:, k) = hx(x_state_ini,fXi(:,k),xy1,xy2,h_0,G_t_1,G_t_2); 
        end
        
        %Measurement Covariance Update
        R_adpt = alpha_adpt*R_adpt+(1-alpha_adpt)*(eps_adpt*eps_adpt'+H*P_s(:,:,i)*H');
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
    
   %==== Ploting of the Particle Postion in each time step -  Deactivated
    %(Uncomment for special effects :))
%     hold on
%     p2= plot(Po(1,:)/10^3,Po(2,:)/10^3,'r.');
%     pause(0.001)
%     delete(p2)
    % Mean of the System
    Pest = mean(Po');


    %% Equation 6: State update
    
    X_s = Pest;
    X_s=MovAvgFilter(X_s');
    x_state = X_s;

    %% Equation 7: Error covariance update (not updated bc only internal use)
    P_cov=P_cov_ini;

end 

%% ===============================================
%%
function h=hx(x_state_ini,X_s,xy1,xy2,h_0,G_t_1,G_t_2)
six=xy1(1);
siy=xy1(2);
 
skx=xy2(1);
sky=xy2(2);
 
xk=X_s(1);
yk=X_s(2);

% h (1X1)=alpha
h=(G_t_2/G_t_1)*((sqrt((((xk-six)^2)+((yk-siy)^2)+(h_0^2)))^2)/(sqrt((((xk-skx)^2)+((yk-sky)^2)+(h_0^2)))^2));

end
%%
function [xPts wPts] = SigmaPoints(x, P, kappa)
% Number of sigma points and scaling terms
n    = size(x(:),1);
nPts = 2*n+1;            % we're using the symmetric SUT

%Design Parameters
alpha =0.000001;
beta =0.5;
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
