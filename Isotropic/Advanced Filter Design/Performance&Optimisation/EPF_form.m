
function [x_state,P_cov,K_EKF_gain]=EPF_form(xy1,xy2,h_0,P_r_filt_ratio,x_state_ini,P_cov_ini,F_KF,G_KF,Q_KF,R_KF,M)

    persistent firstRun
    persistent X_s P_s x_hat N
    persistent Po
    

    % initialize the variables
    alpha=P_r_filt_ratio;   % receiver power ratio
    F=F_KF;                 % phai in equations
    G=G_KF;                 % Noise matrix: (unity for pure additive gaussian noise)
    Q=Q_KF;                 % Process noise matrix: better to be small std for position and power
    R=R_KF;                 % noise on alpha: enabled if wanted 
    x_N = Q_KF;             % Noise covariance in the system (i.e. process noise in the state update, here, we'll use a gaussian.)
    x_R = R;                % Noise covariance in the measurement (i.e. the Quail creates complex illusions in its trail!)
    N = 500;                 % The number of particles the system generates. The larger this is, the better your approximation, but the more computation you need.
    
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
        
        firstRun = 1;
    end
    
    %%========================================
    %% Particle Filter Loop
    % Current Target Measurment
    Z = alpha;
    
    % Baseline estimation
    for i = 1 : N
        
        %Po(:, i) = F*Po(:, i) + G*sqrt(Q) * [randn; randn];
        % Error covariance extrapolation
        P_s(:,:,i) = F*(P_s(:,:,i))*F' + Q; 
        % Measurment Value of Particle
        Zhat = hk(xy1,xy2,Po(:, i),x_state_ini,h_0)+ sqrt(R) * randn;
        % Distance to observation
        diff =  Z - Zhat;                                                      
  
        % Jacobian of nonlinear measurement
        H = JH(xy1,xy2,Po(:, i),h_0);
        % Calculate the Kalman gain
        S = (H*P_s(:,:,i)*H' +R);
        %Calculating Kalman Gain
        K_EKF_gain =P_s(:,:,i)*H'*inv(S) ;

        % to predict a new landmark position
        Po_pr(:, i) = Po(:, i) + K_EKF_gain*(diff);

        % Update the covariance of this landmark
        P_s(:,:,i) = (eye(2)-K_EKF_gain*H)*P_s(:,:,i);
        
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
%% h(X): Nonlinear measurement equation
function h=hk(xy1, xy2, X_s,x_state_ini,h_0)
 
six=xy1(1);
siy=xy1(2);
 
skx=xy2(1);
sky=xy2(2);
 
xk=X_s(1);
yk=X_s(2);

% h with h_0,G_t_2 and G_t_1
h=((sqrt((((xk-six)^2)+((yk-siy)^2)+(h_0^2)))^2)/(sqrt((((xk-skx)^2)+((yk-sky)^2)+(h_0^2)))^2));
end

function H=JH(xy1,xy2,X_s,h_0)
six=xy1(1);
siy=xy1(2);
 
skx=xy2(1);
sky=xy2(2);
 
xk=X_s(1);
yk=X_s(2);

% H without h_0
% H(1)= ((2*skx - 2*xk)*((six - xk)^2 + (siy - yk)^2))/((skx - xk)^2 + (sky - yk)^2)^2 - (2*six - 2*xk)/((skx - xk)^2 + (sky - yk)^2);
% H(2)= ((2*sky - 2*yk)*((six - xk)^2 + (siy - yk)^2))/((skx - xk)^2 + (sky - yk)^2)^2 - (2*siy - 2*yk)/((skx - xk)^2 + (sky - yk)^2);

% H with h_0 
H(1)= ((2*skx - 2*xk)*((six - xk)^2 + (siy - yk)^2 + h_0^2))/((skx - xk)^2 + (sky - yk)^2 + h_0^2)^2 - (2*six - 2*xk)/((skx - xk)^2 + (sky - yk)^2 + h_0^2) ;
H(2)= ((2*sky - 2*yk)*((six - xk)^2 + (siy - yk)^2 + h_0^2))/((skx - xk)^2 + (sky - yk)^2 + h_0^2)^2 - (2*siy - 2*yk)/((skx - xk)^2 + (sky - yk)^2 + h_0^2);
end
%Moving Average for smoothing 
function avg = MovAvgFilterX(x)
    %
    %
    persistent prevAvg n xbuf
    persistent firstRun


    if isempty(firstRun)
      n    = 100;
      xbuf = x*ones(n+1, 1);

      k = 1;
      prevAvg = x;

      firstRun = 1;  
    end


    for m=1:n
      xbuf(m) = xbuf(m+1);
    end
    xbuf(n+1) = x;

    avg = prevAvg + (x - xbuf(1)) / n;


    prevAvg = avg;
end

function avg = MovAvgFilterY(x)
    %
    %
    persistent prevAvg n xbuf
    persistent firstRun


    if isempty(firstRun)
      n    = 200;
      xbuf = x*ones(n+1, 1);

      k = 1;
      prevAvg = x;

      firstRun = 1;  
    end


    for m=1:n
      xbuf(m) = xbuf(m+1);
    end
    xbuf(n+1) = x;

    avg = prevAvg + (x - xbuf(1)) / n;


    prevAvg = avg;
end

function avg = MovAvgFilter(x)

    persistent prevAvg n xbuf
    persistent firstRun


    if isempty(firstRun)
      n    = 50;
      
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

