
function [x_state,P_cov,K_EKF_gain]=PF_form(xy1,xy2,h_0,alpha,X_s,P_s,F,G,Q,R,N,G_t_1,G_t_2)

    persistent firstRun
    persistent Po
    
    %% initilize
    %initilize our initial, prior particle distribution as a gaussian around
    %the true initial value
    %N = 1000;
    if isempty(firstRun)
        
        firstRun = 1;
        
        %========================================
        %InitializeParticle Filter
        for i = 1 : N
            
            % This is the region in which the target initialy is assumed
            % The target is randomly placed in this area (information from code extracted)
            P_init(1,:)=4000*rand(1,N) + 4000; 
            P_init(2,:)=4000*rand(1,N) + 4000;            
        end
        
        %Particle Array
        Po = P_init;  % Particle starts off
        %========================================
    end
    
    %%========================================
    %% Particle Filter Loop
    % Current Target Measurment
    Z = alpha;
    
    % Baseline estimation
    for i = 1 : N
        Po_pr(:, i) = Po(:, i) + G*sqrt(Q) * [randn; randn];
        Zhat = hk(xy1,xy2,Po_pr(:, i),h_0,G_t_1,G_t_2)+ sqrt(R) * randn;    % Measurment Value of Particle
        diff =  Z - Zhat;                                                   % Distance to observation
        w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(diff)^2 / 2 / R);       % Finding the weight
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
    x_state = X_s;

    %% Equation 7: Error covariance update
    P_s = eye(2);
    K_EKF_gain = [0; 0];
    P_cov=P_s;

end 

%% ===============================================
%% h(X): Nonlinear measurement equation
function h=hk(uav_init_pos, uav_actual_pos,X_s,h_0,G_t_1,G_t_2)

uav_init_pos = [uav_init_pos, h_0];
uav_actual_pos = [uav_actual_pos, h_0];

X_predicted = [X_s; h_0];

h=(G_t_2/G_t_1)*norm(X_predicted - uav_init_pos')^2 / norm(X_predicted - uav_actual_pos')^2;
end

