
function [x_state,P_cov,K_EKF_gain]=PF_form_SYS(xy1,xy2,h_0,alpha,x_state_ini,P_cov_ini,F,G,Q,R)

    persistent firstRun
    persistent Po X_s P_s
    
    %% initilize
    N = 1000; % Number of Particles
    K = 100;  % Gitter Factor
    if isempty(firstRun)
        
        X_s = x_state_ini;
        P_s = P_cov_ini;    %not used
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
        Po_pr(:, i) = Po(:, i) + K*G*sqrt(Q) * [randn; randn];
        Zhat = hk(xy1,xy2,Po_pr(:, i),h_0)+ sqrt(R) * randn;            % Measurment Value of Particle
        diff =  Z - Zhat;                                               % Distance to observation
%         R or 2*pi
        w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(diff)^2 / 2 / R);   % Finding the weight
    end

    % normalise weights & form cumulative distribution 
    w=cumsum(w/sum(w));      

    %%=========================
    %% Resample - Systematic Method
    add=1/N; 
    stt=add*rand(1);     
    selection=[stt : add : stt+(N-1)*add];
    j=1;    
    
    x_post=zeros(size(Po_pr));
    for i=1:N  
        while selection(i) >= w(j);        
            j=j+1;      
        end     
        x_post(:,i)=Po_pr(:,j);              
    end
    Po = x_post;

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
    x_state = X_s;

    %% Equation 7: Error covariance update
    P_s = eye(2);
    K_EKF_gain = [0; 0];
    P_cov=P_s;

end 

%% ===============================================
%% h(X): Nonlinear measurement equation
function h=hk(uav_init_pos, uav_actual_pos,X_s,h_0)

uav_init_pos = [uav_init_pos, h_0];
uav_actual_pos = [uav_actual_pos, h_0];

X_predicted = [X_s; h_0];

h=norm(X_predicted - uav_init_pos')^2 / norm(X_predicted - uav_actual_pos')^2;
end


