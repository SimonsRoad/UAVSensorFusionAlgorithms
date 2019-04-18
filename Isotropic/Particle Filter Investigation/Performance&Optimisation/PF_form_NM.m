
function [x_state,P_cov,K_EKF_gain]=PF_form_NM(xy1,xy2,h_0,alpha,x_state_ini,P_cov_ini,F,G,Q,R,N)

    persistent firstRun
    persistent Po X_s P_s
    
    %% initilize
    %initilize our initial, prior particle distribution as a gaussian around
    %the true initial value
%     N = 1000; % Number of Particles
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
        Zhat = hk(xy1,xy2,Po_pr(:, i),h_0)+ sqrt(R) * randn;       % Measurment Value of Particle
        diff =  Z - Zhat;                                               % Distance to observation
%         R or 2*pi
        w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(diff)^2 / 2 / R);   % Finding the weight
    end

    % Normalize to get the sampling weight
    wsum = sum(w);
    
    for i = 1 : N
        w(i) = w(i) / wsum;
    end
    %=========================
    
    % Resample -  Multinomial
    for i = 1 : N
        u = rand; % uniform random number between 0 and 1
        qtempsum = 0;
        for j = 1 : N
            qtempsum = qtempsum + w(j);
            if qtempsum >= u
                Po(:, i) = Po_pr(:, j);
                break;
            end
        end
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
% function h=hk(xy1, xy2, X_s,x_state_ini,h_0)
%  
% six=xy1(1);
% siy=xy1(2);
%  
% skx=xy2(1);
% sky=xy2(2);
%  
% xk=X_s(1);
% yk=X_s(2);
% 
% h=((sqrt((((xk-six)^2)+((yk-siy)^2)+(h_0^2)))^2)/(sqrt((((xk-skx)^2)+((yk-sky)^2)+(h_0^2)))^2));
% end

function h=hk(uav_init_pos, uav_actual_pos,X_s,h_0)

uav_init_pos = [uav_init_pos, h_0];
uav_actual_pos = [uav_actual_pos, h_0];

X_predicted = [X_s; h_0];

h=norm(X_predicted - uav_init_pos')^2 / norm(X_predicted - uav_actual_pos')^2;
end


