
function [x_state,P_cov,K_HINF_gain]=AHINF_form(xy1,xy2,h_0,alpha,x_state_ini,P_cov_ini,F,G,Q,R)
 
%% Adaptive H-Infinity  Filter
%% =========================

persistent firstRun
persistent  Q_adpt R_adpt alpha_adpt eps_adpt 
persistent Po X_s P_s sigma L S

%% initilize
%initilize our initial, prior particle distribution as a gaussian around
%the true initial value
if isempty(firstRun)
     X_s = x_state_ini;
     P_s = P_cov_ini;    %not used
     L = eye(2);
     S = eye(2);
     sigma = 0;
     
     Q_adpt = Q;
     R_adpt=R;
     alpha_adpt = 0.9;
     eps_adpt = 0;    

    firstRun = 1;
end

%% Prediction
%% =========================
%% Equation 2
P_s = F*P_s*F' + G *Q*G'; % Error covariance extrapolation
%% =========================
%% Correction
%% =========================
%% nonlinear measurement eq
 
h = hk(xy1,xy2,X_s,h_0); 
%% Jacobian of nonlinear measurement eq.
H = JH(xy1,xy2,X_s,h_0); 
%% ===============================================
%% Equation 4: Innovation
v = (alpha-h);
S_bar = L'*S*L;

%Measurement Covariance Update
R_adpt = alpha_adpt*R_adpt+(1-alpha_adpt)*(eps_adpt*eps_adpt'+H*P_s*H');

%% Equation 5: Kalman gain

% k =P_s*H'*inv(S) ;
k = P_s*inv(eye(2)-sigma.*S_bar*P_s+H'*inv(R_adpt)*H*P_s)*H'*inv(R_adpt);
K_HINF_gain = k;
%% Equation 6: State update
X_s = X_s +k*(v);
x_state = X_s ;

%computing of residual
eps_adpt = alpha - hk(xy1,xy2,X_s,h_0); 

%% Equation 7: Error covariance update
P_s = P_s*inv(eye(2)-sigma.*S_bar*P_s+H'*inv(R_adpt)*H*P_s)+Q_adpt;

P_cov=P_s;

%% Equation 8: State covariance update
Q_adpt = alpha_adpt*Q_adpt+(1-alpha_adpt)*(K_HINF_gain*v*v'*K_HINF_gain');

end
%% ===============================================
%% h(X): Nonlinear measurement eq

function h=hk(uav_init_pos, uav_actual_pos,X_s,h_0)

uav_init_pos = [uav_init_pos, h_0];
uav_actual_pos = [uav_actual_pos, h_0];

X_predicted = [X_s; h_0];

h=norm(X_predicted - uav_init_pos')^2 / norm(X_predicted - uav_actual_pos')^2;
end

function H=JH(uav_init_pos, uav_actual_pos, X_s,h_0)

X_predicted = [X_s;h_0];

uav_init_pos = [uav_init_pos, h_0];
uav_actual_pos = [uav_actual_pos h_0];

x=X_predicted(1);
y=X_predicted(2);

a = norm(X_predicted - uav_init_pos')^2;
b = norm(X_predicted - uav_actual_pos')^2;
H=[(2 * (x - uav_init_pos(1)) * b - 2 * (x - uav_actual_pos(1)) * a) / b^2,...  
   (2 * (y - uav_init_pos(2)) * b - 2 * (y - uav_actual_pos(2)) * a) / b^2];
end