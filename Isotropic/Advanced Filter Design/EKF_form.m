
function [x_state,P_cov,K_EKF_gain]=KF_form(xy1,xy2,h_0,alpha,x_state_ini,P_cov_ini,F,G,Q,R)
 
%% Extended Kalman Filter
%% =========================

persistent firstRun
persistent Po X_s P_s

%% initilize
%initilize our initial, prior particle distribution as a gaussian around
%the true initial value
if isempty(firstRun)
     X_s = x_state_ini;
     P_s = P_cov_ini;    %not used
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

S = (H*P_s*H' +R);
%% Equation 5: Kalman gain

% k =P_s*H'*inv(S) ;
k = (P_s * H') / (S);

K_EKF_gain=k;
%% Equation 6: State update
X_s = X_s +k*(v);
x_state = X_s ;

%% Equation 7: Error covariance update
P_s = (eye(size(k*H))-k*H)*P_s;
P_cov=P_s;
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