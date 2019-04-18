
function [x_state,P_cov,K_EKF_gain]=AKF_form(xy1,xy2,h_0,alpha,X_s,P_s,F,G,Q,R)
%% Extended Kalman Filter
%% =========================
persistent  Q_adpt R_adpt alpha_adpt eps_adpt 
persistent firstRun
 
 if isempty(firstRun)
     Q_adpt = Q;
     R_adpt=R;
     alpha_adpt = 0.9;
     eps_adpt = 0;
 end

%% initilize
%initilize our initial, prior particle distribution as a gaussian around
%the true initial value

%% Prediction
%% =========================
%% Equation 2
P_s = F*P_s*F' + G *Q_adpt*G'; % Error covariance extrapolation
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
%Measurement Covariance Update
R_adpt = alpha_adpt*R_adpt+(1-alpha_adpt)*(eps_adpt*eps_adpt'+H*P_s*H');
S = (H*P_s*H' +R_adpt);
%% Equation 5: Kalman gain

% k =P_s*H'*inv(S) ;
k = (P_s * H') / (S);

K_EKF_gain=k;
%% Equation 6: State update
X_s = X_s +k*(v);
x_state = X_s ;

%computing of residual
eps_adpt = alpha - hk(xy1,xy2,X_s,h_0); 

%% Equation 7: Error covariance update
P_s = (eye(size(k*H))-k*H)*P_s;
P_cov=P_s;

%% Equation 8: State covariance update
Q_adpt = alpha_adpt*Q_adpt+(1-alpha_adpt)*(K_EKF_gain*v*v'*K_EKF_gain');

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
