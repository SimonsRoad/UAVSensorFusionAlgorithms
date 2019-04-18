function plot_min_circles(x_vec_H,x_vec_dot,r_est_l,r_est_h,peak_intvl_width)
%   This function aims at plotting the circles at the determined peak power
%   However only portion of those circles are plotted as the width of the
%   determined area has been found by the confidence interval
%   This function tries to plot those portions of circle whose center is
%   the peak power position at k=k_H_m and that are delimited by the
%   confidence interval

%   Numerous angles are defined:
%       -   beta: angle between horizontal and flyby path
%       -   gamma: angle between flyby and line joining center and
%       intersection of radius with interval bound
%       -   delta: angle between radiis joining the two interval bounds

hold on

global plot_scaling

beta=atan2(x_vec_dot(:,2),x_vec_dot(:,1));  %   Radians

%   High radius circle: diameter wider than interval width
if ((peak_intvl_width/(2*r_est_h))<=1)
  gamma_h=acos((peak_intvl_width/(2*r_est_h)));
  begin_angle_h=beta+gamma_h;
  end_angle_h=pi+beta-gamma_h;
  theta_data_min_h_1=linspace(begin_angle_h,end_angle_h,50);
  theta_data_min_h_2=linspace(begin_angle_h+pi,end_angle_h+pi,50);
  plot((x_vec_H(:,1)+r_est_h(:,1).*cos(theta_data_min_h_1))/plot_scaling,(x_vec_H(:,2)+r_est_h(:,1).*sin(theta_data_min_h_1))/plot_scaling,'--m','Linewidth',2);
  plot((x_vec_H(:,1)+r_est_h(:,1).*cos(theta_data_min_h_2))/plot_scaling,(x_vec_H(:,2)+r_est_h(:,1).*sin(theta_data_min_h_2))/plot_scaling,'--m','Linewidth',2);
end
%   High radius circle:  diameter narrower than interval width then plot
%   the whole circle
if ((peak_intvl_width/(2*r_est_h))>1)
  theta_data_min_h=linspace(0,2*pi,200);
  plot((x_vec_H(:,1)+r_est_h(:,1).*cos(theta_data_min_h))/plot_scaling,(x_vec_H(:,2)+r_est_h(:,1).*sin(theta_data_min_h))/plot_scaling,'--m','Linewidth',2);
end



%   Low radius circle: diameter wider than interval width
if ((peak_intvl_width/(2*r_est_l))<=1)
  gamma_l=acos((peak_intvl_width/(2*r_est_l)));
  begin_angle_l=beta+gamma_l;
  end_angle_l=pi+beta-gamma_l;
  theta_data_min_l_1=linspace(begin_angle_l,end_angle_l,50);
  theta_data_min_l_2=linspace(begin_angle_l+pi,end_angle_l+pi,50);
  plot((x_vec_H(:,1)+r_est_l(:,1).*cos(theta_data_min_l_1))/plot_scaling,(x_vec_H(:,2)+r_est_l(:,1).*sin(theta_data_min_l_1))/plot_scaling,'--m','Linewidth',2);
  plot((x_vec_H(:,1)+r_est_l(:,1).*cos(theta_data_min_l_2))/plot_scaling,(x_vec_H(:,2)+r_est_l(:,1).*sin(theta_data_min_l_2))/plot_scaling,'--m','Linewidth',2);
end
%   Low radius circle:  diameter narrower than interval width then plot
%   the whole circle
if ((peak_intvl_width/(2*r_est_l))>1)
  theta_data_min_l=linspace(0,2*pi,200);
  plot((x_vec_H(:,1)+r_est_l(:,1).*cos(theta_data_min_l))/plot_scaling,(x_vec_H(:,2)+r_est_l(:,1).*sin(theta_data_min_l))/plot_scaling,'--m','Linewidth',2);
end


