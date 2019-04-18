function plot_animation_search(N_plots,C_loop,x_t_vec,x_vec,psi,r_est_l,r_est_h,centre_geo_circle,radius_geo_circle,x_state,k_obs,N_loops_fb,P_cov,p_e,r_d,psi_jammer)
%   This function aims at updating the plot of the search area by
%   replotting the uav, the jammer and geolocation data (goemetry, filter)
%   Numerous different plotting functions are found below



global plot_scaling x_bnd y_bnd fig_offset

%   All plot handles
persistent h_uav_plot h_jammer_plot h_uav_trace_plot h_uav_jammer_vec_plot h_r_est_plot h_geo_est_plot h_state_pos_est h_loiter_plot h_P_cov_ell

%   ------------------------    At first step         --------------------
if C_loop==1
    figure(N_plots);                                                        %   Create new figure
    %   Plot boundaries, background and scale axes once
    
        %   Area parameters
        bnd_linwid=1.5;                                                     %   Boundary linewidth
    
        bnd_mtx=[0 0; 0 y_bnd; x_bnd y_bnd; x_bnd 0; 0 0]/plot_scaling;   	%   Boundary point matrix
        plot(bnd_mtx(:,1),bnd_mtx(:,2),'-k','LineWidth',bnd_linwid);        %   Plot area boundaries
        grid on
        hold on
        axis([-fig_offset*x_bnd/plot_scaling x_bnd*(1+fig_offset)/plot_scaling -fig_offset*y_bnd/plot_scaling y_bnd*(1+fig_offset)/plot_scaling])     %   scaling: axis([xmin xmax ymin ymax])
        xlabel('x-position [km]','fontsize',12,'color','k');
        ylabel('y-position [km]','fontsize',12,'color','k');
        axis square
        
        h_jammer_plot=plot_jammer(x_t_vec,psi_jammer,[]);                                                  %   Plot jammer, handle is empty: [] at first iteration
        h_uav_trace_plot=plot_uav_trace(x_vec,[]);                                              %   Plot uav trace, handle is empty: [] at first iteration
        h_uav_plot=plot_uav(x_vec(C_loop,:),psi,[]);                                            %   Plot uav, handle is empty: [] at first iteration
        h_uav_jammer_vec_plot=plot_uav_jammer_vec(x_t_vec,x_vec(C_loop,:),[]);                  %   Plot uav-jammer vector, handle is empty: [] at first iteration
        h_r_est_plot=plot_r_est(x_vec(C_loop,:),r_est_l(:,1),r_est_h(:,1),[]);                  %   Plot range circles, handle is empty: [] at first iteration
        h_geo_est_plot=plot_geo_est(C_loop,centre_geo_circle(:,:),radius_geo_circle(:,1),[]);   %   Plot geo circles, handle is empty: [] at first iteration
        h_state_pos_est=plot_state_pos_est(x_state(:,end),[]);                                  %   Plot the position estimation via EKF cross-hair        
        h_P_cov_ell=plot_P_cov_ell(x_state(:,end),P_cov,p_e,[]);                                %   Plot the error covariance
        h_loiter_plot=plot_loiter(x_state(:,end),r_d,[]);                                       %   Plot loiter objective circle
        
 %   -------------------    At all other steps         -------------------- 
else
     plot_jammer(x_t_vec,psi_jammer,h_jammer_plot); 
     plot_uav(x_vec(C_loop,:),psi,h_uav_plot);
     plot_uav_trace(x_vec,h_uav_trace_plot);
     plot_uav_jammer_vec(x_t_vec,x_vec(C_loop,:),h_uav_jammer_vec_plot);
     plot_r_est(x_vec(C_loop,:),r_est_l(:,1),r_est_h(:,1),h_r_est_plot);
     plot_geo_est(C_loop,centre_geo_circle(:,:),radius_geo_circle(:,1),h_geo_est_plot);
     plot_state_pos_est(x_state(:,end),h_state_pos_est);
     plot_state_trace(x_state,C_loop,k_obs,N_loops_fb);
% % %      plot_P_cov_ell(x_state(:,end),P_cov,p_e,h_P_cov_ell);
     plot_loiter(x_state(:,end),r_d,h_loiter_plot);
     
end

%   Halt execustion for n seconds for good animation running time
n_pause=0.0001;
pause(n_pause);

end     %   End function









function [ h_uav_plot  ] = plot_uav(x_vec,psi,h_uav_plot)
%   This function plots the UAV in the desired simple shape for
%   visualisation. The function should be called with graphics window open

global x_bnd y_bnd plot_scaling

x=x_vec(1,1);
y=x_vec(1,2);
%   Geometric constants for the plot
length_ratio_x=1/25;                                                        %   Length ratio related to x_bnd
height_ratio_y=1/15;                                                        %   Height ratio related to y_bnd for the vertical components
height_ratio_y_2=height_ratio_y/3;                                          %   Height ratio related to y_bnd for the diagonal components

uav_linewidth=3;
uav_linstyle='-r';

%   UAV points definition for drawing with respect to origin
uav_p1=[(-x_bnd*length_ratio_x)/2 y_bnd*height_ratio_y_2/2]';
uav_p2=[0 (y_bnd*(height_ratio_y+height_ratio_y_2)/2)]';
uav_p3=[(x_bnd*length_ratio_x)/2 y_bnd*height_ratio_y_2/2]';
uav_p4=[(x_bnd*length_ratio_x)/2 -(y_bnd*(height_ratio_y+height_ratio_y_2)/2)]';
uav_p5=[0 -(y_bnd*(height_ratio_y+height_ratio_y_2)/2)+y_bnd*height_ratio_y_2]';
uav_p6=[(-x_bnd*length_ratio_x)/2 -(y_bnd*(height_ratio_y+height_ratio_y_2)/2)]';

%   Rotation matrix
psi_dec=psi-(pi/2);                                                         %   Shift psi as the UAV is drawn vertical by default
M_rot=[cos(psi_dec) -sin(psi_dec); sin(psi_dec) cos(psi_dec)];              %   Rotation matrix to rotate drawing as heading varies
%   Rotate positions
uav_p1=(M_rot*uav_p1)';
uav_p2=(M_rot*uav_p2)';
uav_p3=(M_rot*uav_p3)';
uav_p4=(M_rot*uav_p4)';
uav_p5=(M_rot*uav_p5)';
uav_p6=(M_rot*uav_p6)';

%   uav crosshair points definition for drawing 
uav_crosshair_p1=[-(x_bnd*length_ratio_x)/3.5 0];
uav_crosshair_p2=[(x_bnd*length_ratio_x)/3.5 0];
uav_crosshair_p3=[0 -(y_bnd*height_ratio_y)/3.5];
uav_crosshair_p4=[0 (y_bnd*height_ratio_y)/3.5];

%   Translate positions
x_vec=[x y];
uav_p1=uav_p1+x_vec;
uav_p2=uav_p2+x_vec;
uav_p3=uav_p3+x_vec;
uav_p4=uav_p4+x_vec;
uav_p5=uav_p5+x_vec;
uav_p6=uav_p6+x_vec;

uav_crosshair_p1=uav_crosshair_p1+x_vec;
uav_crosshair_p2=uav_crosshair_p2+x_vec;
uav_crosshair_p3=uav_crosshair_p3+x_vec;
uav_crosshair_p4=uav_crosshair_p4+x_vec;

%   Link points
uav_drawing_mtx=[uav_p1; uav_p2; uav_p3; uav_p4; uav_p5; uav_p6; uav_p1]/plot_scaling;
uav_crosshair_drawing_mtx=[uav_crosshair_p1; uav_crosshair_p2; uav_crosshair_p3; uav_crosshair_p4]/plot_scaling;

hold on
if isempty(h_uav_plot)
    h_uav_plot=plot(uav_drawing_mtx(:,1),uav_drawing_mtx(:,2),uav_linstyle,'Linewidth',uav_linewidth,'EraseMode','normal');
    h_uav_crosshair_1=plot(uav_crosshair_drawing_mtx(1:2,1),uav_crosshair_drawing_mtx(1:2,2),'-k','LineWidth',uav_linewidth/2,'EraseMode','normal');
    h_uav_crosshair_2=plot(uav_crosshair_drawing_mtx(3:4,1),uav_crosshair_drawing_mtx(3:4,2),'-k','LineWidth',uav_linewidth/2,'EraseMode','normal');
    h_uav_plot=[h_uav_plot h_uav_crosshair_1 h_uav_crosshair_2];
else
    set(h_uav_plot(1),'XData',uav_drawing_mtx(:,1),'YData',uav_drawing_mtx(:,2));
    set(h_uav_plot(2),'XData',uav_crosshair_drawing_mtx(1:2,1),'YData',uav_crosshair_drawing_mtx(1:2,2));
    set(h_uav_plot(3),'XData',uav_crosshair_drawing_mtx(3:4,1),'YData',uav_crosshair_drawing_mtx(3:4,2));   
end


end             %   End function





function [ h_uav_trace_plot  ] = plot_uav_trace(x_vec_all,h_uav_trace_plot)
%   This function plots the trace of the UAV for
%   visualisation. The function should be called with graphics window open

global plot_scaling

hold on
if isempty(h_uav_trace_plot)
    h_uav_trace_plot=plot(x_vec_all(:,1)/plot_scaling,x_vec_all(:,2)/plot_scaling,'-g','Linewidth',3,'EraseMode','normal');
else
    set(h_uav_trace_plot,'XData',x_vec_all(:,1)/plot_scaling,'YData',x_vec_all(:,2)/plot_scaling);   
end
end             %   End function



function [ h_jammer_plot  ] = plot_jammer(x_t_vec,psi_jammer,h_jammer_plot)
%   This function plots the jammer in the desired simple shape for
%   visualisation. The function should be called with graphics window open

%   It returns the handle of the jammer plot and takes into argument the
%   position of the jammer from place_jammer function and the current
%   jammer graphics handle

global x_bnd y_bnd plot_scaling

x_t=x_t_vec(1,1);
y_t=x_t_vec(1,2);
%   Geometric constants for the plot
length_ratio_x=1/20;    %   Length ratio related to x_bnd
height_ratio_y=1/20;    %   Height ratio related to y_bnd for the vertical components
height_ratio_y_2=height_ratio_y/2;  %   Height ratio related to y_bnd for the diagonal components
length_ratio_x_2=length_ratio_x/3;  %   Length ratio related to y_bnd for the diagonal components

jammer_linewidth=3;
jammer_linstyle='-k';

%   Jammer points definition for drawing 
jammer_p1=[(-(x_bnd*length_ratio_x)/2) (-(y_bnd*height_ratio_y)/2)];
jammer_p2=[(-(x_bnd*length_ratio_x)/2) ((y_bnd*height_ratio_y)/2)];
jammer_p3=[(-(x_bnd*length_ratio_x_2)/2) ((y_bnd*height_ratio_y/2)+y_bnd*height_ratio_y_2)];
jammer_p4=[((x_bnd*length_ratio_x_2)/2) ((y_bnd*height_ratio_y/2)+y_bnd*height_ratio_y_2)];
jammer_p5=[((x_bnd*length_ratio_x)/2) ((y_bnd*height_ratio_y)/2)];
jammer_p6=[((x_bnd*length_ratio_x)/2) (-(y_bnd*height_ratio_y)/2)];

%   Rotation matrix
psi_dec=psi_jammer-(pi/2);                                                	%   Shift psi as the UAV is drawn vertical by default
M_rot=[cos(psi_dec) -sin(psi_dec); sin(psi_dec) cos(psi_dec)];              %   Rotation matrix to rotate drawing as heading varies
%   Rotate positions
jammer_p1=(M_rot*(jammer_p1'))';
jammer_p2=(M_rot*(jammer_p2'))';
jammer_p3=(M_rot*(jammer_p3'))';
jammer_p4=(M_rot*(jammer_p4'))';
jammer_p5=(M_rot*(jammer_p5'))';
jammer_p6=(M_rot*(jammer_p6'))';

%   Translate positions
x_t_vec=[x_t y_t];
jammer_p1=jammer_p1+x_t_vec;
jammer_p2=jammer_p2+x_t_vec;
jammer_p3=jammer_p3+x_t_vec;
jammer_p4=jammer_p4+x_t_vec;
jammer_p5=jammer_p5+x_t_vec;
jammer_p6=jammer_p6+x_t_vec;

%   Jammer crosshair points definition for drawing 
jam_crosshair_p1=[x_t-(x_bnd*length_ratio_x)/3 y_t];
jam_crosshair_p2=[x_t+(x_bnd*length_ratio_x)/3 y_t];
jam_crosshair_p3=[x_t y_t-(y_bnd*height_ratio_y)/3];
jam_crosshair_p4=[x_t y_t+(y_bnd*height_ratio_y)/3];

%   Link points
jammer_drawing_mtx=[jammer_p1; jammer_p2; jammer_p3; jammer_p4; jammer_p5; jammer_p6; jammer_p1]/plot_scaling;
jammer_crosshair_drawing_mtx=[jam_crosshair_p1; jam_crosshair_p2; jam_crosshair_p3; jam_crosshair_p4]/plot_scaling;

hold on
if isempty(h_jammer_plot)
    h_jammer_plot=plot(jammer_drawing_mtx(:,1),jammer_drawing_mtx(:,2),jammer_linstyle,'Linewidth',jammer_linewidth,'EraseMode','normal');
    h_jammer_crosshair_1=plot(jammer_crosshair_drawing_mtx(1:2,1),jammer_crosshair_drawing_mtx(1:2,2),'-k','LineWidth',jammer_linewidth/2,'EraseMode','normal');
    h_jammer_crosshair_2=plot(jammer_crosshair_drawing_mtx(3:4,1),jammer_crosshair_drawing_mtx(3:4,2),'-k','LineWidth',jammer_linewidth/2,'EraseMode','normal');
    h_jammer_plot=[h_jammer_plot h_jammer_crosshair_1 h_jammer_crosshair_2];
else
    set(h_jammer_plot(1),'XData',jammer_drawing_mtx(:,1),'YData',jammer_drawing_mtx(:,2));
    set(h_jammer_plot(2),'XData',jammer_crosshair_drawing_mtx(1:2,1),'YData',jammer_crosshair_drawing_mtx(1:2,2));
    set(h_jammer_plot(3),'XData',jammer_crosshair_drawing_mtx(3:4,1),'YData',jammer_crosshair_drawing_mtx(3:4,2));
end


end     %   End function



function [ h_uav_jammer_vec_plot  ] = plot_uav_jammer_vec(x_t_vec,x_vec,h_uav_jammer_vec_plot)
%   This function plots the vector between the UAV and the jammer for
%   visualisation. The function should be called with graphics window open

global plot_scaling

uav_jammer_vec_p(1,:)=x_t_vec;
uav_jammer_vec_p(2,:)=x_vec;

hold on
if isempty(h_uav_jammer_vec_plot)
    h_uav_jammer_vec_plot=plot(uav_jammer_vec_p(:,1)/plot_scaling,uav_jammer_vec_p(:,2)/plot_scaling,'-k','Linewidth',1,'EraseMode','normal');
else
    set(h_uav_jammer_vec_plot,'XData',uav_jammer_vec_p(:,1)/plot_scaling,'YData',uav_jammer_vec_p(:,2)/plot_scaling);   
end
end             %   End function


function [ h_r_est_plot  ] = plot_r_est(x_vec_all,r_est_l,r_est_h,h_r_est_plot)
%   This function plots the range estimation circles for
%   visualisation. The function should be called with graphics window open

global plot_scaling

theta_data=0:(pi/180):2*pi;

hold on
if isempty(h_r_est_plot)
    h_r_est_plot_l=plot((x_vec_all(:,1)+r_est_l(:,1).*cos(theta_data))/plot_scaling,(x_vec_all(:,2)+r_est_l(:,1).*sin(theta_data))/plot_scaling,'--m','Linewidth',2,'EraseMode','normal');
    h_r_est_plot_h=plot((x_vec_all(:,1)+r_est_h(:,1).*cos(theta_data))/plot_scaling,(x_vec_all(:,2)+r_est_h(:,1).*sin(theta_data))/plot_scaling,'--m','Linewidth',2,'EraseMode','normal');
    h_r_est_plot=[h_r_est_plot_l h_r_est_plot_h];
else
    set(h_r_est_plot(1),'XData',(x_vec_all(:,1)+r_est_l(:,1).*cos(theta_data))/plot_scaling,'YData',(x_vec_all(:,2)+r_est_l(:,1).*sin(theta_data))/plot_scaling);   
    set(h_r_est_plot(2),'XData',(x_vec_all(:,1)+r_est_h(:,1).*cos(theta_data))/plot_scaling,'YData',(x_vec_all(:,2)+r_est_h(:,1).*sin(theta_data))/plot_scaling);
end

end             %   End function



function [ h_geo_est_plot  ] = plot_geo_est(C_loop,centre_geo_circle,radius_geo_circle,h_geo_est_plot)
%   This function plots the geolocation circles from iso-(range ratios) for
%   visualisation. The function should be called with graphics window open

global plot_scaling N_loops_fb

theta_data=0:(pi/180):2*pi;

hold on
if isempty(h_geo_est_plot)
    h_geo_est_plot=plot((centre_geo_circle(:,1)+radius_geo_circle(:,1).*cos(theta_data))/plot_scaling,(centre_geo_circle(:,2)+radius_geo_circle(:,1).*sin(theta_data))/plot_scaling,'-m','Linewidth',2,'EraseMode','normal');
    
else
    set(h_geo_est_plot,'XData',(centre_geo_circle(:,1)+radius_geo_circle(:,1).*cos(theta_data))/plot_scaling,'YData',(centre_geo_circle(:,2)+radius_geo_circle(:,1).*sin(theta_data))/plot_scaling);   
    if (rem(C_loop,floor(N_loops_fb/0.01))==0)
        plot((centre_geo_circle(:,1)+radius_geo_circle(:,1).*cos(theta_data))/plot_scaling,(centre_geo_circle(:,2)+radius_geo_circle(:,1).*sin(theta_data))/plot_scaling,'-m','Linewidth',2);
    end
end
end             %   End function



function [ h_state_pos_est  ] = plot_state_pos_est(x_state,h_state_pos_est)
%   This function plots the location of the position estimation provided by
%   

global plot_scaling x_bnd y_bnd

hold on
if isempty(h_state_pos_est)
    h_state_pos_est_1=plot([x_state(1,:)-(x_bnd*3/100) x_state(1,:)+(x_bnd*3/100)]/plot_scaling,[x_state(2,:) x_state(2,:)]/plot_scaling,'-.c','Linewidth',2,'EraseMode','normal');
    h_state_pos_est_2=plot([x_state(1,:) x_state(1,:)]/plot_scaling,[x_state(2,:)-(y_bnd*3/100) x_state(2,:)+(y_bnd*3/100)]/plot_scaling,'-.c','Linewidth',2,'EraseMode','normal');
    h_state_pos_est=[h_state_pos_est_1 h_state_pos_est_2];
else
    set(h_state_pos_est(1),'XData',[x_state(1,:)-(x_bnd*3/100) x_state(1,:)+(x_bnd*3/100)]/plot_scaling,'YData',[x_state(2,:) x_state(2,:)]/plot_scaling);   
    set(h_state_pos_est(2),'XData',[x_state(1,:) x_state(1,:)]/plot_scaling,'YData',[x_state(2,:)-(y_bnd*3/100) x_state(2,:)+(y_bnd*3/100)]/plot_scaling);
end

end             %   End function



function [  ] = plot_state_trace(x_state,C_loop,k_obs,N_loops_fb)
%   This function plots the trace of the position estimation provided by
%   the filter (x_state)
%   This plots the line between previous state (k-1) and current state (k)

global plot_scaling

hold on

k_plot=(k_obs+floor((1/100)*N_loops_fb));

if (C_loop>k_plot)
    plot(x_state(1,C_loop-1:C_loop)/plot_scaling,x_state(2,C_loop-1:C_loop)/plot_scaling,'-.c','Linewidth',4,'EraseMode','normal'); 
end

end             %   End function




function [h_P_cov_ell]=plot_P_cov_ell(x_state,P_cov,p_e,h_P_cov_ell)
%   This function aims at plotting the covariance ellipse
global plot_scaling 


sca=-2*log(1-p_e);
P_cov_eigs=eig(P_cov);
lambda_1=max(P_cov_eigs);
lambda_2=min(P_cov_eigs);

a=sqrt(sca*lambda_1);
b=sqrt(sca*lambda_2);

phi=(1/2)*atan((2*P_cov(1,2))/(P_cov(1,1)-P_cov(2,2)));

theta=linspace(0,2*pi,360);
X=x_state(1,:)+a*cos(theta)*cos(phi)-b*sin(theta)*sin(phi);
Y=x_state(2,:)+a*cos(theta)*sin(phi)+b*sin(theta)*cos(phi);

hold on

if isempty(h_P_cov_ell)
    h_P_cov_ell=plot(X/plot_scaling,Y/plot_scaling,'b-','Linewidth',2);
else
    set(h_P_cov_ell,'XData',X/plot_scaling,'YData',Y/plot_scaling);
end

end



function [ h_loiter_plot  ] = plot_loiter(x_state,r_d,h_loiter_plot)
%   This function plots the geolocation circles from iso-(range ratios) for
%   visualisation. The function should be called with graphics window open

global plot_scaling

theta_data=0:(pi/180):2*pi;

hold on
if isempty(h_loiter_plot)
    h_loiter_plot=plot((x_state(1,:)+r_d*cos(theta_data))/plot_scaling,(x_state(2,:)+r_d*sin(theta_data))/plot_scaling,'--g','Linewidth',2,'EraseMode','normal');
    
else
    set(h_loiter_plot,'XData',(x_state(1,:)+r_d*cos(theta_data))/plot_scaling,'YData',(x_state(2,:)+r_d*sin(theta_data))/plot_scaling);   
end
end             %   End function


























