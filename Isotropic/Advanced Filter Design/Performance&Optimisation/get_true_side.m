function [ true_side ] = get_true_side(x_t_vec,x_vec,psi_0)
%   This function aims at finding whether the jammer is port or starboard
%   of the UAV. It uses true jammer and jammer position. The side is found
%   using a simple 2D crossproduct

%   This function is most useful for straight flyby

    uav_jammer_vec_ini=(((x_t_vec-x_vec)')/(sqrt((x_t_vec-x_vec)*(x_t_vec-x_vec)')));                           %   uav-->jammer vector normalised
    uav_directive_vec=[cos(psi_0) sin(psi_0)]';                                                                 %   similar to x_vec_dot but normalised

    cross_prod=uav_directive_vec(1,1)*uav_jammer_vec_ini(2,1)-uav_directive_vec(2,1)*uav_jammer_vec_ini(1,1);   %   crossproduct of uav_directive_vec by uav_jammer_vec_ini  
    
    tolerance=1e-3;                                                                                             %   Colinearity tolerance
    
    if (cross_prod>tolerance)
        true_side=1;                                                                                            %   Port
    elseif (cross_prod<tolerance)
        true_side=0;                                                                                            %   Starboard
    else
        disp('UAV is passing through jammer within tolerance');
        disp('Simulation stops: check guidance or re-run simulation');
        return
    end
end                                                                                                             %   End function

