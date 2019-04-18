function [ x_vec, psi] = place_uav()
%   This function aims at generating the uav position and heading. This position
%   is set at the lower left-hand corner of the map. The position is
%   set to be fixed but the heading is randomly generated between [0 -
%   90°]. The UAV knows its position accurately (e.g GPS) and heading accurately (e.g
%   digital compass)
global x_bnd y_bnd d2r

        %   Define uav square field in the lower left-hand corner
        uav_x_bnd_l=(x_bnd/20);         %   Low square boundary x-axis
        uav_x_bnd_h=2*uav_x_bnd_l;      %   High square boundary x-axis
        uav_y_bnd_l=(y_bnd/20);
        uav_y_bnd_h=2*uav_y_bnd_l;
        %   Generate random position between position bounds
        x=uav_x_bnd_l+(uav_x_bnd_h-uav_x_bnd_l)*rand(1,1);        %   [m]
        y=uav_y_bnd_l+(uav_y_bnd_h-uav_y_bnd_l)*rand(1,1);        %   [m]
        %   Generate random heading between heading bounds
        heading_bnd_h=90*d2r;
        heading_bnd_l=0*d2r;
        psi=heading_bnd_l+(heading_bnd_h-heading_bnd_l)*rand(1,1); %    [rad]
        %   Vectorise result and return
        x_vec=[x y];
end
