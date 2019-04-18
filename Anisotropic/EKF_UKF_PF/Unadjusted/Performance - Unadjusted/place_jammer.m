function [ x_t_vec ] = place_jammer()
%   This function aims at generating the jammer position. This position
%   should lie within a square centred in the search area. The position is
%   generated randomly. It is not known by the UAV
global x_bnd y_bnd

        %   Define jammer square field: change if needed
        jam_x_bnd_l=(x_bnd/3);
        jam_x_bnd_h=(2*x_bnd/3);
        jam_y_bnd_l=(y_bnd/3);
        jam_y_bnd_h=(2*y_bnd/3);
        %   Generate random position between bounds
        x_t=jam_x_bnd_l+(jam_x_bnd_h-jam_x_bnd_l)*rand(1,1);        %   [m]
        y_t=jam_y_bnd_l+(jam_y_bnd_h-jam_y_bnd_l)*rand(1,1);        %   [m]
        %   Vectorise result and return
        x_t_vec=[x_t y_t];

end

