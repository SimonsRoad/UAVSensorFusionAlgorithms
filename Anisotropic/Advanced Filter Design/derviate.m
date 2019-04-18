
syms xk yk skx sky six siy h_0 G_t_1 G_t_2

jacobian(((G_t_2/G_t_1)*((sqrt((((xk-six)^2)+((yk-siy)^2)+(h_0^2)))^2)/(sqrt((((xk-skx)^2)+((yk-sky)^2)+(h_0^2))))^2)),[xk,yk])

% (sqrt((((xk-skx)^2)+((yk-sky)^2)))^2);
% [ ((2*skx - 2*xk)*((six - xk)^2 + (siy - yk)^2))/((skx - xk)^2 + (sky - yk)^2)^2 - (2*six - 2*xk)/((skx - xk)^2 + (sky - yk)^2),
%     ((2*sky - 2*yk)*((six - xk)^2 + (siy - yk)^2))/((skx - xk)^2 + (sky - yk)^2)^2 - (2*siy - 2*yk)/((skx - xk)^2 + (sky - yk)^2)]
% [ ((2*skx - 2*xk)*((six - xk)^2 + (siy - yk)^2))/((skx - xk)^2 + (sky - yk)^2)^2 - (2*six - 2*xk)/((skx - xk)^2 + (sky - yk)^2)
%     ((2*sky - 2*yk)*((six - xk)^2 + (siy - yk)^2))/((skx - xk)^2 + (sky - yk)^2)^2 - (2*siy - 2*yk)/((skx - xk)^2 + (sky - yk)^2)]
% [ ((2*sxk - 2*xk)*((sx1 - xk)^2 + (sy1 - yk)^2))/((sxk - xk)^2 + (syk - yk)^2)^2 - (2*sx1 - 2*xk)/((sxk - xk)^2 + (syk - yk)^2)
%     
%   ((2*syk - 2*yk)*((sx1 - xk)^2 + (sy1 - yk)^2))/((sxk - xk)^2 + (syk - yk)^2)^2 - (2*sy1 - 2*yk)/((sxk - xk)^2 + (syk - yk)^2)
% 
% % H=JH(xy1,xy2,X_s)
% % six=xy1(1);
% % siy=xy1(2);
% %  
% % skx=xy2(1);
% % sky=xy2(2);
% %  
% % xk=X_s(1);
% yk=X_s(2);
% Equation:
% ((((xk-six)^2)+((yk-siy)^2))/(((xk-skx)^2)+((yk-sky)^2)))
% [ ((2*skx - 2*xk)*((six - xk)^2 + (siy - yk)^2 + h_0^2))/((skx - xk)^2 + (sky - yk)^2 + h_0^2)^2 - (2*six - 2*xk)/((skx - xk)^2 + (sky - yk)^2 + h_0^2) 
%     ((2*sky - 2*yk)*((six - xk)^2 + (siy - yk)^2 + h_0^2))/((skx - xk)^2 + (sky - yk)^2 + h_0^2)^2 - (2*siy - 2*yk)/((skx - xk)^2 + (sky - yk)^2 + h_0^2)]
% 

[ (G_t_2*(2*skx - 2*xk)*((six - xk)^2 + (siy - yk)^2 + h_0^2))/(G_t_1*((skx - xk)^2 + (sky - yk)^2 + h_0^2)^2) - (G_t_2*(2*six - 2*xk))/(G_t_1*((skx - xk)^2 + (sky - yk)^2 + h_0^2))
    (G_t_2*(2*sky - 2*yk)*((six - xk)^2 + (siy - yk)^2 + h_0^2))/(G_t_1*((skx - xk)^2 + (sky - yk)^2 + h_0^2)^2) - (G_t_2*(2*siy - 2*yk))/(G_t_1*((skx - xk)^2 + (sky - yk)^2 + h_0^2))]

