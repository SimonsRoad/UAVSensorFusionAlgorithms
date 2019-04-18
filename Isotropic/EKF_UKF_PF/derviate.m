
%% Code for Differentiade the Measrument Function

syms xk yk skx sky six siy h_0

jacobian(((sqrt((((xk-six)^2)+((yk-siy)^2)+(h_0^2)))^2)/(sqrt((((xk-skx)^2)+((yk-sky)^2)+(h_0^2)))^2)),[xk,yk])

% Output:

% [ ((2*skx - 2*xk)*((six - xk)^2 + (siy - yk)^2 + h_0^2))/((skx - xk)^2 + (sky - yk)^2 + h_0^2)^2 - (2*six - 2*xk)/((skx - xk)^2 + (sky - yk)^2 + h_0^2) 
%     ((2*sky - 2*yk)*((six - xk)^2 + (siy - yk)^2 + h_0^2))/((skx - xk)^2 + (sky - yk)^2 + h_0^2)^2 - (2*siy - 2*yk)/((skx - xk)^2 + (sky - yk)^2 + h_0^2)]
