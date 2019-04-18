function [ OC,r ] = get_geo_data(OM_1,OM_2,alpha)
%   This function aims at calculating the radius of the geolocating circle
%   at time k. The expression is complex hence the need for a subroutine

%   alpha is the ratio of powers
%   OM_1 is the absolute vector to point M1
%   OM_2 is the absolute vector to point M1
%   O is origin (0,0)
%   C is centre position

%   See report for equations

x_1=OM_1(:,1);
y_1=OM_1(:,2);

x_2=OM_2(:,1);
y_2=OM_2(:,2);

OC=zeros(1,2);
OC(1,1)=(x_1-alpha*x_2)/(1-alpha);     %   x of circle
OC(1,2)=(y_1-alpha*y_2)/(1-alpha);     %   y of circle


x_c=OC(:,1);
y_c=OC(:,2);

r_square=x_c*x_c+y_c*y_c+((alpha*(x_2*x_2+y_2*y_2)-(x_1*x_1+y_1*y_1))/(1-alpha));   %   Circle radius

if (r_square<0)
    r=0;
end

r=sqrt(r_square);

end                                                                                                             %   End function

