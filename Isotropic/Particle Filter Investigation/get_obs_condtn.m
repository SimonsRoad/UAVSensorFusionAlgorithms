function [ obs_condtn ] = get_obs_condtn(x_Ca,y_Ca,x_Cb,y_Cb,r_Ca,r_Cb)
%   This function aims at finding whether  2 circles have an
%   intersection at two points. The return of this function is a
%   discriminant that is checked for positivity.
%   
%   This discriminant is obtained after developing algebraically the two
%   circles equation, finding the line passing through intersection
%   (supposing it exists) and check whether this line has indeed
%   intersections with circles

%   This function is used so as to get when subsequent geolocation circles
%   intersect 

    a=2*(x_Ca-x_Cb);
    b=2*(y_Ca-y_Cb);
    c=((r_Cb^2-r_Ca^2)+(x_Ca^2+y_Ca^2)-(x_Cb^2+y_Cb^2));

    if (b~=0)
        A=(1+(a/b)^2);
        B=(((2*a*y_Ca)/b)-(2*x_Ca)-((2*a*c)/((b^2))));
        C=x_Ca^2+y_Ca^2-r_Ca^2+(c/b)^2-2*((y_Ca*c)/(b));
        obs_condtn=B*B-4*A*C;

    elseif (b==0)
        A=1;
        B=-2*y_Ca;
        C=x_Ca^2+y_Ca^2-r_Ca^2+(c/a)^2-2*((x_Ca*c)/(a));
        obs_condtn=B*B-4*A*C; 
    end

end                                                                                                             %   End function

