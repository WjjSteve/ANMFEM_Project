function u0 = u_initial(x1,x2)

x1_ini = 0.3;
x2_ini = 0.0;
r0 = 0.25;

%u0 = 1/2*(1-tanh(((x1-x1_ini)^2+(x2-x2_ini)^2)/r0/r0-1));
if (((x1-x1_ini)^2+(x2-x2_ini)^2)<=r0^2)
    u0 = 1;
else
    u0 = 0;
end

end