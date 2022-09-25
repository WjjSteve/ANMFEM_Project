% This is the code that implement the GFEM
close all;
clear;

g = @circleg;
hmax = 1/16;
[p,e,t] = initmesh(g,'hmax',hmax);
% get the trangular elements of the domain

n = size(p,2);
u = zeros(n,1);
for i = 1:n
    u(i,1) = u_initial(p(1,i),p(2,i));
end
u0 = u;

CFL = 0.5;
T = 1;
nor = 2*pi*max(sqrt((p(1,:).*p(1,:)+p(2,:).*p(2,:))));
dt = CFL*hmax/nor;
% Time steps

bx1 = -2*pi*p(2,:);
bx2 = 2*pi*p(1,:);
% vector field f'(u)

R = ReactMat2D(p,t);
C = ConvectMat2D(p,t,bx1,bx2);

for time = 0:dt:T
    LH = R/dt+C/2;
    RH = (R/dt-C/2)*u;
    I = eye(length(p));
    LH(e(1,:),:) = I(e(1 ,:) ,:); 
    RH(e(1 ,:)) = 0;
    u = LH\RH;
end

e = u0-u;
L2E = sqrt(e'*R*e);

%figure(3)
%plot3(p(1,:),p(2,:),u,'o');
%pdesurf(p,t,u);
%pdeplot(p,e,t,"XYData",u0-u); 
%pbaspect([1 1 1]);
%caxis([-0.06 0.06]);
%colormap turbo;
%title({'e = u-u_h', ['T = ', num2str(T), ',h_{max}=',num2str(rats(hmax))]});

