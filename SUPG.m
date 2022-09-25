% This is the code that implement the SUPG
close all;
clear;

tic
g = @circleg;
hmax = 1/32;
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

nt = size(t,2);
tau = ((2/dt)^2+(2*nor/hmax)^2)^(-1/2);
beta = zeros(nt,1);

R = ReactMat2D(p,t);
C = ConvectMat2D(p,t,bx1,bx2);
S = SDAssembler2D(p,t,bx1,bx2); 

for time = 0:dt:T
    LH = R/dt+C'*tau/dt+1/2*C+S*tau/2;
    RH = R/dt+C'*tau/dt-1/2*C-S*tau/2;
    RH = RH*u;
    I = eye(length(p));
    LH(e(1,:),:) = I(e(1 ,:) ,:); 
    RH(e(1 ,:)) = 0;
    u = LH\RH;
end
toc

%L2E = (u-u0)'*R*(u-u0);
%t = [0.046007 0.193307 1.219309 16.465116]
