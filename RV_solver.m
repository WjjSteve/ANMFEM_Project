% This is the code that implement the RV
close all;
clear;

tic
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
% initial condition

CFL = 0.5;
T = 1;
nor = 2*pi*max(sqrt((p(1,:).*p(1,:)+p(2,:).*p(2,:))));
dt = CFL*hmax/nor;
% Time steps

bx1 = -2*pi*p(2,:);
bx2 = 2*pi*p(1,:);
% vector field f'(u)

nt = size(t,2);
eps = zeros(nt,1);
Cvel = 0.25;
Crv = 1.0;

R = ReactMat2D(p,t);
C = ConvectMat2D(p,t,bx1,bx2);

LH = R/dt+C/2;
RH = (R/dt-C/2)*u;
I = eye(length(p));
LH(e(1,:),:) = I(e(1 ,:) ,:); 
RH(e(1 ,:)) = 0;
up = u;
u = LH\RH;

Res = zeros(n,1);

for time = 0:1:10
    LH = R;
    RH = R*(u-up)/dt+C*u;
    Res = LH\RH; 
    Res = Res/max(abs(u-mean(u)));
    for K = 1:nt 
        beta_K = 2*pi*max(sqrt((p(1,t(1:3,K)).*p(1,t(1:3,K))+p(2,t(1:3,K)).*p(2,t(1:3,K)))));
        Res_K = max(abs(Res(t(1:3,K))));
        h = ComputeDiameter(p(1,t(1:3,K)),p(2,t(1:3,K)));
        eps(K) = min(Cvel*h*beta_K,Crv*h^2*Res_K);
    end
    A = StiffMat2D(p,t,eps);
    LH = R/dt+C/2+A/2;
    RH = (R/dt-C/2-A/2)*u;
    I = eye(length(p));
    LH(e(1,:),:) = I(e(1 ,:) ,:); 
    RH(e(1 ,:)) = 0;
    up = u;
    u = LH\RH;
end

toc
%L2E = (u-u0)'*R*(u-u0);
%t =[0.590081 2.804564 23.235106 404.970036];