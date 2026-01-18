%Transition after unexpected decrease in aggregate productivity ("MIT shock")

clear all; clc; close all;

tic;

ga = 2;
rho = 0.05;
d = 0.05;
al = 1/3;
Aprod = 0.1;
z1 = 1;
z2 = 2*z1;
z = [z1,z2];
la1 = 1/3;
la2 = 1/3;
la = [la1,la2];
z_ave = (z1*la2 + z2*la1)/(la1 + la2);

T = 200;
N = 400;
dt = T/N;
time = (0:N-1)*dt;
max_price_it = 300;
convergence_criterion = 10^(-5);
relax = 0.1;

%construct TFP sequence
corr = 0.8;
nu = 1-corr;
Aprod_t = zeros(N,1);
Aprod_t(1)=.97*Aprod;
for n=1:N-1
Aprod_t(n+1) = dt*nu*(Aprod-Aprod_t(n)) + Aprod_t(n);
end
plot(time,Aprod_t)
xlim([0 40])


I= 1000;
amin = -0.8;
amax = 20;
a = linspace(amin,amax,I)';
da = (amax-amin)/(I-1);

aa = [a,a];
zz = ones(I,1)*z;


maxit= 100;
crit = 10^(-6);
Delta = 1000;

dVf = zeros(I,2);
dVb = zeros(I,2);
c = zeros(I,2);

Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

Ir = 40;
crit_S = 10^(-5);


rmax = 0.049;
r = 0.04;
w = 0.05;

v0(:,1) = (w*z(1) + r.*a).^(1-ga)/(1-ga)/rho;
v0(:,2) = (w*z(2) + r.*a).^(1-ga)/(1-ga)/rho;

r0 = 0.03;
rmin = 0.01;
rmax = 0.99*rho;


%%%%%%%%%%%%%%%%
% STEADY STATE %
%%%%%%%%%%%%%%%%


for ir=1:Ir;

r_r(ir)=r;
rmin_r(ir)=rmin;
rmax_r(ir)=rmax;

KD(ir) = (al*Aprod/(r + d))^(1/(1-al))*z_ave;
w = (1-al)*Aprod*KD(ir).^al*z_ave^(-al);

if ir>1
v0 = V_r(:,:,ir-1);
end

v = v0;

for n=1:maxit
    V = v;
    V_n(:,:,n)=V;
    % forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = (w*z + r.*amax).^(-ga); %will never be used, but impose state constraint a<=amax just in case
    % backward difference
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:) = (w*z + r.*amin).^(-ga); %state constraint boundary condition
    
    %consumption and savings with forward difference
    cf = dVf.^(-1/ga);
    ssf = w*zz + r.*aa - cf;
    %consumption and savings with backward difference
    cb = dVb.^(-1/ga);
    ssb = w*zz + r.*aa - cb;
    %consumption and derivative of value function at steady state
    c0 = w*zz + r.*aa;
    dV0 = c0.^(-ga);
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If = ssf > 0; %positive drift --> forward difference
    Ib = ssb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
    %make sure backward difference is used at amax
    %Ib(I,:) = 1; If(I,:) = 0;
    %STATE CONSTRAINT at amin: USE BOUNDARY CONDITION UNLESS muf > 0:
    %already taken care of automatically
    
    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %important to include third term
    c = dV_Upwind.^(-1/ga);
    u = c.^(1-ga)/(1-ga);
    
    %CONSTRUCT MATRIX
    X = -min(ssb,0)/da;
    Y = -max(ssf,0)/da + min(ssb,0)/da;
    Z = max(ssf,0)/da;
    
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;

    if max(abs(sum(A,2)))>10^(-9)
       disp('Improper Transition Matrix')
       break
    end
    
    B = (1/Delta + rho)*speye(2*I) - A;

    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)];
    
    b = u_stacked + V_stacked/Delta;
    V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
    
    V = [V_stacked(1:I),V_stacked(I+1:2*I)];
    
    Vchange = V - v;
    v = V;

    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOKKER-PLANCK EQUATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%
AT = A';
b = zeros(2*I,1);

%need to fix one value, otherwise matrix is singular
i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
AT(i_fix,:) = row;

%Solve linear system
gg = AT\b;
g_sum = gg'*ones(2*I,1)*da;
gg = gg./g_sum;

g = [gg(1:I),gg(I+1:2*I)];

check1 = g(:,1)'*ones(I,1)*da;
check2 = g(:,2)'*ones(I,1)*da;

g_r(:,:,ir) = g;
adot(:,:,ir) = w*zz + r.*aa - c;
V_r(:,:,ir) = V;

KS(ir) = g(:,1)'*a*da + g(:,2)'*a*da;
S(ir) = KS(ir) - KD(ir);

%UPDATE INTEREST RATE
if S(ir)>crit_S
    disp('Excess Supply')
    rmax = r;
    r = 0.5*(r+rmin);
elseif S(ir)<-crit_S;
    disp('Excess Demand')
    rmin = r;
    r = 0.5*(r+rmax);
elseif abs(S(ir))<crit_S;
    display('Equilibrium Found, Interest rate =')
    disp(r)
    break
end

end


%save some objects
v_st = v;
gg_st = gg;
K_st = KS(ir);
w_st = w;
r_st = r;
C_st = gg'*reshape(c,I*2,1)*da;
Aprod_st = Aprod;

%%%%%%%%%%%%%%%%%%%%%%%
% TRANSITION DYNAMICS %
%%%%%%%%%%%%%%%%%%%%%%%
gg0 = gg_st;

clear Delta r w Aprod gg

%preallocation
gg = cell(N+1,1);
K_t = zeros(N,1);
K_out = zeros(N,1);
r_t = zeros(N,1);
w_t = zeros(N,1);

%initial guess
K_t = K_st*ones(N,1);

for it=1:max_price_it
disp('PRICE ITERATION = ')
disp(it)

w_t = (1-al)*Aprod_t.*K_t.^al*z_ave^(-al);
r_t = al*Aprod_t.*K_t.^(al-1)*z_ave^(1-al) - d;

V = v_st;


for n=N:-1:1
    % forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = (w_t(n)*z + r_t(n).*amax).^(-ga); %will never be used, but impose state constraint a<=amax just in case
    % backward difference
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:) = (w_t(n)*z + r_t(n).*amin).^(-ga); %state constraint boundary condition
    
    %consumption and savings with forward difference
    cf = dVf.^(-1/ga);
    ssf = w_t(n)*zz + r_t(n).*aa - cf;
    %consumption and savings with backward difference
    cb = dVb.^(-1/ga);
    ssb = w_t(n)*zz + r_t(n).*aa - cb;
    %consumption and derivative of value function at steady state
    c0 = w_t(n)*zz + r_t(n).*aa;
    dV0 = c0.^(-ga);
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If = ssf > 0; %positive drift --> forward difference
    Ib = ssb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
    %make sure backward difference is used at amax
    %Ib(I,:) = 1; If(I,:) = 0;
    %STATE CONSTRAINT at amin: USE BOUNDARY CONDITION UNLESS muf > 0:
    %already taken care of automatically
    
    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %important to include third term
    c = dV_Upwind.^(-1/ga);
    u = c.^(1-ga)/(1-ga);
    c_t(:,:,n)=c;
    
    %CONSTRUCT MATRIX
    X = -min(ssb,0)/da;
    Y = -max(ssf,0)/da + min(ssb,0)/da;
    Z = max(ssf,0)/da;
    
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    
    if max(abs(sum(A,2)))>10^(-9)
       disp('Improper Transition Matrix')
       break
    end
    
    A_t{n}=A; %save for future reference
    
    B = (1/dt + rho)*speye(2*I) - A;

    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)];
    
    b = u_stacked + V_stacked/dt;
    V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
    
    V = [V_stacked(1:I),V_stacked(I+1:2*I)];
    

end
toc;

gg{1}=gg0;
for n=1:N
    AT=A_t{n}';
    %Implicit method in Updating Distribution.
    gg{n+1}= (speye(2*I) - AT*dt)\gg{n};
    %gg{n+1}=gg{n}+AT*gg{n}*dt; %This is the explicit method.
    %check(n) = gg(:,n)'*ones(2*I,1)*da;
    K_out(n) = gg{n}(1:I)'*a*da + gg{n}(I+1:2*I)'*a*da;
    C_t(n) = gg{n}'*reshape(c_t(:,:,n),I*2,1)*da;
end

dist_it(it) = max(abs(K_out - K_t));
figure(1)
plot(dist_it)
disp(dist_it(it))
xlabel('Iteration')
title('Convergence Criterion')

K_t = relax.*K_out +(1-relax).*K_t;

if dist_it(it)<convergence_criterion
    disp('Equilibrium Found')
    break
end

end

%%%%%%%%%%%%%%%%%%%%%%%
% INEQUALITY MEASURES %
%%%%%%%%%%%%%%%%%%%%%%%

for n=1:N
    Wealth_neg_t(n) = gg{n}(1:I)'*min(a,0)*da + gg{n}(I+1:2*I)'*min(a,0)*da;
    g_a_cont = gg{n}(1:I)+gg{n}(I+1:2*I);
    %Discrete Gini to check
    g_a = g_a_cont*da;
    S_a = cumsum(g_a.*a)/sum(g_a.*a);
    trapez_a = (1/2)*(S_a(1)*g_a(1) + sum((S_a(2:I) + S_a(1:I-1)).*g_a(2:I)));
    Wealth_Gini_t(n) = 1 - 2*trapez_a;

    %Gini of capital income
    yk = r_t(n).*a;
    g_yk = g_a;

    S_yk = cumsum(g_yk.*yk)/sum(g_yk.*yk);
    trapez_yk = (1/2)*(S_yk(1)*g_yk(1) + sum((S_yk(2:I) + S_yk(1:I-1)).*g_yk(2:I)));
    CapInc_Gini_t(n) = 1 - 2*trapez_yk;

    %Gini of total income
    y = w_t(n)*zz + r_t(n).*aa;
    Ny = 2*I;
    yy =reshape(y,Ny,1);
    [yy,index] = sort(yy);
    g_y = gg{n}(index)*da;

    S_y = cumsum(g_y.*yy)/sum(g_y.*yy);
    trapez_y = (1/2)*(S_y(1)*g_y(1) + sum((S_y(2:Ny) + S_y(1:Ny-1)).*g_y(2:Ny)));
    Income_Gini_t(n) = 1 - 2*trapez_y;

    G_y = cumsum(g_y);
    G_a = cumsum(g_a);

    %Top 10% Income Share
    p1 = 0.1;
    [obj index] = min(abs((1-G_y) - p1));
    top_inc_t(n) = 1-S_y(index);

    %Top 10% Wealth Share
    p1 = 0.1;
    [obj index] = min(abs((1-G_a) - p1));
    top_wealth_t(n) = 1-S_a(index);
    
    neg_frac_t(n) = max(G_a(a<=0));

end


set(gcf,'PaperPosition',[0 0 15 10])
subplot(2,3,1)
set(gca,'FontSize',16)
plot(time,Aprod_t,time,Aprod_st*ones(N,1),'--','LineWidth',2)
xlim([0 100])
xlabel('Year')
title('Aggregate Productivity')

subplot(2,3,2)
set(gca,'FontSize',16)
plot(time,K_t,time,K_st*ones(N,1),'--','LineWidth',2)
xlim([0 100])
xlabel('Year')
title('Aggregate Capital Stock')

subplot(2,3,3)
set(gca,'FontSize',16)
plot(time,w_t,time,w_st*ones(N,1),'--','LineWidth',2)
xlim([0 100])
xlabel('Year')
title('Wage')

subplot(2,3,4)
set(gca,'FontSize',16)
plot(time,r_t,time,r_st*ones(N,1),'--','LineWidth',2)
xlim([0 100])
xlabel('Year')
title('Interest Rate')

subplot(2,3,5)
set(gca,'FontSize',16)
plot(time,Wealth_Gini_t,time,Wealth_Gini_t(N)*ones(N,1),'--','LineWidth',2)
xlim([0 100])
xlabel('Year')
title('Wealth Gini')

subplot(2,3,6)
set(gca,'FontSize',16)
plot(time,Income_Gini_t,time,Income_Gini_t(N)*ones(N,1),'--','LineWidth',2)
xlim([0 100])
xlabel('Year')
title('Income Gini')

print -depsc aiyagari_poisson_MITshock.eps