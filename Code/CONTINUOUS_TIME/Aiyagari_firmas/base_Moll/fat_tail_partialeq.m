clear all; clc; close all;

tic;

ga = 2;
rho = 0.05;
z1 = .01;
z2 = 3*z1;
z = [z1,z2];
la1 = 0.5;
la2 = 0.5;
la = [la1,la2];

w = 3;
R = 0.051;
r = 0.041;
phi = 0.3;

zeta = 1.5;
sig2 = (zeta/ga + 1)/2*(R-r)^2/(rho - r); %pick sig2 so that zeta=1.5
sig = sig2^(1/2);

I= 5000;
amin = -phi;
amax = 1000;

%non-uniform grid
x = linspace(0,1,I)';
coeff = 5; power = 10;
xx  = x + coeff*x.^power;
xmax = max(xx); xmin = min(xx);
a = (amax-amin)/(xmax - xmin)*xx + amin;
daf = ones(I,1);
dab = ones(I,1);
daf(1:I-1) = a(2:I)-a(1:I-1);
dab(2:I) = a(2:I)-a(1:I-1);
daf(I)=daf(I-1); dab(1)=dab(2);

aa = [a,a];
daaf = daf*ones(1,2);
daab = dab*ones(1,2);

%objects for approximation of second derivatives
denom = 0.5*(daaf + daab).*daab.*daaf;
weightf = daab./denom;
weight0 = -(daab + daaf)./denom;
weightb = daaf./denom;


zz = ones(I,1)*z;
kk = ones(I,2);

maxit= 30;
crit = 10^(-6);
Delta = 1000;

dVf = zeros(I,2);
dVb = zeros(I,2);
dV0 = zeros(I,2);
dV2 = zeros(I,2);
c = zeros(I,2);
c0 = zeros(I,2);

Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

%INITIAL GUESS
v0(:,1) = (w*z(1) + r*a).^(1-ga)/(1-ga)/rho;
v0(:,2) = (w*z(2) + r*a).^(1-ga)/(1-ga)/rho;

% v0(:,1) = (w*z(1) + r*a + (R-r)^2/(s*sig2)*a).^(1-s)/(1-s)/rho;
% v0(:,2) = (w*z(2) + r*a + (R-r)^2/(s*sig2)*a).^(1-s)/(1-s)/rho;

v = v0;

for n=1:maxit
    disp(n)
    V = v; 
    V_n(:,:,n)=V;
    % forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))./(aa(2:I,:) - aa(1:I-1,:));
    dVf(I,:) = (w*z + r.*amax + (R-r)^2/(ga*sig2)*amax).^(-ga); %will never be used
    % backward difference
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))./(aa(2:I,:) - aa(1:I-1,:));
    dVb(1,:) = (w*z + r.*amin).^(-ga); %state constraint boundary condition

    %second derivative: approximation only differs at amax
    dV2b(2:I-1,:) = (daab(2:I-1,:).*V(3:I,:) - (daab(2:I-1,:)+daaf(2:I-1,:)).*V(2:I-1,:) + daaf(2:I-1,:).*V(1:I-2,:))./denom(2:I-1,:);
    dV2f(2:I-1,:) = (daab(2:I-1,:).*V(3:I,:) - (daab(2:I-1,:)+daaf(2:I-1,:)).*V(2:I-1,:) + daaf(2:I-1,:).*V(1:I-2,:))./denom(2:I-1,:);
    dV2b(I,:) = -ga*dVb(I,:)/amax;
    dV2f(I,:) = -ga*dVf(I,:)/amax;
     
    %consumption and savings with forward difference
    cf = max(dVf,10^(-10)).^(-1/ga);
    kf = max(- dVf./dV2f.*(R-r)/sig2,0);
    kf = min(kf,aa+phi);
    ssf = w*zz + (R-r).*kf + r.*aa - cf;
    %consumption and savings with backward difference
    cb = max(dVb,10^(-10)).^(-1/ga);
    kb = max(- dVb./dV2b.*(R-r)/sig2,0);
    kb = min(kb,aa+phi);
    ssb = w*zz + (R-r).*kb + r.*aa - cb;
    %consumption and derivative of value function at steady state
    k0 = (kb + kf)/2; %could do something fancier here but the following seems to work well. And more importantly, it's fast
    c0 = w*zz + (R-r).*k0 + r.*aa;
    dV0 = max(c0,10^(-10)).^(-ga);
   
    
    % dV_upwind makes a choice of forward or backward differences based on the sign of the drift    
    If = real(ssf) > 10^(-12); %positive drift --> forward difference
    Ib = (real(ssb) < -10^(-12)).*(1-If); %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
   
    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %important to include third term
    c = max(dV_Upwind,10^(-10)).^(-1/ga);
    u = c.^(1-ga)/(1-ga);
    k = max(-dV_Upwind./dV2b.*(R-r)/sig2,0);
    k = min(k,aa+phi);
    
   
    %CONSTRUCT MATRIX
    X = -Ib.*ssb./daab + sig2/2.*k.^2.*weightb;
    Y = - If.*ssf./daaf + Ib.*ssb./daab + sig2/2.*k.^2.*weight0;
    Z = If.*ssf./daaf + sig2/2.*k.^2.*weightf; 
    
    xi = -amax*(R-r)^2/(2*ga*sig2);
    X(I,:) = -min(ssb(I,:),0)./daab(I,:) - xi./daab(I,:);
    Y(I,:) = -max(ssf(I,:),0)./daaf(I,:) + min(ssb(I,:),0)./daab(I,:) + xi./daab(I,:);
    Z(I,:) = max(ssf(I,:),0)./daaf(I,:);
    
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A1(I,I) = Y(I,1) + Z(I,1);
    A2(I,I) = Y(I,2) + Z(I,2);
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

%RECOMPUTE TRANSITION MATRIX WITH REFLECTING BARRIER AT amax
X = -min(ssb,0)./daab + sig2/2.*k.^2.*weightb;
Y = -max(ssf,0)./daaf + min(ssb,0)./daab + sig2/2.*k.^2.*weight0;
Z = max(ssf,0)./daaf + sig2/2.*k.^2.*weightf; 

A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
A1(I,I) = Y(I,1) + Z(I,1);
A2(I,I) = Y(I,2) + Z(I,2);
A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;

%WORK WITH RESCALED DENSITY \tilde{g} BELOW
da_tilde = 0.5*(dab + daf);
da_tilde(1) = 0.5*daf(1); da_tilde(I) = 0.5*dab(I);
da_stacked = [da_tilde;da_tilde];
grid_diag = spdiags(da_stacked,0,2*I,2*I);

AT = A';
b = zeros(2*I,1);

%need to fix one value, otherwise matrix is singular
i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
AT(i_fix,:) = row;

%Solve linear system
g_tilde = AT\b;

%rescale \tilde{g} so that it sums to 1
g_sum = g_tilde'*ones(2*I,1);
g_tilde = g_tilde./g_sum;

gg = grid_diag\g_tilde; %convert from \tilde{g} to g

g = [gg(1:I),gg(I+1:2*I)];

check1 = g(:,1)'*da_tilde;
check2 = g(:,2)'*da_tilde;

%CALCULATE THEORETICAL POWER LAW EXPONENT
zeta = ga*(2*sig2*(rho - r)/(R-r)^2 -1);

adot = w*zz + (R-r)*k + r.*aa - c;

risky_share = (R-r)/(ga*sig2);
cslope = (rho - (1-ga)*r)/ga - (1-ga)/(2*ga)*(R-r)^2/(ga*sig2);
sbar = (r-rho)/ga + (1+ga)/(2*ga)*(R-r)^2/(ga*sig2);

plot(a,k,a,risky_share.*a)
plot(a,c,a,cslope.*a)

%COMPUTE DISTRIBUTION OF x=log(a)
for i=1:I
    G(i,1) = sum(g(1:i,1).*da_tilde(1:i));
    G(i,2) = sum(g(1:i,2).*da_tilde(1:i));
end

f = zeros(I,2);
x = log(max(a,0));
dx = zeros(I,1);
for i =2:I
dx(i) = x(i)-x(i-1);
f(i,1) = (G(i,1)-G(i-1,1))/dx(i);
f(i,2) = (G(i,2)-G(i-1,2))/dx(i);
end
f(1)=0;

xmin = log(1); xmax = log(amax);
ga = g(1:i,1) + g(1:i,2);
fx = f(:,1)+f(:,2);


amax1 = 40;
amin1 = amin-2;
figure(1)
h1 = plot(a,adot(:,1),'b',a,adot(:,2),'r',linspace(amin,amax,I),zeros(1,I),'k--','LineWidth',2)
set(gca,'FontSize',16)
legend(h1,'s_1(a)','s_2(a)','Location','NorthEast')
text(-0.555,-.085,'$\underline{a}$','FontSize',16,'interpreter','latex')
line([amin amin], [-.08 .06],'Color','Black','LineStyle','--')
xlabel('Wealth, $a$','interpreter','latex')
ylabel('Savings, $s_i(a)$','interpreter','latex')
xlim([amin1 amax1])
set(gca, 'XTick', []);
set(gca, 'YTick', []);
print -depsc sav_fat.eps

amax1 = 10;
amin1 = amin-0.3;
figure(2)
h1 = plot(a,k(:,1),'b',a,k(:,2),'r','LineWidth',2)
set(gca,'FontSize',16)
legend(h1,'k_1(a)','k_2(a)','Location','NorthWest')
text(-0.2,-.2,'$\underline{a}$','FontSize',16,'interpreter','latex')
line([amin amin], [0 7],'Color','Black','LineStyle','--')
xlabel('Wealth, $a$','interpreter','latex')
ylabel('Risky Assets, $k_i(a)$','interpreter','latex')
xlim([amin1 amax1])
set(gca, 'XTick', []);
set(gca, 'YTick', []);
print -depsc k_fat.eps

figure(3)
h1 = plot(a,g(:,1),'b',a,g(:,2),'r','LineWidth',2)
set(gca,'FontSize',16)
legend(h1,'g_1(a)','g_2(a)')
text(-0.555,-.015,'$\underline{a}$','FontSize',16,'interpreter','latex')
line([amin amin], [0 max(max(g))],'Color','Black','LineStyle','--')
xlabel('Wealth, $a$','interpreter','latex')
ylabel('Densities, $g_i(a)$','interpreter','latex')
xlim([amin1 amax1])
ylim([0 0.3])
set(gca, 'XTick', []);
set(gca, 'YTick', []);
print -depsc densities_fat.eps

figure(4)
h1 = plot(x,log(f(:,1)),'b',x,log(f(:,2)),'r',x,-zeta*x','k--','LineWidth',2)
set(gca,'FontSize',16)
legend(h1,'log(f_1(x))','log(f_2(x))','Slope = -\zeta')
xlabel('Log Wealth, $x=\log \ (a)$','interpreter','latex')
ylabel('$\log \ (f_i(x))$','interpreter','latex')
%xlim([1 log(13.5)])
ylim([-10 0])
set(gca, 'XTick', []);
set(gca, 'YTick', []);
print -depsc log_densities_fat.eps
