% codigo de pontus para Aiyagari base
clear;

% First solve the Aiyagari model

% Parameters

alpha = 1/3;
beta = 1.05^(-1/4);
delta = 0.025;
gamma = 2;
mu = 0.5;             
phi = 0;  
rho = 1-beta;
Gamma = inf;
opts.disp = 0;

Na = 2000;

tic

a_max = 500;        
a_min = phi;
% agrd = (exp(linspace(log(a_min+1-a_min),log(a_max+1-a_min),Na))-1+a_min)';
agrd = linspace(a_min,a_max,Na)';

daf = agrd(2:Na)-agrd(1:Na-1);
dab = -(agrd(1:Na-1)-agrd(2:Na));

% Transition matrix in discrete time: g_{t+1}=T*g_{t}

T = [0.95, 1-0.95; 0.8, 1-0.8];

% Convert to continuous time such that dg_t=T*g_t

T = [-T(1,2) T(2,1);T(1,2) -T(2,1)];

% Convert to a transition matrix used in HJB equation

B = [-speye(Na)*T(2,1),speye(Na)*T(2,1);speye(Na)*T(1,2),-speye(Na)*T(1,2)];

% The steady state values of g satisfy 0=T*g. That is g is an eigenvector
% associated with a zero eigenvalue normalized to sum to one. But since T
% has a zero eigenvalue it is also singular. Solution, add a unit
% eigenvalue

[g,l] = eigs(T+eye(2),1);

g=g/sum(g);

% Steady state employment rate

e = g(1);

% Tax rate to balance budget: (1-e)*mu*w=e*tau*w.

tau = (1-e)/e*mu;

% Utility and consumption functions

u = @(c) (c.^(1-gamma)-1)./(1-gamma);
up = @(c) c.^(-gamma);
c = @(dv) dv.^(-1/gamma);

% Matrices for derivatives

Df = -spdiags([ones(Na-1,1)./daf;1],0,Na,Na);
Db = spdiags([1;ones(Na-1,1)./dab],0,Na,Na);
Df(1:end-1,2:end)=Df(1:end-1,2:end)+spdiags(ones(Na-1,1)./daf,0,Na-1,Na-1);
Db(2:end,1:end-1)=Db(2:end,1:end-1)-spdiags(ones(Na-1,1)./dab,0,Na-1,Na-1);

% Ok, problem set up. Let's solve it

rh = 1/beta-1;
rl = 0;

% Guess for value function

r = (rh+rl)/2;

w = (1-alpha)*(((r+delta)/alpha )^(1/(alpha-1)))^(alpha);

ve = u((r*agrd+w*(1-tau)))./rho;
ve(end) = u(r*a_max+w*(1-tau))./rho;
vu = u((r*agrd+w*mu))./rho;
vu(1) = u(r*a_min+w*mu)./rho;

metric_r = 1;

while metric_r>1e-8

    r = (rh+rl)/2;

    w = (1-alpha)*(((r+delta)/alpha )^(1/(alpha-1)))^(alpha);

    metric = 1;

    while metric>1e-8

    % Calculate derivatives

    dvef = Df*ve;
    dveb = Db*ve;

    dvuf = Df*vu;
    dvub = Db*vu;

    dvef(end) = up(w*(1-tau)+r*a_max);
    dvuf(end) = up(w*mu+r*a_max);

    dveb(1) = up(w*(1-tau)+r*a_min);
    dvub(1) = up(w*mu+r*a_min);

    cef = c(dvef);
    cuf = c(dvuf);
    ceb = c(dveb);
    cub = c(dvub);

    sef = r*agrd+w*(1-tau)-cef;
    seb = r*agrd+w*(1-tau)-ceb;
    suf = r*agrd+w*mu-cuf;
    sub = r*agrd+w*mu-cub;

    Ief = sef>0;
    Ieb = seb<0;
    Iuf = suf>0;
    Iub = sub<0;

    I0e = (1-Ief-Ieb);
    I0u = (1-Iuf-Iub);

    ce = Ief.*cef+Ieb.*ceb+I0e.*(w*(1-tau)+r*agrd);
    cu = Iuf.*cuf+Iub.*cub+I0u.*(w*mu+r*agrd);
    
    S = [spdiags(Ief.*sef,0,Na,Na)+spdiags(Ieb.*seb,0,Na,Na),sparse(Na,Na);sparse(Na,Na),spdiags(Iuf.*suf,0,Na,Na)+spdiags(Iub.*sub,0,Na,Na)];
    D = [spdiags(Ief,0,Na,Na)*Df+spdiags(Ieb,0,Na,Na)*Db,sparse(Na,Na);sparse(Na,Na),spdiags(Iuf,0,Na,Na)*Df+spdiags(Iub,0,Na,Na)*Db];
    
    P = S*D+B;

    A = (1/Gamma+rho)*speye(Na*2)-P;
    b = [u(ce);u(cu)]+[ve;vu]./Gamma;

    v = A\b;

    ven = v(1:Na);
    vun = v(Na+1:end);

    metric = max(max(abs([ven-ve vun-vu])));

    ve = ven;
    vu = vun;

    end

    b = zeros(2*Na,1);

    % Solve P'G=0 for the distribution. Solve a linear system instead of
    % using eigs. This works since we can normalize one element of G.
    
    b(1)=.1;
    row = [1,zeros(1,2*Na-1)];
    P(:,1) = row;

    % Normalized, now we can solve the linear system
    
    G = (P')\b;
    G = G./(G'*ones(2*Na,1));
    G = [G(1:Na),G(Na+1:2*Na)];

    K = G(:,1)'*agrd+G(:,2)'*agrd;

    rp = alpha*(K/e)^(alpha-1)-delta;
    frs = rp-r;

    if frs>0
        rl = r;
    else
        rh = r;
    end
    
    metric_r = abs(frs);

end

toc

se = Ief.*sef+Ieb.*seb;
su = Iuf.*suf+Iub.*sub;

xSize = 19; 
ySize = 6.5;

alim = max(find(agrd<101));
alims = max(find(agrd<21));

figure;
subplot(1,2,1);
y1 = plot(agrd(1:alims),ce(1:alims),'linewidth',2);
hold on
y2 = plot(agrd(1:alims),cu(1:alims),'linewidth',2);
xlabel('Assets')
ylabel('Consumption')
legend1 = legend([y1 y2],'Employed','Unemployed');
set(legend1,'Location','southeast');
axis tight
set(gca,'FontSize',10)

subplot(1,2,2);
y1= plot(agrd(1:alims),se(1:alims),'linewidth',2);
hold on
y2=plot(agrd(1:alims),su(1:alims),'linewidth',2);
xlabel('Assets')
ylabel('Savings')
axis tight
set(gca,'FontSize',10)

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
     ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-2.7 ySize+0.5],'PaperPositionMode','auto')
 
print -dpdf hh.pdf

xSize = 10; 
ySize = 6.5;

figure;
y1= plot(agrd(1:alim),G(1:alim,1),'linewidth',2);
hold on
y2=plot(agrd(1:alim),G(1:alim,2),'linewidth',2);
xlabel('Assets')
ylabel('Density')
legend1 = legend([y1 y2],'Employed','Unemployed');
set(legend1,'Location','northeast');
axis tight
set(gca,'FontSize',10)

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
     ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-0 ySize+0.5],'PaperPositionMode','auto')
 
print -dpdf dist.pdf






