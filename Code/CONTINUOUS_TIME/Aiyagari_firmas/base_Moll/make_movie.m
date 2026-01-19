load('aiyagari_transition');

%%%%% CLEAN FIGURE TO GET READY TO RECORD%%%%%
figure(1);
clf;
close 1;

%%%%% Set Parameters%%%%%
plotlength=40;  % End time of the video
icut = 50;      % Index of maximum wealth being plotted.

h=figure('Renderer','zbuffer');
%M(round(plotlength/dt)) = struct('cdata',[],'colormap',[]);
axis tight manual
%set(gca,'NextPlot','replaceChildren');
set(gca,'FontSize',14)
for k=1:plotlength/dt
    g=reshape(gg{k},I,J);
    acut = a(1:icut);
    gcut = g(1:icut,:);
    surf(acut,z,gcut')
    xlim([amin max(acut)])
    ylim([zmin zmax])
    zlim([0,0.5]);
    view([45 25])
    xlabel('Wealth, a','FontSize',14)
    ylabel('Productivity, z','FontSize',14)
    zlabel('Density g(a,z,t)','FontSize',14)
    title(sprintf('Evolution of Distribution (t=%0.1f)',k*dt));
    M(k)=getframe(gcf);
end
movie2avi(M(1:plotlength/dt),'distribution.avi','fps',50);
