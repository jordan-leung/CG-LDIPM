clc
clear all
close all

% Load
load('data_10p0')
colorMatrix = [0 0 0; 228 26 28; 55 126 184; 77 175 74]./255;



% Loop through num iter to see if there are outliers
NUM_ITER_FILTER = NUM_ITER;
% for i = 1:NSample
%     for j = 1:NSample
%         if NUM_ITER(i,j) > 10000
%             NUM_ITER_FILTER(i,j) = 10000;
%         end
%     end
% end


% Plot the surf plot
figure
figSize = [0 0 0.25 0.3]*.75;
set(gcf,'units','normalized','position',figSize)
h1 = surf(BETA,YY,NUM_ITER_FILTER');
set(gca,'zscale','log')
ZLIM = zlim(gca); % query zlim 
xlim([xmin(1), xmax(1)])
ylim([xmin(2), xmax(2)])
hold on

% Get the vertices of gamma, add the Z component, and plot as a polyhedron
% in 3D
Gamma_vert = con2vert(GammaSlice_A,GammaSlice_b);
VV = [Gamma_vert, ZLIM(1)*ones(size(Gamma_vert,1),1)]';
rep.V = VV;
P = polyh(rep,'v');
OPT.color = colorMatrix(2,:);
OPT.edgewidth = 2;
OPT.vertexsize = 0.0001;
plot(P,OPT);

% Do the same with O_infty
Oinf_vert = con2vert(OinfSlice_A, OinfSlice_b);
VV = [Oinf_vert, 1.0001*ZLIM(1)*ones(size(Gamma_vert,1),1)]';
rep.V = VV;
P = polyh(rep,'v');
OPT.color = colorMatrix(4,:);
OPT.edgewidth = 2;
OPT.vertexsize = 0.0001;
plot(P,OPT);
xlabel('$\beta$','interpreter','latex','fontsize',15)
ylabel('$y$','interpreter','latex','fontsize',15)
zlabel('Number of iterations','interpreter','latex','fontsize',15)
% legend('Number of iterations','$\Gamma_N$','','','$O_\infty$','interpreter','latex','fontsize',12,'location','northeast')
hf = gcf;
lines = hf.CurrentAxes.Children;                                    % Get °line' Objects From Structure
% legend([lines(end-1) lines(3)],'$\Gamma_N$','$O_\infty$',...
%     'Location','best','interpreter','latex','fontsize',12)    % Create 'legend' Object
load('cam')
campos(cam)
% title('$\epsilon = 0.1$','interpreter','latex','fontsize',19)
title('$\epsilon = 10$','interpreter','latex','fontsize',19)


% % Plot the initial conditions for which feasibile solutions were found on
% % top of figure 1
% % h1 = GammaSlice.plot('Alpha',.2,'LineWidth',1,'EdgeColor',colorMatrix(2,:));
% xlabel('$\beta$')
% ylabel('$y$')
% hold on 
% % OinfSlice.plot('Color',colorMatrix(4,:),'LineWidth',1,'EdgeColor',colorMatrix(4,:),'alpha',.4);
% for i = 1:NSample
%     for j = 1:NSample
%         if ~isnan(NUM_ITER(i,j))
%             plot3(betaSample(i),ySample(j),100,'k.')
%         end
%     end
% end
