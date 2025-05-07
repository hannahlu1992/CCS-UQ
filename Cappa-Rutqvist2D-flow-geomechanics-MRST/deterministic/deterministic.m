clear all;
close all;

addpath("../../mrst-developer");
addpath("../");
run startup.m

warning('off','all')
mrstModule add diagnostics ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui 

%% write intro text for each case
writeIntroText = @(opt)(fprintf('\n*** Start new simulation\n* fluid model : %s\n* method : %s\n* injection rate : %.2f\n* fault permeability : %.2f\n* bc ratio : %.2f\n* reservoir permeability : %.2f\n* reservoir rock compressibility : %s\n* reservoir porosity : %.2f\n\n', opt.fluid_model, opt.method, opt.inj_rate,log10(opt.kf(1)), opt.lambda,log10(opt.k_reservoir),opt.reservoir_compressibility,opt.poro_reservoir));


%% Summary of the options
clear opt
opt.nonlinearTolerance = 1e-6;
opt.splittingTolerance = 1e-6;
opt.verbose            = false;
opt.splittingVerbose   = false;
opt.mech_BC = 'left bottom roller';
opt.inj_rate = 0.02;
opt.fluid_model = 'blackoil-ccs';
opt.method = 'fully coupled';
opt.reportTimes = [(1:4:24)*hour,(2:1:10)*day];
opt.lambda = 0.7;
opt.k_reservoir = 1e-13;
opt.poro_reservoir = 0.1;
opt.reservoir_compressibility = 1.45e-4;
opt.kf = 1e-16; 
opt.initialization = [];

writeIntroText(opt);
optvals = cellfun(@(x) opt.(x), fieldnames(opt), 'uniformoutput', false);
optlist = reshape(vertcat(fieldnames(opt)', optvals'), [], 1);
[model, states,initState,reservoir_indx,fault_indx_collect] = runCappaRutqvist2D(optlist{:});

S = cell(numel(states)+1,1);
S{1} = initState;
S(2:end) = states;
[P,T,tau,sigma_n_prime] = computingQoIs(states,fault_indx_collect,80);

figure
width = 6.5;     % Width in inches
height = 5;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

plotCellData(model.G, (states{end}.pressure-initState.pressure)/1e6, 'edgealpha', 0);
colormap(jet), c = colorbar;
c.Label.Interpreter = 'latex'; c.Label.FontSize = 20;
c.Label.String = '$\Delta p$ [MPa]';
set(gca, 'YDir','reverse');
xlim([0,1000]);
ylim([1000,2000]);
xticks([0,1000]);
yticks([1000,2000]);
xlabel('$x[m]$','Interpreter','latex','FontSize',20);
ylabel('$z[m]$','Interpreter','latex','FontSize',20);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',20)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',20)

exportgraphics(gcf, 'dp10days.pdf', 'ContentType', 'vector');

figure
width = 6.5;     % Width in inches
height = 5;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

plotCellData(model.G, states{end}.s(:,3), 'edgealpha', 0);
colormap(jet), c = colorbar;
c.Label.Interpreter = 'latex'; c.Label.FontSize = 20;
c.Label.String = '$S_{CO_2}$';
set(gca, 'YDir','reverse');
xlim([0,1000]);
ylim([1000,2000]);
xticks([0,1000]);
yticks([1000,2000]);
xlabel('$x[m]$','Interpreter','latex','FontSize',20);
ylabel('$z[m]$','Interpreter','latex','FontSize',20);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',20)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',20)
exportgraphics(gcf, 'S10days.pdf', 'ContentType', 'vector');

figure
width = 6.5;     % Width in inches
height = 5;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

plotNodeData(model.G, (states{end}.uu(:,1)-initState.uu(:,1)), 'edgealpha', 0);
colormap(jet), c = colorbar;
c.Label.Interpreter = 'latex'; c.Label.FontSize = 20;
c.Label.String = '$\Delta u_x$ [m]';
set(gca, 'YDir','reverse');
xlim([0,1000]);
ylim([1000,2000]);
xticks([0,1000]);
yticks([1000,2000]);
xlabel('$x[m]$','Interpreter','latex','FontSize',20);
ylabel('$z[m]$','Interpreter','latex','FontSize',20);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',20)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',20)

exportgraphics(gcf, 'dux10days.pdf', 'ContentType', 'vector');

figure
width = 6.5;     % Width in inches
height = 5;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

plotNodeData(model.G, (states{end}.uu(:,2)-initState.uu(:,2)), 'edgealpha', 0);
colormap(jet), c = colorbar;
c.Label.Interpreter = 'latex'; c.Label.FontSize = 20;
c.Label.String = '$\Delta u_z$ [m]';
set(gca, 'YDir','reverse');
xlim([0,1000]);
ylim([1000,2000]);
xticks([0,1000]);
yticks([1000,2000]);
xlabel('$x[m]$','Interpreter','latex','FontSize',20);
ylabel('$z[m]$','Interpreter','latex','FontSize',20);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',20)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',20)

exportgraphics(gcf, 'duz10days.pdf', 'ContentType', 'vector');



figure
width = 6.5;     % Width in inches
height = 5;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

hold on;
t = tiledlayout(1,1);
ax1 = axes(t);
hold on;
plot(ax1,P(6,:)/1e6,500+1.25:2.5:2500-1.25,'b-','LineWidth',2);
plot(ax1,P(9,:)/1e6,500+1.25:2.5:2500-1.25,'k-.','LineWidth',2);
plot(ax1,P(14,:)/1e6,500+1.25:2.5:2500-1.25,'r--','LineWidth',2);
xlim([-0.03,4])
xlabel({'$\Delta p$ [MPa]'},'FontUnits','points','interpreter','latex',...
'FontSize',20);
ylabel({'$z$ [m]'},'FontUnits','points','interpreter','latex',...
'FontSize',20);
 legend({'$2$ days','$5$ days','$10$ days'},'interpreter','latex',...
'FontSize',20,'location','best','box','off');

ax2 = axes(t);
hold on;
v1 = [0 1300;0 1450;508.2654 1450;535.2654 1300];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);
v2 = [0 1700; 0 1550; 491.1837 1550; 464.736 1700];
f2 = [1 2 3 4];
patch('Faces',f2,'Vertices',v2,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);
v3 = [486.7755 1575;513.2245 1425; 2000 1425; 2000 1575];
f3 = [1 2 3 4];
patch('Faces',f3,'Vertices',v3,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);
v4 = [530.8572 1325;557.3063 1175;2000 1175; 2000 1325];
f4 = [1 2 3 4];
patch('Faces',f4,'Vertices',v4,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);
v5 = [323.6730-1.25 2500;676.3270-1.25 500; 676.3270+1.25 500; 323.6730+1.25 2500];
f5 = [1 2 3 4];
patch('Faces',f5,'Vertices',v5,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);


ax1.Box = 'off';
ax1.YDir = "reverse";
ax1.FontSize = 20;
ax1.XAxisLocation = 'bottom';
ax1.YAxisLocation = 'left';


ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax2.Box = 'off';
ax2.YDir = "reverse";
ax2.XLim = [0,2000];
ax2.XTick = [];
ax2.YTick = [];
exportgraphics(gcf, 'fault_dp10days.pdf', 'ContentType', 'vector');


figure
width = 6.5;     % Width in inches
height = 5;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

hold on;
t = tiledlayout(1,1);
ax1 = axes(t);
hold on;
plot(ax1,(tau(7,:)-tau(1,:))/1e6,500+1.25:2.5:2500-1.25,'b-','LineWidth',2);
plot(ax1,(tau(10,:)-tau(1,:))/1e6,500+1.25:2.5:2500-1.25,'k-.','LineWidth',2);
plot(ax1,(tau(15,:)-tau(1,:))/1e6,500+1.25:2.5:2500-1.25,'r--','LineWidth',2);
xlabel({'$\Delta \tau$ [MPa]'},'FontUnits','points','interpreter','latex',...
'FontSize',20);
ylabel({'$z$ [m]'},'FontUnits','points','interpreter','latex',...
'FontSize',20);
 legend({'$2$ days','$5$ days','$10$ days'},'interpreter','latex',...
'FontSize',20,'location','best','box','off');

ax2 = axes(t);
hold on;
v1 = [0 1300;0 1450;508.2654 1450;535.2654 1300];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);
v2 = [0 1700; 0 1550; 491.1837 1550; 464.736 1700];
f2 = [1 2 3 4];
patch('Faces',f2,'Vertices',v2,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);
v3 = [486.7755 1575;513.2245 1425; 2000 1425; 2000 1575];
f3 = [1 2 3 4];
patch('Faces',f3,'Vertices',v3,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);
v4 = [530.8572 1325;557.3063 1175;2000 1175; 2000 1325];
f4 = [1 2 3 4];
patch('Faces',f4,'Vertices',v4,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);
v5 = [323.6730-1.25 2500;676.3270-1.25 500; 676.3270+1.25 500; 323.6730+1.25 2500];
f5 = [1 2 3 4];
patch('Faces',f5,'Vertices',v5,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);


ax1.Box = 'off';
ax1.YDir = "reverse";
ax1.FontSize = 20;
ax1.XAxisLocation = 'bottom';
ax1.YAxisLocation = 'left';


ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax2.Box = 'off';
ax2.YDir = "reverse";
ax2.XLim = [0,2000];
ax2.XTick = [];
ax2.YTick = [];
exportgraphics(gcf, 'fault_tau10days.pdf', 'ContentType', 'vector');

figure
width = 6.5;     % Width in inches
height = 5;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

hold on;
t = tiledlayout(1,1);
ax1 = axes(t);
hold on;
plot(ax1,(sigma_n_prime(7,:)-sigma_n_prime(1,:))/1e6,500+1.25:2.5:2500-1.25,'b-','LineWidth',2);
plot(ax1,(sigma_n_prime(10,:)-sigma_n_prime(1,:))/1e6,500+1.25:2.5:2500-1.25,'k-.','LineWidth',2);
plot(ax1,(sigma_n_prime(15,:)-sigma_n_prime(1,:))/1e6,500+1.25:2.5:2500-1.25,'r--','LineWidth',2);
xlabel({'$\Delta \sigma_n^\prime$ [MPa]'},'FontUnits','points','interpreter','latex',...
'FontSize',20);
ylabel({'$z$ [m]'},'FontUnits','points','interpreter','latex',...
'FontSize',20);
 legend({'$2$ days','$5$ days','$10$ days'},'interpreter','latex',...
'FontSize',20,'location','best','box','off');

ax2 = axes(t);
hold on;
v1 = [0 1300;0 1450;508.2654 1450;535.2654 1300];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);
v2 = [0 1700; 0 1550; 491.1837 1550; 464.736 1700];
f2 = [1 2 3 4];
patch('Faces',f2,'Vertices',v2,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);
v3 = [486.7755 1575;513.2245 1425; 2000 1425; 2000 1575];
f3 = [1 2 3 4];
patch('Faces',f3,'Vertices',v3,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);
v4 = [530.8572 1325;557.3063 1175;2000 1175; 2000 1325];
f4 = [1 2 3 4];
patch('Faces',f4,'Vertices',v4,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);
v5 = [323.6730-1.25 2500;676.3270-1.25 500; 676.3270+1.25 500; 323.6730+1.25 2500];
f5 = [1 2 3 4];
patch('Faces',f5,'Vertices',v5,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);


ax1.Box = 'off';
ax1.YDir = "reverse";
ax1.FontSize = 20;
ax1.XAxisLocation = 'bottom';
ax1.YAxisLocation = 'left';


ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax2.Box = 'off';
ax2.YDir = "reverse";
ax2.XLim = [0,2000];
ax2.XTick = [];
ax2.YTick = [];
exportgraphics(gcf, 'fault_sigma_n10days.pdf', 'ContentType', 'vector');

figure
width = 6.5;     % Width in inches
height = 5;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

hold on;
t = tiledlayout(1,1);
ax1 = axes(t);
hold on;
plot(ax1,T(7,:),500+1.25:2.5:2500-1.25,'b-','LineWidth',2);
plot(ax1,T(10,:),500+1.25:2.5:2500-1.25,'k-.','LineWidth',2);
plot(ax1,T(15,:),500+1.25:2.5:2500-1.25,'r--','LineWidth',2);
xlabel({'$T_s$'},'FontUnits','points','interpreter','latex',...
'FontSize',20);
ylabel({'$z$ [m]'},'FontUnits','points','interpreter','latex',...
'FontSize',20);
 legend({'$2$ days','$5$ days','$10$ days'},'interpreter','latex',...
'FontSize',20,'location','best','box','off');

ax2 = axes(t);
hold on;
v1 = [0 1300;0 1450;508.2654 1450;535.2654 1300];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);
v2 = [0 1700; 0 1550; 491.1837 1550; 464.736 1700];
f2 = [1 2 3 4];
patch('Faces',f2,'Vertices',v2,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);
v3 = [486.7755 1575;513.2245 1425; 2000 1425; 2000 1575];
f3 = [1 2 3 4];
patch('Faces',f3,'Vertices',v3,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);
v4 = [530.8572 1325;557.3063 1175;2000 1175; 2000 1325];
f4 = [1 2 3 4];
patch('Faces',f4,'Vertices',v4,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);
v5 = [323.6730-1.25 2500;676.3270-1.25 500; 676.3270+1.25 500; 323.6730+1.25 2500];
f5 = [1 2 3 4];
patch('Faces',f5,'Vertices',v5,'FaceColor',[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',.3);


ax1.Box = 'off';
ax1.YDir = "reverse";
ax1.FontSize = 20;
ax1.XAxisLocation = 'bottom';
ax1.YAxisLocation = 'left';


ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax2.Box = 'off';
ax2.YDir = "reverse";
ax2.XLim = [0,2000];
ax2.XTick = [];
ax2.YTick = [];
exportgraphics(gcf, 'Ts10days.pdf', 'ContentType', 'vector');





