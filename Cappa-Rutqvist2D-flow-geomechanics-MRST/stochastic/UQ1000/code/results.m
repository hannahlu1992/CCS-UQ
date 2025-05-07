clear all;
close all;

load("../data/QoI3.mat");
load("../data/miniMax100samples_percentile.mat");
load("../FML/Tmax/06tanh_logTmax_20mem_10rec_miniMax_100budget_1/Tmax_pred.mat");
load("../FMLc/Pmax/06tanh_logPmax_20mem_10rec_miniMax_100budget_1/Pmax_pred.mat");
load("../FML/W/06tanh_logW_20mem_10rec_miniMax_100budget_1/W_pred.mat");

for m = 1:1000
    Tmax(L06(m)+1:end,:,m) = nan;
    T_pred(L06(m)+1:end,:,m) = nan;
    dPmax(L06(m)+1:end,:,m) = nan;
    P_pred(L06(m)+1:end,:,m) = nan;
    dMw(L06(m)+1:end,:,m) = nan;
    W_pred(L06(m)+1:end,:,m) = nan;
end
T_pred_p10(max(L06):end) = nan;
T_pred_p50(max(L06):end) = nan;
T_pred_p90(max(L06):end) = nan;
P_pred_p10(max(L06):end) = nan;
P_pred_p50(max(L06):end) = nan;
P_pred_p90(max(L06):end) = nan;
W_pred_p10(max(L06):end) = nan;
W_pred_p50(max(L06):end) = nan;
W_pred_p90(max(L06):end) = nan;
T_p10(max(L06):end) = nan;
T_p50(max(L06):end) = nan;
T_p90(max(L06):end) = nan;
P_p10(max(L06):end) = nan;
P_p50(max(L06):end) = nan;
P_p90(max(L06):end) = nan;
W_p10(max(L06):end) = nan;
W_p50(max(L06):end) = nan;
W_p90(max(L06):end) = nan;
T100_p10(max(L06):end) = nan;
T100_p50(max(L06):end) = nan;
T100_p90(max(L06):end) = nan;
P100_p10(max(L06):end) = nan;
P100_p50(max(L06):end) = nan;
P100_p90(max(L06):end) = nan;
W100_p10(max(L06):end) = nan;
W100_p50(max(L06):end) = nan;
W100_p90(max(L06):end) = nan;


nbins = 15;
fh = figure();
width = 13;     % Width in inches
height = 2.5;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile(1)
a = min(min(Tmax(31,1,:)),min(T_pred(31,1,:)));
b = max(max(Tmax(31,1,:)),max(T_pred(31,1,:)));
edges = linspace(a, b, nbins);
histogram(Tmax(31,1,:),edges, 'Normalization', 'probability', ...
    'FaceColor', 'b', 'FaceAlpha', 0.5);
hold on;
histogram(T_pred(31,1,:),edges, 'Normalization', 'probability', ...
    'FaceColor', 'r', 'FaceAlpha', 0.5);
xlabel('$5$ days - $T_{s,max}$', 'Interpreter','latex','FontSize',12);
ylabel('P [-]', 'Interpreter','latex','FontSize',12);
xlim([0.2 0.6]);
ylim([0 0.2]);
grid on
% title('$5$ day', 'Interpreter','latex','FontSize',12);

nexttile(2)
a = min(min(dPmax(31,1,:)),min(P_pred(31,1,:)));
b = max(max(dPmax(31,1,:)),max(P_pred(31,1,:)));
edges = linspace(a/1e6, b/1e6, nbins);
histogram(dPmax(31,1,:)/1e6,edges, 'Normalization', 'probability', ...
    'FaceColor', 'b', 'FaceAlpha', 0.5);
hold on;
histogram(P_pred(31,1,:)/1e6,edges, 'Normalization', 'probability', ...
    'FaceColor', 'r', 'FaceAlpha', 0.5);
xlabel('$5$ days - $\Delta p_{max}$ [MPa]', 'Interpreter','latex','FontSize',12);
ylabel('P [-]', 'Interpreter','latex','FontSize',12);
xlim([0 5]);
grid on
set(gca, 'YScale', 'log')
ylim([10^(-3.5),1]);
% title('$5$ day', 'Interpreter','latex','FontSize',12,'Color',[1 1 1]);

nexttile(3)
a = min(min(dMw(31,1,:)),min(W_pred(31,1,:)));
b = max(max(dMw(31,1,:)),max(W_pred(31,1,:)));
edges = linspace(a/1e3, b/1e3, nbins);
histogram(dMw(31,1,:)/1e3,edges, 'Normalization', 'probability', ...
    'FaceColor', 'b', 'FaceAlpha', 0.5);
hold on;
histogram(W_pred(31,1,:)/1e3,edges, 'Normalization', 'probability', ...
    'FaceColor', 'r', 'FaceAlpha', 0.5);
xlabel('$5$ days - $\Delta M_{b}$ [tonne]', 'Interpreter','latex','FontSize',12);
ylabel('P [-]', 'Interpreter','latex','FontSize',12);
xlim([0 30]);
ylim([0 0.15]);
grid on
% title('$5$ day', 'Interpreter','latex','FontSize',12,'Color',[1 1 1]);
set(fh, 'position', [200, 200, 600, 200]);
exportgraphics(gcf, sprintf('../histogram_day5.pdf'));

nbins = 15;
fh = figure();
width = 13;     % Width in inches
height = 2.5;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile(1)
a = min(min(Tmax(61,1,:)),min(T_pred(61,1,:)));
b = max(max(Tmax(61,1,:)),max(T_pred(61,1,:)));
edges = linspace(a, b, nbins);
histogram(Tmax(61,1,:),edges, 'Normalization', 'probability', ...
    'FaceColor', 'b', 'FaceAlpha', 0.5);
hold on;
histogram(T_pred(61,1,:),edges, 'Normalization', 'probability', ...
    'FaceColor', 'r', 'FaceAlpha', 0.5);
xlabel('$10$ days - $T_{s,max}$', 'Interpreter','latex','FontSize',12);
ylabel('P [-]', 'Interpreter','latex','FontSize',12);
xlim([0.2 0.6]);
ylim([0 0.2]);
grid on
% title('$10$ days', 'Interpreter','latex','FontSize',12);

nexttile(2)
a = min(min(dPmax(61,1,:)),min(P_pred(61,1,:)));
b = max(max(dPmax(61,1,:)),max(P_pred(61,1,:)));
edges = linspace(a/1e6, b/1e6, nbins);
histogram(dPmax(61,1,:)/1e6,edges, 'Normalization', 'probability', ...
    'FaceColor', 'b', 'FaceAlpha', 0.5);
hold on;
histogram(P_pred(61,1,:)/1e6,edges, 'Normalization', 'probability', ...
    'FaceColor', 'r', 'FaceAlpha', 0.5);
xlabel('$10$ days - $\Delta p_{max}$ [MPa]', 'Interpreter','latex','FontSize',12);
ylabel('P [-]', 'Interpreter','latex','FontSize',12);
xlim([0 5]);
grid on
set(gca, 'YScale', 'log')
ylim([10^(-3.5),1]);
% title('$10$ days', 'Interpreter','latex','FontSize',12,'Color',[1 1 1]);

nexttile(3)
a = min(min(dMw(61,1,:)),min(W_pred(61,1,:)));
b = max(max(dMw(61,1,:)),max(W_pred(61,1,:)));
edges = linspace(a/1e3, b/1e3, nbins);
histogram(dMw(61,1,:)/1e3,edges, 'Normalization', 'probability', ...
    'FaceColor', 'b', 'FaceAlpha', 0.5);
hold on;
histogram(W_pred(61,1,:)/1e3,edges, 'Normalization', 'probability', ...
    'FaceColor', 'r', 'FaceAlpha', 0.5);
xlabel('$10$ days - $\Delta M_{b}$ [tonne]', 'Interpreter','latex','FontSize',12);
ylabel('P [-]', 'Interpreter','latex','FontSize',12);
xlim([0 30]);
ylim([0 0.15]);
grid on
% title('$10$ days', 'Interpreter','latex','FontSize',12,'Color',[1 1 1]);
set(fh, 'position', [200, 200, 600, 200]);
exportgraphics(gcf, sprintf('../histogram_day10.pdf'));

nbins = 15;
fh = figure();
width = 13;     % Width in inches
height = 2.5;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile(1)
a = min(min(Tmax(121,1,:)),min(T_pred(121,1,:)));
b = max(max(Tmax(121,1,:)),max(T_pred(121,1,:)));
edges = linspace(a, b, nbins);
histogram(Tmax(121,1,~isnan(Tmax(121,1,:))),edges, 'Normalization', 'probability', ...
    'FaceColor', 'b', 'FaceAlpha', 0.5);
hold on;
histogram(T_pred(121,1,~isnan(T_pred(121,1,:))),edges, 'Normalization', 'probability', ...
    'FaceColor', 'r', 'FaceAlpha', 0.5);
xlabel('$20$ days - $T_{s,max}$', 'Interpreter','latex','FontSize',12);
ylabel('P [-]', 'Interpreter','latex','FontSize',12);
xlim([0.2 0.6]);
ylim([0 0.2]);
grid on
% title('$20$ days', 'Interpreter','latex','FontSize',12);


nexttile(2)
a = min(min(dPmax(121,1,:)),min(P_pred(121,1,:)));
b = max(max(dPmax(121,1,:)),max(P_pred(121,1,:)));
edges = linspace(a/1e6, b/1e6, nbins);
histogram(dPmax(121,1,~isnan(dPmax(121,1,:)))/1e6,edges, 'Normalization', 'probability', ...
    'FaceColor', 'b', 'FaceAlpha', 0.5);
hold on;
histogram(P_pred(121,1,~isnan(P_pred(121,1,:)))/1e6,edges, 'Normalization', 'probability', ...
    'FaceColor', 'r', 'FaceAlpha', 0.5);
xlabel('$20$ days - $\Delta p_{max}$ [MPa]', 'Interpreter','latex','FontSize',12);
ylabel('P [-]', 'Interpreter','latex','FontSize',12);
grid on
xlim([0 5]);
set(gca, 'YScale', 'log')
ylim([10^(-3.5),1]);
% title('$20$ days', 'Interpreter','latex','FontSize',12,'Color',[1 1 1]);


nexttile(3)
a = min(min(dMw(121,1,:)),min(W_pred(121,1,:)));
b = max(max(dMw(121,1,:)),max(W_pred(121,1,:)));
edges = linspace(a/1e3, b/1e3, nbins);
histogram(dMw(121,1,~isnan(dMw(121,1,:)))/1e3,edges, 'Normalization', 'probability', ...
    'FaceColor', 'b', 'FaceAlpha', 0.5);
hold on;
histogram(W_pred(121,1,~isnan(W_pred(121,1,:)))/1e3,edges, 'Normalization', 'probability', ...
    'FaceColor', 'r', 'FaceAlpha', 0.5);
xlabel('$20$ days - $\Delta M_{b}$ [tonne]', 'Interpreter','latex','FontSize',12);
ylabel('P [-]', 'Interpreter','latex','FontSize',12);
xlim([0 30]);
ylim([0 0.15]);
grid on
% title('$20$ days', 'Interpreter','latex','FontSize',12,'Color',[1 1 1]);
set(fh, 'position', [200, 200, 600, 200]);
exportgraphics(gcf, sprintf('../histogram_day20.pdf'));


t = 1:4:800;
fh = figure();
width = 13;     % Width in inches
height = 2.5;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile(1)
hold on;
plot(t,T_p10/1e3,'k-','LineWidth',1);
plot(t,T100_p10/1e3,'b.','MarkerSize',4);
plot(t,T_pred_p10/1e3,'r--','LineWidth',1);
plot(t,T_p50/1e3,'k-','LineWidth',1);
plot(t,T100_p50/1e3,'b.','MarkerSize',4);
plot(t,T_pred_p50/1e3,'r--','LineWidth',1);
plot(t,T_p90/1e3,'k-','LineWidth',1);
plot(t,T100_p90/1e3,'b.','MarkerSize',4);
plot(t,T_pred_p90/1e3,'r--','LineWidth',1);
xlabel('$t$ [hour]', 'Interpreter','latex','FontSize',12);
ylabel('$T_{s,max}$', 'Interpreter','latex','FontSize',12);
xlim([0,t(max(L06))]);


nexttile(2)
hold on;
plot(t,P_p10/1e6,'k-','LineWidth',1);
plot(t,P100_p10/1e6,'b.','MarkerSize',4);
plot(t,P_pred_p10/1e6,'r--','LineWidth',1);
plot(t,P_p50/1e6,'k-','LineWidth',1);
plot(t,P100_p50/1e6,'b.','MarkerSize',4);
plot(t,P_pred_p50/1e6,'r--','LineWidth',1);
plot(t,P_p90/1e6,'k-','LineWidth',1);
plot(t,P100_p90/1e6,'b.','MarkerSize',4);
plot(t,P_pred_p90/1e6,'r--','LineWidth',1);
xlabel('$t$ [hour]', 'Interpreter','latex','FontSize',12);
ylabel('$\Delta p_{max}$ [MPa]', 'Interpreter','latex','FontSize',12);
xlim([0,t(max(L06))]);

nexttile(3)
hold on;
plot(t,W_p10/1e3,'k-','LineWidth',1);
plot(t,W100_p10/1e3,'b.','MarkerSize',4);
plot(t,W_pred_p10/1e3,'r--','LineWidth',1);
plot(t,W_p50/1e3,'k-','LineWidth',1);
plot(t,W100_p50/1e3,'b.','MarkerSize',4);
plot(t,W_pred_p50/1e3,'r--','LineWidth',1);
plot(t,W_p90/1e3,'k-','LineWidth',1);
plot(t,W100_p90/1e3,'b.','MarkerSize',4);
plot(t,W_pred_p90/1e3,'r--','LineWidth',1);
xlabel('$t$ [hour]', 'Interpreter','latex','FontSize',12);
ylabel('$\Delta M_{b}$ [tonne]', 'Interpreter','latex','FontSize',12);
xlim([0,t(max(L06))]);

set(fh, 'position', [200, 200, 600, 200]);
exportgraphics(gcf, sprintf('../percentile.pdf'));
