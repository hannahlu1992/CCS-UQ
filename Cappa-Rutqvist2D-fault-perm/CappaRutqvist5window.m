%% Cappa-Rutqvist Model
% 5 throw windows of depth 120m.

clear
close all force
%% 1. Load Required MRST Modules and PREDICT Folder
% First, navigate to the mrst folder and run |startup.m|. We can then load the 
% appropriate modules for generating MRST grids and upscale the permeability:
addpath("../mrst-developer");
run startup.m
mrstModule add mrst-gui coarsegrid upscaling incomp mpfa
mrstVerbose off
% Navigate to PREDICT folder
addpath(genpath("../predict-main/"));

%% 2. Define Throw Windows
clay_fraction_clay = 0.4;
clay_fraction_sand = 0.2;
Nthrow = 5;
Nsim = 1000;       % Number of simulations/realizations

%% throw window 1
faultDip  = 80; 
thickness = {[125]/cosd(90-faultDip), [125]/cosd(90-faultDip)};                % [m]
vcl       = {[clay_fraction_clay], ...
             [clay_fraction_sand]};                        % fraction [-]
dip       = [0, 0];                                                        % [deg.]                                                            % [deg.]
zf        = [200, 200];                                                     % [FW, HW], [m]
zmax      = {repelem(1150, numel(vcl{1})), repelem(1150, numel(vcl{2}))};   % {FW, HW}
dim       = 2;                    % dimensions (2 = 2D, 3 = 3D)
unit_plot = 'm';                  % 'm' or 'cm' depending on fault dimensions
totalThickness = [150, ...    % [FW
              nan]/cosd(90-faultDip);  % HW] same criterion for Tap or trueT as thickness above

throw1_footwall = Stratigraphy(thickness{1}, vcl{1}, 'Dip', dip(1), ...
                        'DepthFaulting', zf(1), 'DepthBurial', zmax{1});
throw1_hangingwall = Stratigraphy(thickness{2}, vcl{2}, 'Dip', dip(2), 'IsHW', 1, ...
                           'NumLayersFW', throw1_footwall.NumLayers, ...
                           'DepthFaulting', zf(2), 'DepthBurial', zmax{2});
[throw1_mySect,throw1_faults,throw1_smears,throw1_G0] = singleStrati(throw1_footwall,throw1_hangingwall,faultDip,totalThickness,Nsim);
throw1_thick = cellfun(@(x) x.MatProps.thick, throw1_faults);

%% throw window 2
faultDip  = 80; 
thickness = {[100,25]/cosd(90-faultDip), [125]/cosd(90-faultDip)};                % [m]
vcl       = {[clay_fraction_sand,clay_fraction_clay], ...
             [clay_fraction_clay]};                        % fraction [-]
dip       = [0, 0];                                                        % [deg.]                                                            % [deg.]
zf        = [200, 200];                                                     % [FW, HW], [m]
zmax      = {repelem(1275, numel(vcl{1})), repelem(1275, numel(vcl{2}))};   % {FW, HW}
dim       = 2;                    % dimensions (2 = 2D, 3 = 3D)
unit_plot = 'm';                  % 'm' or 'cm' depending on fault dimensions
totalThickness = [nan, 150, 150]/cosd(90-faultDip);  % HW] same criterion for Tap or trueT as thickness above
throw2_footwall = Stratigraphy(thickness{1}, vcl{1}, 'Dip', dip(1), ...
                        'DepthFaulting', zf(1), 'DepthBurial', zmax{1});
throw2_hangingwall = Stratigraphy(thickness{2}, vcl{2}, 'Dip', dip(2), 'IsHW', 1, ...
                           'NumLayersFW', throw2_footwall.NumLayers, ...
                           'DepthFaulting', zf(2), 'DepthBurial', zmax{2});
[throw2_mySect,throw2_faults,throw2_smears,throw2_G0] = singleStrati(throw2_footwall,throw2_hangingwall,faultDip,totalThickness,Nsim,throw1_thick);

%% throw window 3
faultDip  = 80;      
thickness = {[125]/cosd(90-faultDip), [100,25]/cosd(90-faultDip)};                % [m]
vcl       = {[clay_fraction_clay], ...
             [clay_fraction_sand,clay_fraction_clay]};                        % fraction [-]
dip       = [0, 0];                                                        % [deg.]                                          % [deg.]
zf        = [200, 200];                                                     % [FW, HW], [m]
zmax      = {repelem(1400, numel(vcl{1})), repelem(1400, numel(vcl{2}))};   % {FW, HW}
dim       = 2;                    % dimensions (2 = 2D, 3 = 3D)
unit_plot = 'm';                  % 'm' or 'cm' depending on fault dimensions
totalThickness = [150, ...    % [FW
              nan, 150]/cosd(90-faultDip);  % HW] same criterion for Tap or trueT as thickness above
throw3_footwall = Stratigraphy(thickness{1}, vcl{1}, 'Dip', dip(1), ...
                        'DepthFaulting', zf(1), 'DepthBurial', zmax{1});
throw3_hangingwall = Stratigraphy(thickness{2}, vcl{2}, 'Dip', dip(2), 'IsHW', 1, ...
                           'NumLayersFW', throw3_footwall.NumLayers, ...
                           'DepthFaulting', zf(2), 'DepthBurial', zmax{2});
[throw3_mySect,throw3_faults,throw3_smears,throw3_G0] = singleStrati(throw3_footwall,throw3_hangingwall,faultDip,totalThickness,Nsim,throw1_thick);

%% throw window 4
faultDip  = 80; 
thickness = {[100,25]/cosd(90-faultDip), [125]/cosd(90-faultDip)};                % [m]
vcl       = {[clay_fraction_sand,clay_fraction_clay], ...
             [clay_fraction_clay]};                        % fraction [-]
dip       = [0, 0];                                                        % [deg.]                                       % [deg.]
zf        = [200, 200];                                                     % [FW, HW], [m]
zmax      = {repelem(1525, numel(vcl{1})), repelem(1525, numel(vcl{2}))};   % {FW, HW}
dim       = 2;                    % dimensions (2 = 2D, 3 = 3D)
unit_plot = 'm';                  % 'm' or 'cm' depending on fault dimensions
totalThickness = [nan,150, ...    % [FW
              150]/cosd(90-faultDip);  % HW] same criterion for Tap or trueT as thickness above
throw4_footwall = Stratigraphy(thickness{1}, vcl{1}, 'Dip', dip(1), ...
                        'DepthFaulting', zf(1), 'DepthBurial', zmax{1});
throw4_hangingwall = Stratigraphy(thickness{2}, vcl{2}, 'Dip', dip(2), 'IsHW', 1, ...
                           'NumLayersFW', throw4_footwall.NumLayers, ...
                           'DepthFaulting', zf(2), 'DepthBurial', zmax{2});
[throw4_mySect,throw4_faults,throw4_smears,throw4_G0] = singleStrati(throw4_footwall,throw4_hangingwall,faultDip,totalThickness,Nsim,throw1_thick);

%% throw window 5
faultDip  = 80;  
thickness = {[125]/cosd(90-faultDip), [100,25]/cosd(90-faultDip)};                % [m]
vcl       = {[clay_fraction_sand], ...
             [clay_fraction_sand,clay_fraction_clay]};                        % fraction [-]
dip       = [0, 0];                                                        % [deg.]                                    % [deg.]
zf        = [200, 200];                                                     % [FW, HW], [m]
zmax      = {repelem(1650, numel(vcl{1})), repelem(1650, numel(vcl{2}))};   % {FW, HW}
dim       = 2;                    % dimensions (2 = 2D, 3 = 3D)
unit_plot = 'm';                  % 'm' or 'cm' depending on fault dimensions
totalThickness = [nan, ...    % [FW
              nan, 150]/cosd(90-faultDip);  % HW] same criterion for Tap or trueT as thickness above
throw5_footwall = Stratigraphy(thickness{1}, vcl{1}, 'Dip', dip(1), ...
                        'DepthFaulting', zf(1), 'DepthBurial', zmax{1});
throw5_hangingwall = Stratigraphy(thickness{2}, vcl{2}, 'Dip', dip(2), 'IsHW', 1, ...
                           'NumLayersFW', throw5_footwall.NumLayers, ...
                           'DepthFaulting', zf(2), 'DepthBurial', zmax{2});
[throw5_mySect,throw5_faults,throw5_smears,throw5_G0] = singleStrati(throw5_footwall,throw5_hangingwall,faultDip,totalThickness,Nsim,throw1_thick);


% %% 3. Output Analysis
% outputs_analysis(throw1_mySect,throw1_faults,throw1_smears,faultDip,unit_plot,1,dim,Nsim,throw1_G0);
% outputs_analysis(throw2_mySect,throw2_faults,throw2_smears,faultDip,unit_plot,1,dim,Nsim,throw2_G0);
% outputs_analysis(throw3_mySect,throw3_faults,throw3_smears,faultDip,unit_plot,1,dim,Nsim,throw3_G0);
% outputs_analysis(throw4_mySect,throw4_faults,throw4_smears,faultDip,unit_plot,1,dim,Nsim,throw4_G0);
% outputs_analysis(throw5_mySect,throw5_faults,throw5_smears,faultDip,unit_plot,1,dim,Nsim,throw5_G0);

% %% Save the output
% save(sprintf('%dsamples.mat',Nsim),'-v7.3');

%% 4. Compute upscaled fault permeability
%% Basic options for upscaling
U.useAcceleration = 1;          % 1 requires MEX setup, 0 otherwise (slower for MPFA).
U.method          = 'mpfa';     % 'tpfa' recommended if useAcceleration = 0
U.outflux         = 0;          % compare outflux of fine and upscaled model
U.ARcheck         = 0;          % check if Perm obtained with grid with aspect ratio of 
                                % only 5 gives same output.
U.coarseDims      = [1 1 1];    % Mandatory one cell if 2D

kf_ups_5 = zeros(3,5,Nsim);
for n = 1:Nsim
    % Generate grid
    G = computeGeometry(cartGrid([100, 100], [throw1_thick(n), throw1_faults{1}.Disp]));
    % Upscale permeability 
    permy = throw1_faults{1}.Perm(2); % needs updating if you actually want to use kyy
    kf_ups_5(:,1,n) = computeCoarsePerm(G, throw1_faults{n}.Grid.perm, permy, U)';
    kf_ups_5(:,2,n) = computeCoarsePerm(G, throw2_faults{n}.Grid.perm, permy, U)';
    kf_ups_5(:,3,n) = computeCoarsePerm(G, throw3_faults{n}.Grid.perm, permy, U)';
    kf_ups_5(:,4,n) = computeCoarsePerm(G, throw4_faults{n}.Grid.perm, permy, U)';
    kf_ups_5(:,5,n) = computeCoarsePerm(G, throw5_faults{n}.Grid.perm, permy, U)';
end
save(sprintf('%dkf.mat',Nsim), "kf_ups_5","-v7.3");


%% 5. Visualization
figure(1)
width = 6.5;     % Width in inches
height = 5;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);


n = randi(Nsim, 1, 1);
kf = [throw5_faults{n}.Grid.perm;throw4_faults{n}.Grid.perm;...
    throw3_faults{n}.Grid.perm;throw2_faults{n}.Grid.perm;throw1_faults{n}.Grid.perm];
G = computeGeometry(cartGrid([100, 100*Nthrow], [throw1_thick(n), throw1_faults{1}.Disp*Nthrow]));
subplot(1,2,1)
plotCellData(G,log10(kf(:,1)/(milli*darcy)),'EdgeAlpha', 0.1); colormap(copper);
plotCellData(G,log10(kf(:,1)/(milli*darcy)),'EdgeAlpha', 0.1); colormap(copper);
c = colorbar;c.Label.FontSize = 20;
c.FontSize = 20;
hold on;
plot([0,throw1_thick(n)],[throw1_mySect.Tap(1),throw1_mySect.Tap(1)],'r-','LineWidth',2);
plot([0,throw1_thick(n)],[2*throw1_mySect.Tap(1),2*throw1_mySect.Tap(1)],'r-','LineWidth',2);
plot([0,throw1_thick(n)],[3*throw1_mySect.Tap(1),3*throw1_mySect.Tap(1)],'r-','LineWidth',2);
plot([0,throw1_thick(n)],[4*throw1_mySect.Tap(1),4*throw1_mySect.Tap(1)],'r-','LineWidth',2);
xlim([0,throw1_thick(n)]);
ylim([0,5*throw1_mySect.Tap(1)]);
xticks([]);
yticks([]);
title('$\log_{10}(k_{xx})$ [mD]','Interpreter','latex','FontSize',20);
pbaspect([0.1 1 1]);

subplot(1,2,2)
plotCellData(G,log10(kf(:,3)/(milli*darcy)),'EdgeAlpha', 0.1); colormap(copper);
c = colorbar;c.Label.FontSize = 20;c.FontSize = 20;
hold on;
plot([0,throw1_thick(n)],[throw1_mySect.Tap(1),throw1_mySect.Tap(1)],'r-','LineWidth',2);
plot([0,throw1_thick(n)],[2*throw1_mySect.Tap(1),2*throw1_mySect.Tap(1)],'r-','LineWidth',2);
plot([0,throw1_thick(n)],[3*throw1_mySect.Tap(1),3*throw1_mySect.Tap(1)],'r-','LineWidth',2);
plot([0,throw1_thick(n)],[4*throw1_mySect.Tap(1),4*throw1_mySect.Tap(1)],'r-','LineWidth',2);
xlim([0,throw1_thick(n)]);
ylim([0,5*throw1_mySect.Tap(1)]);
xticks([]);
yticks([]);
title('$\log_{10}(k_{zz})$ [mD]','Interpreter','latex','FontSize',20);
pbaspect([0.1 1 1]);
exportgraphics(gcf, 'placement.pdf');


for k = 1:5 %thorw window
    fh = figure(k);
    width = 13;     % Width in inches
    height = 2.5;    % Height in inches
    set(gcf,'InvertHardcopy','on');
    set(gcf,'PaperUnits', 'inches');
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(gcf,'PaperPosition', myfiguresize);

    [R, P] = corrcoef(log10(reshape(kf_ups_5(:,k,:),[3,Nsim])'/(milli*darcy)));           % corrcoeff and pval matrices
    r = R'; r = r(3);
    nbins = 25;
    logMinP = -6; %min(min(K));
    logMaxP = 2; %max(max(K));
    edges = linspace(fix(logMinP)-1, fix(logMaxP)+1, nbins);

    tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile(1)
    histogram(log10(kf_ups_5(1,k,:)/(milli*darcy)), edges, 'Normalization', 'probability', ...
        'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1);
    xlabel('$\log_{10}(k_{xx})$ [mD]', 'Interpreter','latex','FontSize',12);
    ylabel('P [-]', 'Interpreter','latex','FontSize',12);
    xlim([fix(logMinP)-1 fix(logMaxP)+1])
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12)
    a = get(gca,'YTickLabel');
    set(gca,'YTickLabel',a,'fontsize',12)
    grid on

    nexttile(2)
    histogram(log10(kf_ups_5(3,k,:)/(milli*darcy)), edges, 'Normalization', 'probability', ...
        'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1);
    xlabel('$\log_{10}(k_{zz})$ [mD]', 'Interpreter','latex','FontSize',12);
    ylabel('P [-]', 'Interpreter','latex','FontSize',12);
    xlim([fix(logMinP)-1 fix(logMaxP)+1])
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12)
    a = get(gca,'YTickLabel');
    set(gca,'YTickLabel',a,'fontsize',12)
    grid on

    nexttile(3)
    colormap(turbo);
    histogram2(log10(kf_ups_5(1,k,:)/(milli*darcy)), log10(kf_ups_5(3,k,:)/(milli*darcy)), edges, edges, ...
    'Normalization', 'Probability', 'DisplayStyle','tile', ...
    'ShowEmptyBins','off');
    c = colorbar;
    c.Label.Interpreter = 'latex';
    c.Label.String = 'P [-]';
    c.Label.FontSize = 12;
    c.FontSize = 12;
    text(edges(1) + 0.05*(edges(end) - edges(1)), ...
        edges(end) - 0.05*(edges(end) - edges(1)), ...
        ['\rho = ' num2str(round(r, 3))],  'color', 'm', 'fontSize', 12, 'fontWeight', 'bold');
    xlim([fix(logMinP)-1 fix(logMaxP)+1])
    ylim([fix(logMinP)-1 fix(logMaxP)+1])
    xlabel('$\log_{10}(k_{xx})$ [mD]', 'Interpreter','latex','FontSize',12);
    ylabel('$\log_{10}(k_{zz})$ [mD]', 'Interpreter','latex','FontSize',12);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12)
    a = get(gca,'YTickLabel');
    set(gca,'YTickLabel',a,'fontsize',12)
    grid on
    sgtitle(sprintf('N_{sim} = %d, throw window %d',Nsim,k),'FontSize',12);
    set(fh, 'position', [200, 200, 600, 200]);
    exportgraphics(gcf, sprintf('%d_ThrowWindow%d_kf_distribution.pdf',Nsim,k));
end
