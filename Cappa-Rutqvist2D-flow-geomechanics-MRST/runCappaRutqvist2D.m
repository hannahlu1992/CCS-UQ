function [model, states,initState,reservoir_indx,fault_indx_collect] = runCappaRutqvist2D(varargin)
    mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui

    % setup default option values
    opt = struct('method'             , 'fully coupled' , ...
                 'fluid_model'        , 'blackoil-ccs', ...
                 'inj_rate'           , 0.02, ...
                 'mech_BC'            , 'left bottom roller', ...
                 'reportTimes'        , (12:12:60)*minute,...
                 'nonlinearTolerance' , 1e-6            , ...
                 'splittingTolerance' , 1e-6            , ...
                 'verbose'            , false           , ...
                 'splittingVerbose'   , false           , ...
                 'kf'                 , 1e-16           , ...
                 'lambda'             , 0.7             , ...
                 'k_reservoir'        ,1e-13            , ...
                 'poro_reservoir'     ,0.1              , ...
                 'reservoir_compressibility', 1.45e-4   , ...
                 'initialization',    []);
    opt = merge_options(opt, varargin{:});

    gravity reset y on;
    g = norm(gravity);

    %% 1. Mesh
    depth = 500; % depth of aquifer top surface
    Lx = 2000;Ly = 2000;
    location_f = [500,1500];% [fault center x, fault center z]
    dip = 80;% fault dip (m)
    fW = 2.5;% fault width (m)
    f_D = 125; % fault displacement/throw window length (m)
    h_reservoir = 100;% reservoir thickness (m)
    h_caprock = 150;% caprock thickness (m)
    x_bottom = 500-cotd(dip)*Ly/2-fW/2;
    x_top = 500+cotd(dip)*Ly/2-fW/2;
    left_StretchFact_bottom = x_bottom/(location_f(1)-fW/2);
    left_StretchFact_top = x_top/(location_f(1)-fW/2);
    right_StretchFact_bottom = (Lx-x_bottom-fW)/(Lx-location_f(1)-fW/2);
    right_StretchFact_top = (Lx-x_top-fW)/(Lx-location_f(1)-fW/2);
     
    x1 = 500-fW/2-fW*5:fW:500+fW/2+fW*4;
    y1 = 500:fW:500+Ly;
    G1 = tensorGrid(x1,y1);
    
    x2 = 500-fW/2-fW*5-(fW*2)*5:(fW*2):500-fW/2-fW*5;
    y2 = 500:(fW*2):500+Ly;
    G2 = tensorGrid(x2,y2);
    
    x3 = 500-fW/2-fW*5-(fW*2)*5-(fW*5)*5:(fW*5):500-fW/2-fW*5-(fW*2)*5;
    y3 = 500:(fW*5):500+Ly;
    G3 = tensorGrid(x3,y3);
    
    x4 = linspace(0,500-fW/2-fW*5-(fW*2)*5-(fW*5)*5,17);
    y4 = 500:(fW*10):500+Ly;
    G4 = tensorGrid(x4,y4);
    
    x5 = linspace(500+fW/2+fW*4+(fW*2)*5+(fW*5)*5,Lx,57);
    y5 = 500:(fW*10):500+Ly;
    G5 = tensorGrid(x5,y5);
    
    G = glue2DGrid(G2, G1);
    G = glue2DGrid(G,translateGrid(G2,[fW*20 0]));
    G = glue2DGrid(G3, G);
    G = glue2DGrid(G,translateGrid(G3,[fW*55 0]));
    G = glue2DGrid(G4, G);
    G = glue2DGrid(G, G5);
    
    G = computeGeometry(G);
    
    fault_indx = find(G.cells.centroids(:,1)<location_f(1)+fW/2 & G.cells.centroids(:,1)>location_f(1)-fW/2);
    reservoir_indx = [find(G.cells.centroids(:,1)<location_f(1)-fW/2 &G.cells.centroids(:,2)>location_f(2)-h_reservoir/2 &G.cells.centroids(:,2)<location_f(2)+h_reservoir/2);...
                      find(G.cells.centroids(:,1)>location_f(1)+fW/2 &G.cells.centroids(:,2)>location_f(2)-h_reservoir/2-f_D &G.cells.centroids(:,2)<location_f(2)+h_reservoir/2-f_D)];
    caprock_indx = [find(G.cells.centroids(:,1)<location_f(1)-fW/2 &G.cells.centroids(:,2)>location_f(2)-h_reservoir/2-h_caprock &G.cells.centroids(:,2)<location_f(2)-h_reservoir/2);...
                    find(G.cells.centroids(:,1)<location_f(1)-fW/2 &G.cells.centroids(:,2)>location_f(2)+h_reservoir/2 &G.cells.centroids(:,2)<location_f(2)+h_reservoir/2+h_caprock);...
                    find(G.cells.centroids(:,1)>location_f(1)+fW/2 &G.cells.centroids(:,2)>location_f(2)-h_reservoir/2-h_caprock-f_D &G.cells.centroids(:,2)<location_f(2)-h_reservoir/2-f_D);...
                    find(G.cells.centroids(:,1)>location_f(1)+fW/2 &G.cells.centroids(:,2)>location_f(2)+h_reservoir/2-f_D &G.cells.centroids(:,2)<location_f(2)+h_reservoir/2+h_caprock-f_D)];
    basal_aquifer_indx = [find(G.cells.centroids(:,1)<location_f(1)-fW/2 &G.cells.centroids(:,2)>location_f(2)+h_reservoir/2+h_caprock);...
                      find(G.cells.centroids(:,1)>location_f(1)+fW/2 &G.cells.centroids(:,2)>location_f(2)+h_reservoir/2+h_caprock-f_D)];
            
    
    fault_top_indx = intersect(fault_indx,find(G.cells.centroids(:,2)<location_f(2)-h_caprock-h_reservoir/2-f_D));
    fault_bot_indx = intersect(fault_indx,find(G.cells.centroids(:,2)>location_f(2)-h_caprock-h_reservoir/2-f_D+5*f_D));
    fault_throw1_indx = intersect(fault_indx,find(G.cells.centroids(:,2)>location_f(2)-h_caprock-h_reservoir/2-f_D &...
                        G.cells.centroids(:,2)<location_f(2)-h_caprock-h_reservoir/2-f_D +f_D));
    fault_throw2_indx = intersect(fault_indx,find(G.cells.centroids(:,2)>location_f(2)-h_caprock-h_reservoir/2-f_D +f_D &...
                        G.cells.centroids(:,2)<location_f(2)-h_caprock-h_reservoir/2-f_D +2*f_D));
    fault_throw3_indx = intersect(fault_indx,find(G.cells.centroids(:,2)>location_f(2)-h_caprock-h_reservoir/2-f_D +2*f_D &...
                        G.cells.centroids(:,2)<location_f(2)-h_caprock-h_reservoir/2-f_D +3*f_D));
    fault_throw4_indx = intersect(fault_indx,find(G.cells.centroids(:,2)>location_f(2)-h_caprock-h_reservoir/2-f_D +3*f_D &...
                        G.cells.centroids(:,2)<location_f(2)-h_caprock-h_reservoir/2-f_D +4*f_D));
    fault_throw5_indx = intersect(fault_indx,find(G.cells.centroids(:,2)>location_f(2)-h_caprock-h_reservoir/2-f_D +4*f_D &...
                        G.cells.centroids(:,2)<location_f(2)-h_caprock-h_reservoir/2-f_D +5*f_D));
    assert(length(fault_throw1_indx)+length(fault_throw2_indx)+length(fault_throw3_indx)+length(fault_throw4_indx)...
        +length(fault_throw5_indx)+length(fault_top_indx)+length(fault_bot_indx)==length(fault_indx), 'sum of throw windows must equal fault length');
    fault_indx_collect = struct('top', fault_top_indx, 'bot', fault_bot_indx, 'throw1', fault_throw1_indx, 'throw2', fault_throw2_indx,...
                        'throw3', fault_throw3_indx,'throw4', fault_throw4_indx,'throw5', fault_throw5_indx);
           
    makeSkewL = @(c) ((left_StretchFact_top-left_StretchFact_bottom)/Ly*(Ly+depth-c(:,2))+left_StretchFact_bottom).*c(:,1);
    left_indx = find(G.nodes.coords(:,1)<= location_f(1)-fW/2);
    
    makeSkewR = @(c) Lx-((right_StretchFact_top-right_StretchFact_bottom)/Ly*(Ly+depth-c(:,2))+right_StretchFact_bottom).*(Lx-c(:,1));
    right_indx = find(G.nodes.coords(:,1)>= location_f(1)+fW/2);
    
        
    G.nodes.coords(left_indx,1) = makeSkewL(G.nodes.coords(left_indx,:));
    G.nodes.coords(right_indx,1) = makeSkewR(G.nodes.coords(right_indx,:));
    G = computeGeometry(G);

%% Grid visualization    
%     figure
%     plotGrid(G,'FaceColor','none');
%     plotGrid(G,fault_indx,'FaceColor','g');
%     plotGrid(G,reservoir_indx,'FaceColor','b');
%     plotGrid(G,caprock_indx,'FaceColor','y');
% %     plotGrid(G,well_indx,'FaceColor','r');
%     set(gca, 'YDir','reverse');

    %% 2. Rock Property
    % Here we define the porosity and permeability [L^2] of each grid cell. 
    poro = 0.1;
    rock.poro = repmat(poro, G.cells.num, 1);
    rock.poro(caprock_indx) = 0.01;
    rock.poro(reservoir_indx) = opt.poro_reservoir; 
    
    perm = 1e-14;
    rock.perm = repmat(perm, G.cells.num, 2);
    rock.perm(reservoir_indx,:) = repmat(opt.k_reservoir*ones(1,2),size(reservoir_indx));
    rock.perm(caprock_indx,:) = repmat(1e-20*ones(1,2),size(caprock_indx));
    rock.perm(basal_aquifer_indx,:) = repmat(1e-15*ones(1,2),size(basal_aquifer_indx));
    
    if max(size(opt.kf)) == 1
        rock.perm(fault_indx,:) = repmat([opt.kf opt.kf],size(fault_indx)); 
    else
        opt.kf = reshape(opt.kf,[2,5]);opt.kf = opt.kf';
        rock.perm(fault_top_indx,:) = repmat(1e-14*ones(1,2),size(fault_top_indx));
        rock.perm(fault_bot_indx,:) = repmat(1e-16*ones(1,2),size(fault_bot_indx));
        rock.perm(fault_throw1_indx,:) = repmat(opt.kf(1,:),size(fault_throw1_indx));
        rock.perm(fault_throw2_indx,:) = repmat(opt.kf(2,:),size(fault_throw2_indx));
        rock.perm(fault_throw3_indx,:) = repmat(opt.kf(3,:),size(fault_throw3_indx));
        rock.perm(fault_throw4_indx,:) = repmat(opt.kf(4,:),size(fault_throw4_indx));
        rock.perm(fault_throw5_indx,:) = repmat(opt.kf(5,:),size(fault_throw5_indx));
    end

    %% Visualize porosity and perm
    % figure
    % plotCellData(G,rock.poro,'EdgeColor','None');
    % colorbar;
    % set(gca, 'YDir','reverse');
    % ylabel('z(m)')
    % xlabel('x(m)')
    % 
    % figure
    % plotCellData(G,log10(rock.perm(:,2)),'EdgeColor','k');
    % colorbar;
    % set(gca, 'YDir','reverse');


    %% 3. Setup fluid parameters
    switch opt.fluid_model
      case 'blackoil-ccs'
        % Define 3 Saturation regions
        rock.regions.saturation = ones(G.cells.num, 1);        % upper aquifer
        rock.regions.saturation(caprock_indx) = 2;              % caprock
        rock.regions.saturation(fault_indx) = 3;                % fault
        rock.regions.rocknum = ones(G.cells.num,1); 
        fn = 'fPcFault_co2brine_modif_3ph.DATA';
        deck = convertDeckUnits(readEclipseDeck(fn));
        deck.PROPS.ROCK(2) = opt.reservoir_compressibility*1e-6;
        deck.REGIONS.ROCKNUM = rock.regions.rocknum;
        fluid = initDeckADIFluid(deck);
      case 'blackoil-wg'
        pth = getDatasetPath('spe1');
        fn  = fullfile(pth, 'BENCH_SPE1.DATA');
        deck = readEclipseDeck(fn);
        deck = convertDeckUnits(deck);
        fluid = initDeckADIFluid(deck);
        if isfield(fluid, 'pcOW')
            fluid = rmfield(fluid, 'pcOW');
        end
        if isfield(fluid, 'pcOG')
            fluid = rmfield(fluid, 'pcOG');
        end
        % Setup quadratic relative permeabilities, since SPE1 relperm are a bit rough.
        fluid.krW = @(s) s.^2;
        fluid.krG = @(s) s.^2;
        fluid.krOW = @(s) s.^2;
        fluid.krOG = @(s) s.^2;
      case 'oil water'
        fluid = initSimpleADIFluid('phases', 'WO', 'mu', [1, 10]*centi*poise, ...
                                   'n',  [1, 1], 'rho', [1000, 700]*kilogram/ ...
                                   meter^3, 'c', 1e-10*[1, 1], 'cR', 4e-10);

      case 'water'
        fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                                   1000, 'c', 0*1e-10, 'cR',4e-10);
      otherwise
        error('fluid_model  not recognized.');
    end

    %% 4. Setup material parameters for Biot and mechanics
    E          = 10 * giga * Pascal; % Young's module
    nu         = 0.25;               % Poisson's ratio
    alpha      = 1;                 % Biot's coefficient
    rock.alpha = alpha * ones(G.cells.num, 1);
    bulk_density = 2260*kilogram/meter^3;

    %% 5.1. ICs and BCs for flow
    f_top = find(G.faces.centroids(:,2) == min(G.faces.centroids(:,2)));
    f_bot = find(G.faces.centroids(:,2) == max(G.faces.centroids(:,2)));
    f_right = find(G.faces.centroids(:,1) == max(G.faces.centroids(:,1)));
    
    f = [f_top;f_bot;f_right];
    [y_top, y_bot] = deal(min(G.faces.centroids(f,2)), max(G.faces.centroids(f,2)));

    clear initState;
    switch opt.fluid_model
      case 'blackoil-ccs'
        rho_or = fluid.rhoOS*kilogram/meter^3;
        p_top = 1*barsa + g*rho_or*y_top; % p at y_top
        equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS(1), [y_top, y_bot], p_top);
        fp_val = reshape(deval(equil, G.faces.centroids(f,2)), [], 1);
        bc = addBC([], f, 'pressure', fp_val, 'sat', [0, 1, 0]);
        if ~isempty(opt.initialization)
            initState = opt.initialization;
        else
            initState.pressure = reshape(deval(equil, G.cells.centroids(:,2)), [], 1);
            initState.s  = repmat([0, 1, 0], [G.cells.num, 1]);
            initState.rs  = 0*fluid.rsSat(initState.pressure);
        end
        clear equil;
      case 'blackoil-wg'
        rho_wr = fluid.rhoWS*kilogram/meter^3;
        p_top = 1*barsa + g*rho_wr*y_top; % p at y_top
        equil  = ode23(@(z,p) g *fluid.rhoWS(1), [y_top, y_bot], p_top);
        fp_val = reshape(deval(equil, G.faces.centroids(f,2)), [], 1);
        
        bc = addBC([], f, 'pressure', fp_val, 'sat', [1, 0, 0]);
        initState.pressure = reshape(deval(equil, G.cells.centroids(:,2)), [], 1);
        initState.s  = repmat([1, 0, 0], [G.cells.num, 1]);
        initState.rs  = fluid.rsSat(initState.pressure);
        clear equil;
      case 'oil water'
        rho_or = fluid.rhoOS*kilogram/meter^3;
        p_top = 1*barsa + g*rho_or*y_top; % p at y_top
        equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS(1), [y_top, y_bot], p_top);
        fp_val = reshape(deval(equil, G.faces.centroids(f,2)), [], 1);
        bc = addBC([], f, 'pressure', fp_val, 'sat', [0, 1]);
        initState.pressure = reshape(deval(equil, G.cells.centroids(:,2)), [], 1);
        initState.s  = repmat([0, 1], [G.cells.num, 1]); 
        clear equil;
      case 'water'
        rho_wr = fluid.rhoWS*kilogram/meter^3;
        p_top = 1*barsa + g*rho_wr*y_top; % p at y_top
        equil  = ode23(@(z,p) g .* fluid.bW(p,0,false)*fluid.rhoWS(1), [y_top, y_bot], p_top);
        fp_val = reshape(deval(equil, G.faces.centroids(f,2)), [], 1);
        bc = addBC([], f, 'pressure', fp_val, 'sat', 1);
        initState.pressure = reshape(deval(equil, G.cells.centroids(:,2)), [], 1);
        initState.s  = ones(G.cells.num, 1); 
        clear equil;
      otherwise
        error('fluid_model not recognized.')
    end
    %% Visualize initial pressure
    % figure 
    % plotCellData(G, initState.pressure/barsa, 'edgealpha', 0.2);
    % colormap(jet), c = colorbar;
    % c.Label.Interpreter = 'latex'; c.Label.FontSize = 11;
    % c.Label.String = '$p_0 $ [bar]';
    % set(gca, 'YDir','reverse');


    %% 5.2 BCs for mechanics
    switch opt.mech_BC
      case 'left bottom roller'
        bottom_nodes = find(G.nodes.coords(:,2) ==max(G.nodes.coords(:,2)));
        left_nodes = find(G.nodes.coords(:,1) ==min(G.nodes.coords(:,1)));
        
        left_bottom_nodes = intersect(left_nodes,bottom_nodes);
        left_nodes = setdiff(left_nodes,left_bottom_nodes);
        bottom_nodes = setdiff(bottom_nodes,left_bottom_nodes);
        
        disp_bc.nodes = [bottom_nodes;left_nodes;left_bottom_nodes];
        disp_bc.uu = repmat([0,0],G.nodes.num,1);
        disp_bc.mask = [repmat([false,true],numel(bottom_nodes),1);...% roller btm
                        repmat([true,false],numel(left_nodes),1);... % roller left
                        repmat([true,true],numel(left_bottom_nodes),1)];% left bottom corner fixed
        
        %% impose a given pressure at top and right, right sigma_h = 0.7*sigma_v
        pb_top = 1*barsa + g*bulk_density*y_top;
        equil  = ode23(@(z,p) g*bulk_density, [y_top, y_bot], pb_top); 
        pb_right_val = reshape(deval(equil, G.faces.centroids(f_right,2)), [], 1);
        clear equil;
        force = [repmat([0, pb_top],numel(f_top),1);[-opt.lambda*pb_right_val,zeros(length(f_right),1)]];
        force_bc = struct('faces', [f_top;f_right], 'force', force);
    %% other cases not implemented yet    
    otherwise
        error('mech_BC not recognized.')
    end
    el_bc = struct('disp_bc' , disp_bc,'force_bc', force_bc);   

    %% Setup load for mechanics
    loadfun = @(x) repmat(bulk_density* [0 g], size(x, 1), 1);

    %% 6. Gather all the mechanical parameters in a struct
    % Transform these global properties (uniform) to cell values.
    E_matrix          = repmat(E, G.cells.num, 1);
    E_matrix(fault_indx) = E_matrix(fault_indx)/2;   %5 * giga * Pascal Young's module 
    C = Enu2C(E_matrix, nu*ones(G.cells.num,1),G);
    mech = struct('E', E_matrix, 'nu', repmat(nu, G.cells.num, 1), 'el_bc', el_bc, 'load', loadfun);

    %% 7. Setup model
    modeltype = [opt.method, ' and ', opt.fluid_model];
    fullycoupledOptions = {'verbose', opt.verbose};
    splittingOptions = {'splittingTolerance', opt.splittingTolerance, ...
                        'splittingVerbose', opt.splittingVerbose};
    switch modeltype
      case {'fully coupled and water'}
        model = MechWaterModel(G, rock, fluid, mech, fullycoupledOptions{:});
      case 'fully coupled and oil water'
        model = MechOilWaterModel(G, rock, fluid, mech, fullycoupledOptions{:});
      case {'fully coupled and blackoil-ccs','fully coupled and blackoil-wg'}
        model = MechBlackOilModel(G, rock, fluid, mech, fullycoupledOptions{:});
      case {'fixed stress splitting and blackoil-ccs','fixed stress splitting and blackoil-wg'}
        model = MechFluidFixedStressSplitModel(G, rock, fluid, mech, ...
                                               'fluidModelType', 'blackoil', ...
                                               splittingOptions{:});
      case 'fixed stress splitting and oil water'
        model = MechFluidFixedStressSplitModel(G, rock, fluid, mech, ...
                                               'fluidModelType', 'oil water', ...
                                               splittingOptions{:});
      case {'fixed stress splitting and water'}
        model = MechFluidFixedStressSplitModel(G, rock, fluid, mech, ...
                                               'fluidModelType', 'water', ...
                                               splittingOptions{:});
      otherwise
        error('modeltype not recognized.');
    end


    %% 8. Setup wells, well cell indices in 'global' grid
    wc = find(G.cells.centroids(:,1)< 25 & G.cells.centroids(:,2)>1500-10 & G.cells.centroids(:,2)<1500+20);
    switch opt.fluid_model
      case 'blackoil-wg'
            W = addWell([ ], G, rock, wc, 'Name', 'I1', 'Dir', 'z', ...
                        'Type', 'rate', 'Val', opt.inj_rate/fluid.rhoGS, 'compi', [0, 0, 1], ...    
                        'refDepth', G.cells.centroids(wc, G.griddim), ... 
                         'Radius', 0.2);
      case 'blackoil-ccs'
            W = addWell([ ], G, rock, wc, 'Name', 'I1', 'Dir', 'z', ...
                        'Type', 'rate', 'Val', opt.inj_rate/fluid.rhoGS, 'compi', [0, 0, 1], ...    
                        'refDepth', G.cells.centroids(wc, G.griddim), ... 
                         'Radius', 0.2);
      case 'oil water'
          W = addWell([ ], G, rock, wc, 'Name', 'I1', 'Dir', 'z', ...
                'Type', 'rate', 'Val', opt.inj_rate/fluid.rhoWS, 'compi', [1, 0], ...    
                'refDepth', G.cells.centroids(wc, G.griddim), ... 
                 'Radius', 0.2);
      case 'water'
            W = addWell([ ], G, rock, wc, 'Name', 'I1', 'Dir', 'z', ...
                'Type', 'rate', 'Val', opt.inj_rate/fluid.rhoWS, 'compi', 1, ...    
                'refDepth', G.cells.centroids(wc, G.griddim), ... 
                 'Radius', 0.2);
      otherwise
        error('fluid_model not recognized.')
    end

    facilityModel = FacilityModel(model.fluidModel);
    facilityModel = facilityModel.setupWells(W);
    model.FacilityModel = facilityModel;

    %% 9. Setup schedule
    timesteps = [opt.reportTimes(1) diff(opt.reportTimes)];
    assert(sum(timesteps)==opt.reportTimes(end), 'sum of timesteps must equal simTime');
    schedule.control    = struct('W', W,'bc',bc);
    schedule.step.val = timesteps;
    schedule.step.control = [ones(numel(timesteps), 1)];
    % initial nodal displacements and pore pressure values are all zero
    % only need to determine the number of displacement values that must be
    % provided.
    if ~isempty(opt.initialization)
        initState = opt.initialization;
    else
        [uu,~] = VEM_linElast(createAugmentedGrid(G),C,el_bc,loadfun,'pressure',initState.pressure);
        initState.uu = uu;
        initState.u = reshape(initState.uu',[G.nodes.num*2,1]);
        initState.xd = initState.u(~model.mechModel.operators.isdirdofs);
        initState = addDerivedQuantities(model.mechModel, initState);
    end
   
    %% 10. Simulate
    solver = NonLinearSolver('maxIterations', 100);
    [~, states, ~] = simulateScheduleAD(initState, model, schedule, 'nonlinearsolver', ...
                                        solver);
    fm = model.fluidModel; % to avoid excess typing
    new_model = GenericBlackOilModel(fm.G,fm.rock,fm.fluid,'gas',fm.gas,...
        'oil',fm.oil,'water',fm.water,'disgas',fm.disgas,'vapoil',fm.vapoil);
    new_model = new_model.validateModel();
    model = model.validateModel();
    initState.ComponentTotalMass = new_model.getProp(initState,'ComponentTotalMass');
    %% at this pot, new model should be endowed with separate oil, gas and water Components.
    %% Compute total mass:
    for n=1:numel(states)
        states{n}.ComponentTotalMass = new_model.getProp(states{n},'ComponentTotalMass');
%         states{n}.Dp = states{n}.pressure - states{n-1}.pressure;
%         states{n}.Dxd = states{n}.xd -states{n-1}.xd;
%         states{n}.Du = states{n}.u-states{n-1}.u; % 
%         states{n}.Duu = states{n}.uu- states{n-1}.uu; % full lists of
%                                                    % displacement values, including the fixed displacement values associated
%                                                    % with imposed boundary conditions
%         states{n}.Dvdiv = states{n}.vdiv- states{n-1}.vdiv; % cell-wise volumetric strain
%         states{n}.Dstress = states{n}.stress- states{n-1}.stress; % stress tensor
%         states{n}.Dstrain = states{n}.strain- states{n-1}.strain; % strain tensor
%     pc = model.getProp(states{n}, 'CapillaryPressure');
%     states{n}.FlowProps.CapillaryPressure = pc{1,2};
%     states{n}.FlowProps.RelativePermeability = model.getProp(states{n}, ...
%                                                              'RelativePermeability');
    end
end