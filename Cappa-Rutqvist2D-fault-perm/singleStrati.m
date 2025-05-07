function [mySect,faults,smears,G0] = singleStrati(footwall,hangingwall,faultDip,totalThickness,Nsim,thick)
% 1. Optional input parameters
% In this case, we indicate a maximum fault material permeability of and a correlation 
% coefficient for dependent variables:
maxPerm = [];                 % [mD]
rho     = 0.6;                % Corr. coeff. for multivariate distributions

% 2. Flow upscaling options and number of simulations
U.useAcceleration = 1;          % 1 requires MEX setup, 0 otherwise (slower for MPFA).
U.method          = 'mpfa';     % 'tpfa' recommended if useAcceleration = 0
U.outflux         = 0;          % compare outflux of fine and upscaled model
U.ARcheck         = 0;          % check if Perm obtained with grid with aspect ratio of 
                                % only 5 gives same output.
U.coarseDims      = [1 1 1];    % Mandatory one cell if 2D
% Instantiate FaultedSection object (Strati in Faulted Section)
mySect = FaultedSection(footwall, hangingwall, faultDip, 'maxPerm', maxPerm, 'totThick', totalThickness);


% 3. Get material distributions
% We use the inputs to constrain the ranges and distributions for each of the 
% intermediate variables.
% Get material property distributions
mySect = mySect.getMatPropDistr();

% 4. Get base grid
% We generate a base grid with arbitrary thickness, to be modified at each
% realization (faster than generating n grids from scratch)
D = sum(mySect.Tap(mySect.FW.Id));
T0 = 1;
G0 = makeFaultGrid(T0, D);


% 5. Generate intermediate variable samples, calculate smear dimensions and upscale permeability
% We create two container variables (faults and smears) where we'll save all 
% data for each realization. For each realization, the code defines a Fault object, 
% generates intermediate variable samples, calculates the smear dimensions, and, 
% within upscaleSmearPerm, generates a fault material distribution consistent 
% with the inputs and upscales the permeability.
% Generate fault object with properties for each realization
faults = cell(Nsim, 1);
smears = cell(Nsim, 1);
if nargin<6
    thick = [];
end
for n=1:Nsim    % parfor allowed if you have the parallel computing toolbox
    myFault = Fault2D(mySect, faultDip);
    
    % Get material property (intermediate variable) samples, and fix 
    % along-strike thickness of current realization if 3D.
    myFault = myFault.getMaterialProperties(mySect, 'corrCoef', rho);
    
    % Update grid dimensions with sampled fault thickness
    if ~isempty(thick)
        G = updateGrid(G0, thick(n));
        myFault.MatProps.thick = thick(n);
    else
        G = updateGrid(G0, myFault.MatProps.thick);
    end
    
    % Generate smear object with T, Tap, L, Lmax
    smear = Smear(mySect, myFault, G, 1);
    
    % Place fault materials and assign cell-based properties
    myFault = myFault.placeMaterials(mySect, smear, G);
    
    % Compute upscaled permeability distribution
    myFault = myFault.upscaleProps(G, U);
    
    % Save result
    faults{n} = myFault;
    smears{n} = smear;
    if mod(n, 100) == 0
        disp(['Simulation ' num2str(n) ' / ' num2str(Nsim) ' completed.'])
    end
end
end