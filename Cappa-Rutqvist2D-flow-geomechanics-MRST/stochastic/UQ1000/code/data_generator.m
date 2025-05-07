clear all;
close all;

addpath("../../../../mrst-developer");
addpath("../../../");
run startup.m

warning('off','all')
mrstModule add diagnostics ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui 

%% write intro text for each case
writeIntroText = @(opt)(fprintf('\n*** Start new simulation\n* fluid model : %s\n* method : %s\n* injection rate : %.2f\n* fault permeability : %.2f\n* bc ratio : %.2f\n* reservoir permeability : %.2f\n* reservoir rock compressibility : %s\n* reservoir porosity : %.2f\n\n', opt.fluid_model, opt.method, opt.inj_rate,log10(opt.kf(1)), opt.lambda,log10(opt.k_reservoir),opt.reservoir_compressibility,opt.poro_reservoir));

%% generate parameters
% rng(1);
% params = rand(4,1000);
load('../data/params.mat');
lambda = 0.6+0.2*params(1,:);
k_reservoir = 10.^(log10(5e-14)+(log10(5e-13)-log10(5e-14))*params(2,:));
poro_reservoir = 0.1+0.2*params(3,:);
cr_reservoir =  10.^(log10(3e-5)+(log10(3e-4)-log10(3e-5))*params(4,:));


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
opt.reportTimes = [(1:4:800)*hour];

load("../data/1000kf.mat",'kf_ups_5');
kf_ups_5 = permute(kf_ups_5,[3,1,2]);
M = size(kf_ups_5,1);
Tmax = zeros(200,1,M);
dPmax = zeros(200,1,M);
dMw = zeros(200,1,M);

for m = 1:M
    opt.lambda = lambda(m);
    opt.k_reservoir = k_reservoir(m);
    opt.poro_reservoir = poro_reservoir(m);
    opt.reservoir_compressibility = cr_reservoir(m);
    opt.initialization = [];
    opt.kf = kf_ups_5(m,[1,3],:);
    writeIntroText(opt);
    optvals = cellfun(@(x) opt.(x), fieldnames(opt), 'uniformoutput', false);
    optlist = reshape(vertcat(fieldnames(opt)', optvals'), [], 1);
    [model, states,initState,reservoir_indx,fault_indx_collect]= runCappaRutqvist2D(optlist{:});
    filename = sprintf('../data/states_%d.mat',m);
    save(filename,'states','initState');
    S = cell(numel(states)+1,1);
    S{1} = initState;
    S(2:end) = states;
    %% Compute QoIs
    [Tmax(:,:,m),dPmax(:,:,m),dMw(:,:,m)] = computing3QoIs(states,initState,reservoir_indx,fault_indx_collect,80);
end

%% Postprocessing
L = zeros(M,1);
for m = 1:M
    if max(Tmax(:,1,m))<=1
        L(m) = 200;
    else
        L(m) = find(Tmax(:,1,m)>1,1)-1;
        Tmax(L(m)+1:end,1,m) = NaN;
        dPmax(L(m)+1:end,1,m) = NaN;
        dMw(L(m)+1:end,1,m) = NaN;
    end
end

L06 = zeros(M,1);
for m = 1:M
    if max(Tmax(:,1,m))<=0.6
        L06(m) = 200;
    else
        L06(m) = find(Tmax(:,1,m)>0.6,1)-1;
    end
end
save('../data/QoI3.mat','L','L06','Tmax','dPmax','dMw');