%% find the most isolated parameter points
clear all;
close all;
load("../data/1000kf.mat");
load("../data/params.mat");
load("../data/QoI3.mat");

kf = kf_ups_5([1,3],:,:);
kf = reshape(log10(kf),[10,1000]);
params = [(kf-min(min(kf)))./(max(max(kf))-min(min(kf)))-0.5;params-0.5];

[Np,N_mc] = size(params);

N_budget = 100;
idx_sparse = sparse_filling(params,N_budget);
idx_miniMax = miniMax(params,N_budget);
idx_Maxmin = Maxmin(params,N_budget);
s = RandStream('mlfg6331_64'); 
idx_random = randsample(s,N_mc,N_budget);

save('../data/training_idx.mat','idx_random',"idx_Maxmin","idx_miniMax","idx_sparse");

figure
subplot(2,2,1)
for m = 1:1000
hold on;
plot(Tmax(1:L(m),1,m),'k-');
end
for i = 1:N_budget
plot(Tmax(1:L(idx_random(i)),1,idx_random(i)),'r','LineWidth',2);
end
title('random')
subplot(2,2,2)
for m = 1:1000
hold on;
plot(Tmax(1:L(m),1,m),'k-');
end
for i = 1:N_budget
plot(Tmax(1:L(idx_sparse(i)),1,idx_sparse(i)),'r','LineWidth',2);
end
title('sparse filling')
subplot(2,2,3)
for m = 1:1000
hold on;
plot(Tmax(1:L(m),1,m),'k-');
end
for i = 1:N_budget
plot(Tmax(1:L(idx_miniMax(i)),1,idx_miniMax(i)),'r','LineWidth',2);
end
title('miniMax')
subplot(2,2,4)
for m = 1:1000
hold on;
plot(Tmax(1:L(m),1,m),'k-');
end
for i = 1:N_budget
plot(Tmax(1:L(idx_Maxmin(i)),1,idx_Maxmin(i)),'r','LineWidth',2);
end
title('Maxmin')



%% sparse-filling
function idx = sparse_filling(params,N_budget)
    idx = zeros(N_budget,1);
    for i = 1:N_budget
        d = zeros(1,size(params,2));
        for m = 1:size(params,2)
            d(m) = sum((params-repmat(params(:,m),[1,size(params,2)])).^2,"all");
        end
        [~,idx(i)] = max(d);
        params(:,idx(i)) = [];
    end
end

%% miniMax
function idx = miniMax(params,N_budget)
    idx = zeros(N_budget,1);
    d = zeros(1,size(params,2));
    for m = 1:size(params,2)
        d(m) = sum(vecnorm(params-params(:,m)));
    end
    [~,idx(1)] = max(d);
    params_selected = params(:,idx(1));
    params(:,idx(1)) = [];
    
    for i = 2:N_budget
        d = zeros(1,size(params,2));
        for m = 1:size(params,2)
            d(m) = min(vecnorm(params_selected-params(:,m)));
        end
        [~,idx(i)] = max(d);
        params_selected = [params_selected,params(:,idx(i))];
        params(:,idx(i)) = [];
    end
end

%% Maxmin
function idx = Maxmin(params,N_budget)
    idx = zeros(N_budget,1);
    d = zeros(1,size(params,2));
    for m = 1:size(params,2)
        d(m) = sum(vecnorm(params-params(:,m)));
    end
    [~,idx(1)] = max(d);
    params_selected = params(:,idx(1));
    params(:,idx(1)) = [];
    
    for i = 2:N_budget
        d = zeros(size(params_selected,2),size(params,2));
        dmin = zeros(1,size(params,2));
        for m = 1:size(params,2)
            for n = 1:size(params_selected,2)
                d(n,m) = vecnorm(params_selected(n)-params(:,m));
            end
            dmin(m) = min(d(:,m));
        end
        [~,idx(i)] = max(dmin);
        params_selected = [params_selected,params(:,idx(i))];
        params(:,idx(i)) = [];
    end
end




