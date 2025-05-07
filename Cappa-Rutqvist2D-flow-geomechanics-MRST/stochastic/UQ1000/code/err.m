clear all;
close all;

load("../data/QoI3.mat");
load("../FML/Tmax/06tanh_logTmax_20mem_10rec_random_100budget_1/Tmax_pred.mat");
load("../FML/Pmax/06tanh_logPmax_20mem_10rec_random_100budget_1/Pmax_pred.mat");
load("../FML/W/06tanh_logW_20mem_10rec_random_100budget_1/W_pred.mat");

for m = 1:1000
    Tmax(L06(m)+1:end,:,m) = nan;
    T_pred(L06(m)+1:end,:,m) = nan;
    dPmax(L06(m)+1:end,:,m) = nan;
    P_pred(L06(m)+1:end,:,m) = nan;
    dMw(L06(m)+1:end,:,m) = nan;
    W_pred(L06(m)+1:end,:,m) = nan;
end

e_T = zeros(1000,1);
e_P = zeros(1000,1);
e_W = zeros(1000,1);
for m = 1:1000
    e_T(m) = nansum((T_pred(:,1,m)-Tmax(:,1,m)).^2)/nansum((Tmax(:,1,m)).^2);
    e_P(m) = nansum((P_pred(:,1,m)-dPmax(:,1,m)).^2)/nansum((dPmax(:,1,m)).^2);
    e_W(m) = nansum((W_pred(:,1,m)-dMw(:,1,m)).^2)/nansum((dMw(:,1,m)).^2);
end

mean(e_T)

mean(e_P)

mean(e_W)





