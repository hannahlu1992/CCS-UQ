function [Dfault_p,slip_tendency,fault_tau,fault_sigma_n_prime] = computingQoIs(states,fault_indx_collect,dip)
fault_indx = [fault_indx_collect.top;fault_indx_collect.throw1;fault_indx_collect.throw2;...
fault_indx_collect.throw3;fault_indx_collect.throw4;fault_indx_collect.throw5;fault_indx_collect.bot];

%% pressure buildup and slip tendency
delta = 90-dip;
outputs = states;
for n = 1:numel(states)
    %% outputs{n}.Stress is stress, outputs{n}.stress is effective stress
    outputs{n}.Stress = outputs{n}.stress -...
        [outputs{n}.pressure, outputs{n}.pressure, 0*outputs{n}.pressure];
    outputs{n}.pstress_max = zeros(size(outputs{n}.stress,1),1);% sigma1   %> 0, extensive; <0 compressive
    outputs{n}.pstress_min = zeros(size(outputs{n}.stress,1),1);% sigma3
    outputs{n}.ShearStress = zeros(size(outputs{n}.stress,1),1);% tau
    outputs{n}.NormalStress = zeros(size(outputs{n}.stress,1),1);% sigma_n
    outputs{n}.EffectiveNormalStress = zeros(size(outputs{n}.stress,1),1);% sigma_n'
    outputs{n}.StabilityIndicator = zeros(size(outputs{n}.stress,1),1);% T
    for i = 1:size(outputs{n}.stress,1)
        [~,D] = eigs([outputs{n}.Stress(i,1),outputs{n}.Stress(i,3);...
            outputs{n}.Stress(i,3),outputs{n}.Stress(i,2)]);
        outputs{n}.pstress_max(i) = D(1,1);
        outputs{n}.pstress_min(i) = D(2,2);
        outputs{n}.ShearStress(i) = (outputs{n}.pstress_max(i)-outputs{n}.pstress_min(i))/2*sind(2*delta);
        outputs{n}.NormalStress(i) = (outputs{n}.pstress_max(i)+outputs{n}.pstress_min(i))/2-...
            (outputs{n}.pstress_max(i)-outputs{n}.pstress_min(i))/2*cosd(2*delta);
        outputs{n}.EffectiveNormalStress(i) = outputs{n}.NormalStress(i) + outputs{n}.pressure(i);
        outputs{n}.StabilityIndicator(i) = outputs{n}.ShearStress(i)/outputs{n}.EffectiveNormalStress(i);
    end
end
fault_sigma_n_prime = zeros(numel(outputs),numel(fault_indx));
fault_tau = fault_sigma_n_prime;
fault_p = fault_sigma_n_prime;
for n = 1:numel(outputs)
    fault_sigma_n_prime(n,:) = outputs{n}.EffectiveNormalStress(fault_indx)';
    fault_tau(n,:) = outputs{n}.ShearStress(fault_indx)';
    fault_p(n,:) = outputs{n}.pressure(fault_indx)';
end
Dfault_p = fault_p(2:end,:)-fault_p(1,:);
slip_tendency = fault_tau./fault_sigma_n_prime;
end