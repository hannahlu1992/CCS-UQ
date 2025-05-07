function [Tmax,dPmax,dMw] = computing3QoIs(states,initState,reservoir_indx,fault_indx_collect,dip)
fault_indx = [fault_indx_collect.top;fault_indx_collect.throw1;fault_indx_collect.throw2;...
fault_indx_collect.throw3;fault_indx_collect.throw4;fault_indx_collect.throw5;fault_indx_collect.bot];

delta = 90-dip;
for n = 1:numel(states)
    %% states{n}.Stress is stress, states{n}.stress is effective stress
    states{n}.Stress = states{n}.stress -...
        [states{n}.pressure, states{n}.pressure, 0*states{n}.pressure];
    states{n}.pstress_max = zeros(size(states{n}.stress,1),1);% sigma1   %> 0, extensive; <0 compressive
    states{n}.pstress_min = zeros(size(states{n}.stress,1),1);% sigma3
    states{n}.ShearStress = zeros(size(states{n}.stress,1),1);% tau
    states{n}.NormalStress = zeros(size(states{n}.stress,1),1);% sigma_n
    states{n}.EffectiveNormalStress = zeros(size(states{n}.stress,1),1);% sigma_n'
    states{n}.StabilityIndicator = zeros(size(states{n}.stress,1),1);% sigma_n'
    for i = 1:numel(fault_indx)
        [~,D] = eigs([states{n}.Stress(fault_indx(i),1),states{n}.Stress(fault_indx(i),3);...
            states{n}.Stress(fault_indx(i),3),states{n}.Stress(fault_indx(i),2)]);
        states{n}.pstress_max(fault_indx(i)) = D(1,1);
        states{n}.pstress_min(fault_indx(i)) = D(2,2);
        states{n}.ShearStress(fault_indx(i)) = (states{n}.pstress_max(fault_indx(i))-states{n}.pstress_min(fault_indx(i)))/2*sind(2*delta);
        states{n}.NormalStress(fault_indx(i)) = (states{n}.pstress_max(fault_indx(i))+states{n}.pstress_min(fault_indx(i)))/2-...
            (states{n}.pstress_max(fault_indx(i))-states{n}.pstress_min(fault_indx(i)))/2*cosd(2*delta);
        states{n}.EffectiveNormalStress(fault_indx(i)) = states{n}.NormalStress(fault_indx(i)) + states{n}.pressure(fault_indx(i));
        states{n}.StabilityIndicator(fault_indx(i)) = states{n}.ShearStress(fault_indx(i))/states{n}.EffectiveNormalStress(fault_indx(i));
    end
end
fault_sigma_n_prime = zeros(numel(states),numel(fault_indx));
fault_tau = fault_sigma_n_prime;
fault_p = fault_sigma_n_prime;
fault_pL = fault_sigma_n_prime;
for n = 1:numel(states)
    fault_sigma_n_prime(n,:) = states{n}.EffectiveNormalStress(fault_indx)';
    fault_tau(n,:) = states{n}.ShearStress(fault_indx)';
    fault_p(n,:) = states{n}.pressure(fault_indx)';
    fault_pL(n,:) = states{n}.pressure(fault_indx-1)';
end
dP = fault_p-initState.pressure(fault_indx)';
T = fault_tau./fault_sigma_n_prime;
Tmax = T(:,420);
dPmax = dP(:,420);

Mw0 = sum(initState.ComponentTotalMass{2}(reservoir_indx));
dMw = zeros(numel(states),1);
for m = 1:numel(states)
    dMw(m,1) = Mw0 - sum(states{m}.ComponentTotalMass{2}(reservoir_indx));
end

end