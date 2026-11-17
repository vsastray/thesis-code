function [target_exp, esti_exp] = sampling(state_in, N_sample, JN, J1, J2, overhead, H)
%% sampling process
noisy_state = ApplyMap(state_in, JN);
target_exp = real(trace(H * state_in));

[vec, val] = eig(H);
out_J1 = ApplyMap(noisy_state, J1);
out_J2 = ApplyMap(noisy_state, J2);

p_J1 = real(trace(out_J1));
p_J2 = real(trace(out_J2));
rho_J1 = out_J1 / p_J1;
rho_J2 = out_J2 / p_J2;
prob_H1_J1 = real(trace(rho_J1 * vec(:,1) * vec(:,1)'));
prob_H1_J2 = real(trace(rho_J2 * vec(:,1) * vec(:,1)'));

total = 0;
for i = 1:N_sample
    x=rand;
    if x < p_J1
        y = rand;
        if y < prob_H1_J1
            total = total + val(1,1);
        else
            total = total + val(2,2);
        end
    else
        y = rand;
        if y < prob_H1_J2
            total = total - val(1,1);
        else
            total = total - val(2,2);
        end
    end
end

esti_exp = overhead * real(total / N_sample);
