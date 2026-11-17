%% Required packages:
% CVX http://cvxr.com/cvx/download/
% QETLAB http://www.qetlab.com/Main_Page
% http://www.dr-qubit.org/download.php?file=matlab/quantinf.tar.gz
clear;
dim_id = 2;
dim_target_in = 4;
dim_target_out = 4;

dim_target = dim_target_in * dim_target_out;
noise_sys = [dim_target_in, dim_target_out];

cost_budget = 1.05;

% noise_channel = PermuteSystems(kron(DephasingChannel(2, 0.9), DephasingChannel(2, 0.9)), [1 3 2 4]);
noise_channel = PermuteSystems(kron(DepolarizingChannel(2, 0.9), DepolarizingChannel(2, 0.9)), [1 3 2 4]);
% amplitude damping
eps = 0.1;
A{1}=[1 0; 0 sqrt(1-eps)];
A{2}=[0 sqrt(eps); 0 0];
AmplitudeDamping = chanconv(A, 'kraus', 'choi', [2 2]);  %choi matrix of the AD channel
% noise_channel = PermuteSystems(kron(AmplitudeDamping, AmplitudeDamping), [1 3 2 4]);


%% simplified SDP for the minimum error of shadow simulation
cvx_begin sdp quiet
    variable Y_pos(dim_target, dim_target) hermitian
    variable Y_neg(dim_target, dim_target) hermitian
    variable V_pos(dim_target_out, dim_target_out) hermitian
    variable V_neg(dim_target_out, dim_target_out) hermitian
    variable Z(dim_target, dim_target) hermitian
    variable epsilon
    eps_simplified = epsilon;
    minimize eps_simplified
    subject to
        Z >= 0; Z >= noise_channel - Y_pos + Y_neg; PartialTrace(Z, 2, noise_sys) <= (2 * epsilon + 1 - trace(V_pos - V_neg) / (dim_id^2)) * eye(dim_target_in) / 2;
        trace(V_pos + V_neg) / (dim_id^2) <= cost_budget;
        Y_pos >= 0; Y_neg >= 0; Y_pos <= kron(eye(dim_target_in), V_pos); Y_neg <= kron(eye(dim_target_in), V_neg);
        PartialTrace(Y_pos, 2, noise_sys) == trace(V_pos) * eye(dim_target_in) / (dim_id^2);
        PartialTrace(Y_neg, 2, noise_sys) == trace(V_neg) * eye(dim_target_in) / (dim_id^2);
cvx_end
cvx_status_simplified = cvx_status


%% SDP for the minimum error of quantum simulation
cvx_begin sdp quiet
    variable JM(dim_target, dim_target) hermitian
    variable V(dim_target_out, dim_target_out) hermitian
    variable Z(dim_target, dim_target) hermitian
    variable epsilon
    eps_quantum = epsilon;
    minimize eps_quantum
    subject to
        Z >= 0; Z >= JM - noise_channel; PartialTrace(Z, 2, noise_sys) <= epsilon * eye(dim_target_in);
        JM >= 0; PartialTrace(JM, 2, noise_sys) == eye(dim_target_in); % CPTP
        JM <= kron(eye(dim_target_in), V); trace(V) == dim_id^2;
cvx_end
cvx_status_quantum = cvx_status


%% output
eps_simplified
eps_quantum
