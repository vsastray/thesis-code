%% Required packages:
% CVX http://cvxr.com/cvx/download/
% QETLAB http://www.qetlab.com/Main_Page
% http://www.dr-qubit.org/download.php?file=matlab/quantinf.tar.gz
clear;
dim_target = 2;
dim_noise_in = 2;
dim_noise_out = 2;

dim_noise = dim_noise_in * dim_noise_out;

MES = MaxEntangled(dim_target, 1) * MaxEntangled(dim_target, 1)';
JI = MES * dim_target;  %choi matrix of the target noiseless channel

cost_budget = 1.05;

% noise_channel = DephasingChannel(dim_noise_in, 0.9);
noise_channel = DepolarizingChannel(dim_noise_in, 0.9);
% amplitude damping
eps = 0.1;
A{1}=[1 0; 0 sqrt(1-eps)];
A{2}=[0 sqrt(eps); 0 0];
AmplitudeDamping = chanconv(A, 'kraus', 'choi', [2 2]);  %choi matrix of the AD channel
% noise_channel = AmplitudeDamping;


%% simplified SDP for the minimum error of shadow communication
noise_sys = [dim_noise_in, dim_noise_out];
cvx_begin sdp quiet
    variable T_pos(dim_noise, dim_noise) hermitian
    variable T_neg(dim_noise, dim_noise) hermitian
    variable V_pos(dim_noise_in, dim_noise_in) hermitian
    variable V_neg(dim_noise_in, dim_noise_in) hermitian
    variable Z(dim_target^2, dim_target^2) hermitian
    variable epsilon
    JM = trace(noise_channel.' * (T_pos - T_neg)) * (dim_target^2 * MES - eye(dim_target^2)) / (dim_target * (dim_target^2 - 1)) + trace(V_pos - V_neg) * (eye(dim_target^2) - MES) / (dim_target * (dim_target^2 - 1));
    eps_simplified = epsilon;
    minimize eps_simplified
    subject to
        Z >= 0; Z >= JI - JM; PartialTrace(Z, 2) <= (2 * epsilon + 1 - trace(V_pos - V_neg) / (dim_target^2)) * eye(dim_target) / 2;
        trace(V_pos + V_neg) / (dim_target^2) <= cost_budget;
        T_pos >= 0; T_neg >= 0; T_pos <= kron(V_pos, eye(dim_noise_out)); T_neg <= kron(V_neg, eye(dim_noise_out));
        PartialTrace(T_pos, 1, noise_sys) == trace(V_pos) * eye(dim_noise_out) / (dim_target^2);
        PartialTrace(T_neg, 1, noise_sys) == trace(V_neg) * eye(dim_noise_out) / (dim_target^2);
cvx_end
cvx_status_simplified = cvx_status


%% SDP for the minimum error of quantum communication
superop_sys = [dim_noise_in, dim_noise_out, dim_target, dim_target];
cvx_begin sdp quiet
    variable JPI(dim_target^2 * dim_noise, dim_target^2 * dim_noise) hermitian
    variable Z(dim_target^2, dim_target^2) hermitian
    variable epsilon
    JE = PartialTrace(JPI, [2, 4], superop_sys) / dim_noise_out;
    JD = PartialTrace(JPI, [1, 3], superop_sys) / dim_target;
    JM = PartialTrace(kron(noise_channel.', eye(dim_target^2)) * JPI, 1, [dim_noise, dim_target^2]);  % effective final map M = PI circ N
    eps_quantum = epsilon;
    minimize eps_quantum
    subject to
        Z >= 0; Z >= JI - JM; PartialTrace(Z, 2) <= epsilon * eye(dim_target);
        JPI >= 0; PartialTrace(JPI, [1, 4], superop_sys) == eye(dim_target * dim_noise_out); % CPTP
        PartialTrace(JPI, 4, superop_sys) == PermuteSystems(kron(JE, eye(dim_noise_out)), [1 3 2], [dim_noise_in, dim_target, dim_noise_out]); % non-signalling from Bob to Alice
        PartialTrace(JPI, 1, superop_sys) == PermuteSystems(kron(JD, eye(dim_target)), [1 3 2], [dim_noise_out, dim_target, dim_target]); % non-signalling from Alice to Bob
cvx_end
cvx_status_quantum = cvx_status


%% output
eps_simplified
eps_quantum
