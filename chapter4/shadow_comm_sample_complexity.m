%% Required packages:
% CVX http://cvxr.com/cvx/download/
% QETLAB http://www.qetlab.com/Main_Page
% http://www.dr-qubit.org/download.php?file=matlab/quantinf.tar.gz
clear;
dim_target = 2;
dim_noise_in = 2;
dim_noise_out = 2;

dim_noise = dim_noise_in * dim_noise_out;
noise_sys = [dim_noise_in, dim_noise_out];

% noise_channel = DephasingChannel(dim_noise_in, 0.9);
noise_channel = DepolarizingChannel(dim_noise_in, 0.9);
% amplitude damping
eps = 0.1;
A{1}=[1 0; 0 sqrt(1-eps)];
A{2}=[0 sqrt(eps); 0 0];
AmplitudeDamping = chanconv(A, 'kraus', 'choi', [2 2]);  %choi matrix of the AD channel
% noise_channel = AmplitudeDamping;


%% simplified SDP for the zero-error sample complexity of shadow communication
cvx_begin sdp quiet
    variable T_pos(dim_noise, dim_noise) hermitian
    variable T_neg(dim_noise, dim_noise) hermitian
    variable V_pos(dim_noise_in, dim_noise_in) hermitian
    variable V_neg(dim_noise_in, dim_noise_in) hermitian
    cost_simplified = trace(V_pos + V_neg) / (dim_target^2);
    minimize cost_simplified
    subject to
        trace(noise_channel.' * (T_pos - T_neg)) == trace(V_pos - V_neg); trace(V_pos - V_neg) == dim_target^2;
        T_pos >= 0; T_neg >= 0; T_pos <= kron(V_pos, eye(dim_noise_out)); T_neg <= kron(V_neg, eye(dim_noise_out));
        PartialTrace(T_pos, 1, noise_sys) == trace(V_pos) * eye(dim_noise_out) / (dim_target^2);
        PartialTrace(T_neg, 1, noise_sys) == trace(V_neg) * eye(dim_noise_out) / (dim_target^2);
cvx_end
cvx_status_simplified = cvx_status


%% dual SDP for the zero-error sample complexity of shadow communication
cvx_begin sdp quiet
    variable lambda
    variable mu
    variable M_pos(dim_noise, dim_noise) hermitian
    variable M_neg(dim_noise, dim_noise) hermitian
    variable N_pos(dim_noise_out, dim_noise_out) hermitian
    variable N_neg(dim_noise_out, dim_noise_out) hermitian
    cost_dual = lambda - mu;
    maximize cost_dual
    subject to
        M_pos >= 0; M_neg >= 0;
        M_pos + dim_target^2 * kron(eye(dim_noise_in), N_pos) >= lambda * noise_channel.';
        M_neg + dim_target^2 * kron(eye(dim_noise_in), N_neg) >= -lambda * noise_channel.';
        PartialTrace(M_pos, 2, noise_sys) == (1 + mu - trace(N_pos)) * eye(dim_noise_in);
        PartialTrace(M_neg, 2, noise_sys) == (1 - mu - trace(N_neg)) * eye(dim_noise_in);
cvx_end
cvx_status_dual = cvx_status


%% output
cost_simplified
cost_dual
