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

% noise_channel = PermuteSystems(kron(DephasingChannel(2, 0.9), DephasingChannel(2, 0.9)), [1 3 2 4]);
noise_channel = PermuteSystems(kron(DepolarizingChannel(2, 0.9), DepolarizingChannel(2, 0.9)), [1 3 2 4]);
% amplitude damping
eps = 0.1;
A{1}=[1 0; 0 sqrt(1-eps)];
A{2}=[0 sqrt(eps); 0 0];
AmplitudeDamping = chanconv(A, 'kraus', 'choi', [2 2]);  %choi matrix of the AD channel
% noise_channel = PermuteSystems(kron(AmplitudeDamping, AmplitudeDamping), [1 3 2 4]);


%% primal SDP for the minimum sampling cost of zero-error channel formation
cvx_begin sdp quiet
    variable R(dim_target, dim_target) hermitian
    variable V_pos(dim_target_out, dim_target_out) hermitian
    variable V_neg(dim_target_out, dim_target_out) hermitian
    cost_primal = trace(V_pos + V_neg) / (dim_id^2);
    minimize cost_primal
    subject to
        trace(V_pos - V_neg) == dim_id^2;
        PartialTrace(R, 2, noise_sys) == trace(V_neg) * eye(dim_target_in) / (dim_id^2);
        R >= 0; noise_channel + R <= kron(eye(dim_target_in), V_pos); R <= kron(eye(dim_target_in), V_neg);

cvx_end
cvx_status_primal = cvx_status;


%% dual SDP for the minimum sampling cost of zero-error channel formation
cvx_begin sdp quiet
    variable M(dim_target, dim_target) hermitian
    variable N(dim_target, dim_target) hermitian
    variable K(dim_target_out, dim_target_out) hermitian
    variable lambda
    cost_dual = real(trace(M * noise_channel)) - lambda;
    maximize cost_dual
    subject to
        M >= 0; N >= 0; M + N + kron(K, eye(dim_target_out)) >= 0;
        dim_id^2 * PartialTrace(M, 1, noise_sys) == (1 + lambda) * eye(dim_target_out);
        dim_id^2 * PartialTrace(N, 1, noise_sys) == (1 - lambda - trace(K)) * eye(dim_target_out);
cvx_end
cvx_status_dual = cvx_status;


%% output
cvx_status_primal
cvx_status_dual
cost_primal
cost_dual
