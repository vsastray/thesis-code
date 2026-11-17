% inputs: damp_rate and epsilon, set in the Workspace before running

% amplitude damping
A{1}=[1 0; 0 sqrt(1-damp_rate)];
A{2}=[0 sqrt(damp_rate); 0 0];
AmplitudeDamping = chanconv(A, 'kraus', 'choi', [2 2]);  %choi matrix of the AD channel
init_chan = tensor(AmplitudeDamping, AmplitudeDamping);

dim_sys = [2 2 2 2];


%% ico capacity SDP
cvx_begin sdp quiet
    variable M(16, 16) hermitian % A_0 A_1 B_0 B_1
    variable N(16, 16) hermitian % A_0 A_1 B_0 B_1
    variable m
    % useful definitions part for N
    b1_N = tensor(PartialTrace(N, 4, dim_sys), eye(2) / 2);
    b0b1_N = tensor(PartialTrace(b1_N, 3, dim_sys), eye(2) / 2);
    bpart_N = N - b1_N + b0b1_N;
    a1bpart_N = PermuteSystems(tensor(PartialTrace(bpart_N, 2, dim_sys), eye(2) / 2), [1 4 2 3], dim_sys);
    a0a1bpart_N = tensor(eye(2) / 2, PartialTrace(a1bpart_N, 1, dim_sys));
    N_ns = bpart_N - a1bpart_N + a0a1bpart_N;

    maximize m
    subject to
        real(trace(init_chan * M)) >= m * (1 - epsilon);
        M >= 0; N >= M; PartialTrace(M, [1 3], dim_sys) == eye(4);
        N_ns == m * eye(16) / 4;
cvx_end
capacity_ico = m;
cvx_status_ico = cvx_status;


%% def capacity SDP
cvx_begin sdp quiet
    variable E(16, 16) hermitian % A_0 A_1 B_0 B_1
    variable F(16, 16) hermitian % A_0 A_1 B_0 B_1
    variable G(16, 16) hermitian % A_0 A_1 B_0 B_1
    variable H(16, 16) hermitian % A_0 A_1 B_0 B_1
    variable m
    % useful definitions part for F and H
    b1_F = tensor(PartialTrace(F, 4, dim_sys), eye(2) / 2);
    F_a0a1 = PartialTrace(F, 2, [4 4]);
    a1_F_a0a1 = tensor(PartialTrace(F_a0a1, 2, [2 2]), eye(2) / 2);
    a1_H = PermuteSystems(tensor(PartialTrace(H, 2, dim_sys), eye(2) / 2), [1 4 2 3], dim_sys);
    H_b0b1 = PartialTrace(H, 1, [4 4]);
    b1_H_b0b1 = tensor(PartialTrace(H_b0b1, 2, [2 2]), eye(2) / 2);

    maximize m
    subject to
        real(trace(init_chan * (E + G))) >= m * (1 - epsilon);
        E >= 0; F >= E; G >= 0; H >= G;
        PartialTrace(E + G, [1 3], dim_sys) == eye(4);
        F == b1_F; F_a0a1 == a1_F_a0a1; H == a1_H; H_b0b1 == b1_H_b0b1; trace(F + H) == m * 4;
cvx_end
capacity_def = m;
cvx_status_def = cvx_status;
