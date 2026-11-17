function [t_CPTN, t_TC] = TC_advantage(noise_type, noise_level, H)
%% define some basics
I = eye(2);
JI = MaxEntangled(2,1)*MaxEntangled(2,1)'*2;

if noise_type == 'DE'
    JN = DepolarizingChannel(2, 1 - noise_level) * eye(4);
end

if noise_type == 'AD'
    k0 = [1 0; 0 sqrt(1 - noise_level)];
    k1 = [0 sqrt(noise_level); 0 0];
    JN = kron(I, k0) * JI * kron(I, k0') + kron(I, k1) * JI * kron(I, k1');
end

if noise_type == 'DP'
    JN = DephasingChannel(2, 1 - noise_level) * eye(4);
end


%% information recovery using QPD

cvx_begin sdp quiet
    variable J1(4,4) hermitian % choi matrix of the positive part of the decoding map
    variable J2(4,4) hermitian
    variable p1
    variable p2

    JD = J1 - J2;  % choi matrix of the decoding map
    JPI = PermuteSystems(kron(JI, JD), [2 3 1 4]);
    JF = PartialTrace(JPI*kron(JN.', eye(4)), 1);  % choi matrix of the final channel

    t_CPTN = p1+p2; % cost function (related to the sampling cost)
    minimize t_CPTN
    subject to
        J1 >= 0; PartialTrace(J1, 2) <= p1 * eye(2);
        J2 >= 0; PartialTrace(J2, 2) <= p2 * eye(2);
        ApplyMap(H, PermuteSystems(JF.', [2, 1])) == H;

cvx_end
cvx_status_CPTN = cvx_status;


%% decompose with twisted channel

cvx_begin sdp quiet
    variable J1(4,4) hermitian % choi matrix of the positive part of the decoding map
    variable J2(4,4) hermitian
    variable p

    JD = J1 - J2;  % choi matrix of the decoding map
    JPI = PermuteSystems(kron(JI, JD), [2 3 1 4]);
    JF = PartialTrace(JPI*kron(JN.', eye(4)), 1);  % choi matrix of the final channel

    t_TC = p; % cost function (related to the sampling cost)
    minimize t_TC
    subject to
        J1 >= 0; J2 >= 0;
        PartialTrace(J1+ J2, 2) == p * eye(2);
        ApplyMap(H, PermuteSystems(JF.', [2, 1])) == H;
cvx_end
cvx_status_TC = cvx_status;

fprintf('%s, noise = %.1f: QPD %s; TC %s\n', ...
    noise_type, noise_level, cvx_status_CPTN, cvx_status_TC);
