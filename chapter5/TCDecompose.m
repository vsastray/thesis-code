function [J1, J2, a, cvx_status] = TCDecompose(JF, dim_list)
    % compute a minimum-cost twisted-channel decomposition of JF
    % dim_list = [input_dim, output_dim]
    % JF = J1 - J2, with J1 and J2 positive semidefinite
    % (J1 + J2) / a is the Choi matrix of a CPTP map
    % the twisted-channel cost is a

    dim = numel(JF(:, 1));

    cvx_begin sdp quiet
    variable J1(dim, dim) complex semidefinite
    variable J2(dim, dim) complex semidefinite
    variable a nonnegative

    minimize a
    subject to
        PartialTrace(J1 + J2, 2, dim_list) == a * eye(dim_list(1));
        J1 - J2 == JF;

    cvx_end

    J1 = full(J1);
    J2 = full(J2);
end
