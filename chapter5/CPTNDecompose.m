function [J1, J2, p1, p2, cvx_status] = CPTNDecompose(JF, dim_list)
    % compute a minimum-cost QPD decomposition of the Choi matrix JF
    % dim_list = [input_dim, output_dim]
    % the returned CPTN Choi matrices satisfy JF = p1 * J1 - p2 * J2
    % the QPD cost is p1 + p2

    dim = numel(JF(:, 1));

    cvx_begin sdp quiet
    variable J1(dim, dim) complex semidefinite
    variable J2(dim, dim) complex semidefinite
    variable p1 nonnegative
    variable p2 nonnegative

    minimize p1 + p2
    subject to

        PartialTrace(J1, 2, dim_list) <= p1 * eye(dim_list(1));
        PartialTrace(J2, 2, dim_list) <= p2 * eye(dim_list(1));
        J1 - J2 == JF;

    cvx_end

    J1 = full(J1 / p1);
    J2 = full(J2 / p2);
end
