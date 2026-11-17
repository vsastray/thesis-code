clear;

case_id = 2; % examples 1–7, from left to right in Figure 5.4
input_dim = 6;

% upper-triangular input indices; symmetric entries are added below
selected_pairs = {
    [1 1; 1 2; 1 4; 2 2; 2 3; 2 4; 3 3; 3 4; 3 5; ...
     4 4; 4 5; 5 5];
    [3 5; 5 5];
    [4 4; 4 6];
    [2 3; 3 3; 3 5];
    [2 5; 4 5; 5 6; 6 6];
    [1 2; 1 3; 2 2; 3 3];
    [2 2; 2 4; 3 4; 3 5; 4 5]
};

pairs = selected_pairs{case_id};

% selected input indices and retained output-index pairs
I = unique(pairs(:)).';
[~, A] = ismember(pairs, I);
output_dim = numel(I);

JF = choi_convolution(input_dim, I, A);
dim_list = [input_dim output_dim];

[J1, J2, cost_TC, cvx_status_TC] = TCDecompose(JF, dim_list);
assert(norm(JF - (J1 - J2), 'fro') < 1e-6);

[J1, J2, p1, p2, cvx_status_QPD] = CPTNDecompose(JF, dim_list);
assert(norm(JF - (p1 * J1 - p2 * J2), 'fro') < 1e-6);
cost_QPD = p1 + p2;

fprintf('Example %d\n', case_id);
fprintf('TC cost  = %.12g\n', cost_TC);
fprintf('QPD cost = %.12g\n', cost_QPD);

function choi = choi_convolution(orig_dim, index_map, feature)
    output_dim = numel(index_map);
    feature_adj = feature(:, [2 1]);
    index_set = unique([feature; feature_adj], 'rows');

    choi = zeros(orig_dim * output_dim);

    for n = 1:size(index_set, 1)
        j = index_set(n, 1);
        k = index_set(n, 2);

        mat_1 = zeros(orig_dim, orig_dim);
        mat_1(index_map(j), index_map(k)) = 1;

        mat_2 = zeros(output_dim, output_dim);
        mat_2(j, k) = 1;

        choi = choi + kron(mat_1, mat_2);
    end
end
