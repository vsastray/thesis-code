clear;

noise_level = 0.2;

X = [0 1; 1 0];
Z = [1 0; 0 -1];
Y = [0 1i; -1i 0];
I = eye(2);

% depolarizing channel
JN = DepolarizingChannel(2, 1 - noise_level) * eye(4);

% observable
H = I + X + Y + Z;

% SDP of information recovering using TC

cvx_begin sdp quiet
    variable J1(4,4) hermitian % choi matrix of the positive part of the decoding map
    variable J2(4,4) hermitian
    variable overhead

    % choi matrix of the decoding map
    JD = J1 - J2;
    % choi matrix of the final channel
    JF = PartialTrace(kron(PartialTranspose(JN, 2, [2,2]), eye(2)) * kron(eye(2), JD), 2, [2, 2, 2]);

    minimize overhead
    subject to
        J1 >= 0;
        J2 >= 0;
        PartialTrace(J1 + J2, 2) == overhead * eye(2);
        ApplyMap(H, PermuteSystems(JF.', [2, 1])) == H;
cvx_end

J1 = J1 / overhead;
J2 = J2 / overhead;

%% sampling times
rng(13)
state = RandomDensityMatrix(2);
max_sample_times = 6003;
repeat_time = 300;

sample_ratios = [0.2 0.4 0.6 0.8 1.0];
sample_data = zeros(numel(sample_ratios), repeat_time);
for k = 1:numel(sample_ratios)
    N_sample = round(sample_ratios(k) * max_sample_times);

    for i = 1:repeat_time
        [tar, est] = sampling(state, N_sample, JN, J1, J2, overhead, H);
        sample_data(k, i) = est;
    end

    fprintf('Ratio %.1f: mean = %.12g, std = %.12g\n', ...
        sample_ratios(k), mean(sample_data(k, :)), ...
        std(sample_data(k, :)));

    filename = sprintf('sample_%.1f.csv', sample_ratios(k));
    csvwrite(filename, sample_data(k, :));
end
