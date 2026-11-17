clear;
X = [0 1; 1 0];
Z = [1 0; 0 -1];
Y = [0 -1i; 1i 0];
I = eye(2);

channel_noise = 0:0.1:0.5;
H = I + X + Y + Z;

DE_TC = [];
DE_CPTN = [];
AD_TC = [];
AD_CPTN = [];
DP_TC = [];
DP_CPTN = [];


for i = channel_noise
    [t_CPTN, t_TC] = TC_advantage('DE', i, H);
    DE_CPTN(end+1) = t_CPTN;
    DE_TC(end+1) = t_TC;

    [t_CPTN, t_TC] = TC_advantage('AD', i, H);
    AD_CPTN(end+1) = t_CPTN;
    AD_TC(end+1) = t_TC;

    [t_CPTN, t_TC] = TC_advantage('DP', i, H);
    DP_CPTN(end+1) = t_CPTN;
    DP_TC(end+1) = t_TC;

end
