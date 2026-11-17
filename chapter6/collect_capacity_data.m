clear;

epsilon = 0.0;
capacity_ico_list = zeros(1, 51);
capacity_def_list = zeros(1, 51);
for index = 0:50
    damp_rate = index / 100;
    classical_capacity_ico;
    capacity_ico_list(index+1) = capacity_ico;
    capacity_def_list(index+1) = capacity_def;
end
save("capacity_data.mat", "capacity_ico_list", "capacity_def_list");

epsilon = 0.02;
for index = 0:50
    damp_rate = index / 100;
    classical_capacity_ico;
    capacity_ico_list(index+1) = capacity_ico;
    capacity_def_list(index+1) = capacity_def;
end
save("capacity_with_err_data.mat", "capacity_ico_list", "capacity_def_list");
