function ydot = K_by_T(t, y)
% 参数 %
D = 789; % 236;
k_M = 0.03;
b_1 = 0.009;
b_2 = 0.0317;
a = 0.000365;
d_m = 0.002;
d_c = 0.00028;
d_a = 0.00231;
d_a2 = 0.000231;
d_large = 0.198;
beta = 0.1;
K_t = 90; % MazF cannot reach the threshold of K_t
% phiT_MazE = 1; % T > 37 
phiT_MazE = 0.132; % T <= 37 
phiT_MazF = 1;

% 函数命名 %
% ydot(1) = d[mRNA_MazF] / dt %
% ydot(2) = d[mRNA_MazE] / dt %
% ydot(3) = d[MazF] / dt % 
% yodt(4) = d[MazE] / dt %
% yodt(5) = d[MazEF] / dt %

ydot = zeros(5, 1);

ydot(1) = D * k_M * phiT_MazF - (d_large * (beta * y(3)) ^ 2) / ((beta * y(3)) ^ 2 + K_t ^ 2) * y(1) - d_m * y(1);
ydot(2) = D * k_M * phiT_MazE - (d_large * (beta * y(3)) ^ 2) / ((beta * y(3)) ^ 2 + K_t ^ 2) * y(2) - d_m * y(2);
ydot(3) = b_1 * y(1) - a * y(3) * y(4) - d_c * y(3) + d_a2 * y(5);
ydot(4) = b_2 * y(2) - a * y(3) * y(4) - d_a * y(4);
ydot(5) = a * y(3) * y(4) - d_c * y(5) - d_a2 * y(5);

end

















