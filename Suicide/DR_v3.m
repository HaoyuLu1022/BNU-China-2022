import ProgressBar.*

%% 设置系数值与初始矩阵
ka = 0.5;
KA = 10;
a_araC = 1; % estimated
r_ccdB = 0.1333;
d_mc = 0.00203;
b_ccdB = 0.033;
KRA = 0.2;
KRaraC = 0.1;
d_araC = 0.003683; % 0.1
d_ccdB = 0.00115524;
A0 = 1;
X0 = 1;
C0 = 600; % 9.29 500~700
m = 2;
p = 2;
q = 1;
da = 1.03701937823210e-3; %A的扩散率 1.03701937823210e-9 m^2/s -> 1.03701937823210e-3 mm^2/s 9.24

% 网格的大小
width = 100; %128;

% 5,000个模拟秒，每模拟秒4步
dt = 1; %.25;
stoptime = 5000; %5000;
lifespan = 1380;

%% 主程序
% 初始化
[t, A, araC, mRNAccdB, ccdB, cells, life] = initial_conditions(width, stoptime, lifespan);
cells_num = zeros(1, stoptime);
cells_num(1) = sum(sum(cells(:, :, 1)));
threshold = 1e-2;

% 主程序
tic
updateRateHz = 10;
b1 = ProgressBar(stoptime, 'Title', 'Simulating', 'UpdateRate', updateRateHz);
while t < stoptime
    b1(1, [], []);

    A(:, :, t+dt) = A(:, :, t) + (da .* my_laplacian(A(:, :, t)) - ka .* araC(:, :, t) .* power(A(:, :, t), m) ./ (power(KA, m) + power(A(:, :, t), m))) .* dt;
    araC(:, :, t+dt) = araC(:, :, t) + (a_araC ./ (power(A(:, :, t) ./ KRA, p) + 1) - d_araC .* araC(:, :, t)) .* dt;
    mRNAccdB(:, :, t+dt) = mRNAccdB(:, :, t) + (r_ccdB .* C0 ./ (power(araC(:, :, t) ./ KRaraC, q) + 1) .* ((ka .* araC(:, :, t) .* power(A(:, :, t), m) ./ (power(KA, m) + power(A(:, :, t), m))) > 0)- d_mc .* mRNAccdB(:, :, t)) .* dt;
    ccdB(:, :, t+dt) = ccdB(:, :, t) + (b_ccdB .* mRNAccdB(:, :, t) - d_ccdB .* ccdB(:, :, t)) .* dt;
    for i = 1:width
        for j = 1:width % 遍历
            if(life(i, j) == lifespan) % 分裂周期一般是23min/lifespans
                if(ccdB(i, j, t+dt) < threshold) 
                    cells(i, j, t+dt) = cells(i, j, t) * 2;
                    life(i, j) = 1;
                else
                    cells(i, j, t+1800) = 0; % 达到分裂周期但ccdB浓度过高
                    % 可能考虑引入随机因子normrnd(1800, 4)或randi([1700 1800])：效果不稳定，收敛速度很慢
                    life(i, j) = 1;
                end
            else % 没有达到分裂周期
                if(ccdB(i, j, t+dt) < threshold)
                    cells(i, j, t+dt) = cells(i, j, t);
                    life(i, j) = life(i, j) + 1;
                else
                    cells(i, j, t+1800) = 0; % ccdB浓度过高
                    % 可能考虑引入随机因子normrnd(1800, 4)或randi([1700 1800])：效果不稳定，收敛速度很慢
                    life(i, j) = 0; %
                end
            end
        end
    end
%     if(life(:, :) == lifespan)
%         if(ccdB(:, :, t+dt) < threshold)
%             cells(:, :, t+dt) = cells(:, :, t) * 2;
%             life(:, :) = 1;
%         else
%             cells(:, :, t+ceil(normrnd(1800, 25))) = 0;
%             life(:, :) = 1;
%         end
%     else
%         if(ccdB(:, :, t+dt) < threshold)
%             cells(:, :, t+dt) = cells(:, :, t);
%             life(:, :) = life(:, :) + 1;
%         else
%             cells(:, :, t+ceil(normrnd(1800, 25))) = 0;
%             life(:, :) = 1;
%         end
%     end

    cells_num(t+dt) = sum(sum(cells(:, :, t+dt))); % 计算活细胞+未裂解死细胞总数量
%     cells_num(t+dt) = width*width - sum(sum(ccdB(:, :, t+dt) >= threshold)); % 仅计算活细胞
    t = t+dt;
%     if min(min(ccdB(:, :, t))) >= threshold
%         disp(t); % 可以打断点调试看是从哪一秒开始所有网格的ccdB含量都是大于阈值的
%     end
end

%% 画图

% v = VideoWriter("araSuicide-" + num2str(stoptime) + "s.mp4", 'MPEG-4');
% open(v);
% 
% h = figure('Visible','off');
% set(h, 'Units', 'centimeter', 'Position', [5 5 24 16]);
% 
% updateRateHz = 10;
% b2 = ProgressBar(stoptime, 'Title', 'Video generating...', 'UpdateRate', updateRateHz);
% 
% for i = 1:stoptime
%     b2(1, [], []);
% 
%     sgtitle("Calculating "+num2str(i)+"/"+num2str(stoptime)+" s");
%     drawnow;
% 
%     subplot(2, 2, 1);
%     surf(A(:, :, i), 'EdgeColor','none');
%     axis ([0, width, 0, width, 0, 0.1])
%     colorbar;
%     xlabel('x');
%     ylabel('y');
%     zlabel('concentration, M');
%     title('Arabinose');
% 
%     subplot(2, 2, 2);
%     surf(araC(:, :, i), 'EdgeColor','none');
%     axis ([0, width, 0, width])
%     colorbar;
%     xlabel('x');
%     ylabel('y');
%     zlabel('concentration, M');
%     title('araC');
% 
%     subplot(2, 2, 3);
%     surf(mRNAccdB(:, :, i), 'EdgeColor','none');
%     axis ([0, width, 0, width, 0, 10])
%     colorbar;
%     xlabel('x');
%     ylabel('y');
%     zlabel('concentration, M');
%     title('mRNAccdB');
% 
%     subplot(2, 2, 4);
%     surf(ccdB(:, :, i), 'EdgeColor','none');
%     axis ([0, width, 0, width, 0, 150])
%     colorbar;
%     xlabel('x');
%     ylabel('y');
%     zlabel('concentration, M');
%     title('ccdB');
%     
%     frame = getframe(h);
%     writeVideo(v, frame);
% end
% close(v);

h1 = figure;
set(h1, 'Units', 'centimeter', 'Position', [5 5 24 16]);
title('Concentrations of Multiple Substances')

subplot(2, 2, 1);
surf(A(:, :, stoptime), 'EdgeColor','none');
axis ([0, width, 0, width])
colorbar;
xlabel('x');
ylabel('y');
zlabel('concentration, mol/mL');
title('Arabinose');

subplot(2, 2, 2);
surf(araC(:, :, stoptime), 'EdgeColor','none');
axis ([0, width, 0, width])
colorbar;
xlabel('x');
ylabel('y');
zlabel('concentration, mol/mL');
title('araC');

subplot(2, 2, 3);
surf(mRNAccdB(:, :, stoptime), 'EdgeColor','none');
axis ([0, width, 0, width])
colorbar;
xlabel('x');
ylabel('y');
zlabel('concentration, mol/mL');
title('mRNAccdB');

subplot(2, 2, 4);
surf(ccdB(:, :, stoptime), 'EdgeColor','none');
axis ([0, width, 0, width])
colorbar;
xlabel('x');
ylabel('y');
zlabel('concentration, mol/mL');
title('ccdB');

h2 = figure;
plot(cells_num);
title("Live Cells Number by Time");
xlabel('Time, s');
ylabel('Number');
ylim([0, inf])

%% 构造网格影响方程函数：

function out = my_laplacian(in)
%     out = -in ...
%     + .20.*(circshift(in, [1, 0]) + circshift(in, [-1, 0]) + circshift(in, [ 0, 1]) + circshift(in, [0, -1])) ...
%     + .05.*(circshift(in, [1, 1]) + circshift(in, [-1, 1]) + circshift(in, [-1,-1]) + circshift(in, [1, -1]));    
    sz = size(in, 1);
    in1 = circshift(in, [1, 0]);
    in1(1, :) = zeros(1, sz);
    in2 = circshift(in, [-1, 0]);
    in2(sz, :) = zeros(1, sz);
    in3 = circshift(in, [0, 1]);
    in3(:, 1) = zeros(sz, 1);
    in4 = circshift(in, [0, -1]);
    in4(:, sz) = zeros(sz, 1);

    in5 = circshift(in, [1, 1]);
    in5(1, :) = zeros(1, sz);
    in5(:, 1) = zeros(sz, 1);
    in6 = circshift(in, [-1, 1]);
    in6(sz, :) = zeros(1, sz);
    in6(:, 1) = zeros(sz, 1);
    in7 = circshift(in, [-1, -1]);
    in7(sz, :) = zeros(1, sz);
    in7(:, sz) = zeros(sz, 1);
    in8 = circshift(in, [1, -1]);
    in8(1, :) = zeros(1, sz);
    in8(:, sz) = zeros(sz, 1);
    out = -in ...
        + .20 .* (in1 + in2 + in3 + in4) ...
        + .05 .* (in5 + in6 + in7 + in8);
end

%% 设置初始矩阵：

function [t, A, araC, mRNAccdB, ccdB, cells, life] = initial_conditions(n, tmesh, lifespan)
    t = 1;
    A = zeros(n, n, tmesh);
    araC = zeros(n, n, tmesh);
    mRNAccdB = zeros(n, n, tmesh);
    ccdB = zeros(n, n, tmesh);
    cells = ones(n, n, tmesh);
    life = ceil(rand(n, n) * lifespan);

    A(1, :, 1) = 6.6608938919603010724039166056085e-6; %  mmol/mm^3
    % arabinose相对分子质量150.130
    % 9.19: 目前确定是0.05~0.1% g/mL，换算得到6.6608938919603010724039166056085e-6 mol/mL
end