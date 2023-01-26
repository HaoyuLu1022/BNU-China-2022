import ProgressBar.*

%% 反应速率：propensity function中的常数
par.c1 = 0.0751;
par.k1 = 9.23e-3;
par.k2 = 6.76e-3;
par.k3 = 9.23e-3;
par.c2 = 1e-4;
par.c3 = 0.1;
par.c4 = 0.004;
par.c5 = 6.42e-5;
par.c6 = 6.42e-5;
par.c7 = 6.42e-5;
par.c8 = 1e-3; % KZW fixed c8*c9 = 3e-5
par.c9 = 3e-2; % slow transcription, fast translation
par.c10 = 0.004;
par.c11 = 0.1; 
par.c12 = 0.5;
par.c13 = 6.42e-5;
par.c28 = 0.05; % default 0.05, could be set to 0 for basal expression
par.c29 = 7.9;  % PmrB~P2 --> PmrB2, equation-free那篇将刺激强度定义为c28/c29
% c29 < 7.9, only activated; c29 > 110.4, only basal; in between, both coexists
par.c24 = 1e-4;
par.c25 = 1.6;
par.c26 = 1e-4;
par.c27 = 1.6;
par.c14 = 6.42e-5;
par.c15 = 6.42e-5;
par.c16 = 6.42e-5;
par.c17 = 0.083;
par.c18 = 0.5;
par.c19 = 0.3;
par.c20 = 0.01;
par.c21 = 0.1;
par.c22 = 0.004;
par.c23 = 6.202e-4; % 6.42e-5
par.c24 = 0.04125;
par.Ka = 6.2e7;

prop = @(x,par)([par.c1*par.k2*x(1)*x(2)/(1+par.k2*x(1)*x(2)+par.k1*x(1)+par.k3*x(2)) ...
    par.c2*par.k1*x(1)/(1+par.k2*x(1)*x(2)+par.k1*x(1)+par.k3*x(2)) ...
    par.c3*x(3) ...
    par.c4*x(3) ...
    par.c5*x(4) ...
    par.c6*x(11) ...
    par.c7*x(11) ...
    par.c8 ...
    par.c9*x(5) ...
    par.c10*x(5) ...
    par.c13*x(6) ...
    par.c11*x(6)*(x(6)-1) ...
    par.c12*x(6) ...
    par.c28*x(6) ...
    par.c29*x(9) ...
    par.c24*x(4)*x(9) ...
    par.c25*x(10) ...
    par.c26*x(11)*x(6) ...
    par.c27*x(12) ...
    par.c14*x(7) ...
    par.c15*x(8) ...
    par.c16*x(9) ...
    par.c17*x(11)*(x(11)-1) ...
    par.c18*x(2) ...
    par.c19*par.k2*x(1)*x(2)/(1+par.k2*x(1)*x(2)+par.k1*x(1)+par.k3*x(2)) ...
    par.c20*par.k1*x(1)/(1+par.k2*x(1)*x(2)+par.k1*x(1)+par.k3*x(2)) ...
    par.c21*x(13) ...
    par.c22*x(13) ...
    par.c23*x(14) ...
    ]);

init = traj(:, t1); % 接着traj的末尾
init(15) = 0.1; % THC浓度发生变化
name_elements = ["RNAP", "PmrA~P2", "mRNAa", "PmrA", "mRNAb", "PmrB", "PmrB~P", "PmrB2", "PmrB~P2", "PmrA·PmrB~P2", "PmrA~P", "PmrA~P·PmrB2", "mRNArep", "reporter", "THC"];
num_elements = size(init);

detect_stoch = cell2mat(struct2cell(load("Gillespie/detect-stoch-v3.mat")));
num_reactions = size(detect_stoch, 2);

t = 10800;
tmesh = linspace(t1, t1+t, t);

traj1 = SSA(tmesh, par, prop, detect_stoch, init);
% 然后使用traj2 = cat(2, traj, traj1)拼接两个矩阵

figure(2)
plot(tmesh, traj1(13, :), 'r')
hold on;
plot(tmesh, traj1(14, :), 'b')
hold on;
plot(tmesh, traj1(9, :), 'g')
xlabel('Time, s')
ylabel('#Substance amount, M')
legend("mRNArep", "Reporter", "PmrB~P2")

%% Gillespie算法：朴素版 & vector优化版
function trajectory = SSA(tmesh, par, prop, stoch, init)
    %tmesh - time mesh on which solution should be returned
    %per - parameters of the pathway
    %prop - definition of propensity functions
    %stoch - stochiometric matrix
    %init - initial condition for the pathway
    
    updateRateHz = 10;
    b = ProgressBar([], 'Title', 'Simulating', 'UpdateRate', updateRateHz);
    t = 0;                                           %current time
    state = init(:);                                 %variable with current system state
    trajectory = zeros(length(init), length(tmesh)); %preparing output trajectory
    trajectory(:, 1) = init(:);                      %setting initial value as the first element in trajectory
    cindx = 2;                                       %current trajectory index
    N = length(tmesh);                               %number of time points
    
    while t < tmesh(end)
        Q = feval(prop, state, par);        %calculating propensities of the reactions
        Qs = sum(Q);                        %total propensity
        dt = -log(rand())/Qs;               %generating time to the next reaction
        R = sum(rand > cumsum([0 Q])/Qs);  %selecting reaction
        state = state + stoch(:, R);        %updating state
        t = t + dt;                         %updating time
    
        %writing the output
        while cindx <= N && t > tmesh(cindx)
            trajectory(:, cindx) = state;
            cindx = cindx+1;
        end
%         pause(1e-3);
        b(1, [], []);
    end
%     b.release();
end