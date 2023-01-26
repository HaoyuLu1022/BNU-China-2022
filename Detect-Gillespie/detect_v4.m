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
par.c8 = 1e-4; % KZW fixed c8*c9 = 3e-5
par.c9 = 0.3; % slow transcription, fast translation
par.c10 = 0.004;
par.c11 = 0.1; 
par.c12 = 0.5;
par.c13 = 6.42e-5;
par.c28 = 0.05; % 0.05; % default 0.05, could be set to 0 for basal expression
par.c29 = 150;  % PmrB~P2 --> PmrB2, equation-free那篇将刺激强度定义为c28/c29
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
par.c30 = 0.04125;
par.Ka = 6.2e7;
% 共32个参数

%% propensity function
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


%% 定义反应初始条件
% VNa set to 1e9
% x(1): RNAP, 30, mol
% x(2): PmrA~P2, mol
% x(3): mRNAa, mol
% x(4): PmrA, mol
% x(5): mRNAb, mol
% x(6): PmrB, mol
% x(7): PmrB~P, mol
% x(8): PmrB2, mol
% x(9): PmrB~P2, mol
% x(10): PmrA·PmrB~P2, mol
% x(11): PmrA~P, mol
% x(12): PmrA~P·PmrB2, mol
% x(13): mRNArep, mol
% x(14): reporter, mol
% x(15): THC, mol/L
init = [30; 0; 0; 3500; 0; 100; 0; 0; 0; 0; 0; 0; 0; 0; 0];
name_elements = ["RNAP", "PmrA~P2", "mRNAa", "PmrA", "mRNAb", "PmrB", "PmrB~P", "PmrB2", "PmrB~P2", "PmrA·PmrB~P2", "PmrA~P", "PmrA~P·PmrB2", "mRNArep", "reporter", "THC"];
num_elements = size(init);
% num_reactions = 29;
% num_reactions = size(detect_stoch, 2);

%% 定义化学计量矩阵
% detect_stoch = [];
detect_stoch = cell2mat(struct2cell(load("Gillespie/detect-stoch-v3.mat")));
% 定义反应通道的时间序列
t1 = 86400; % 21600
tmesh = linspace(0, t1, t1);
num_reactions = size(detect_stoch, 2);

%% 反应通道仿真
traj = SSA(tmesh, par, prop, detect_stoch, init);
% nSim = 10;
% traj = SSAv(tmesh, par, prop, detect_stoch, init, nSim, num_reactions);

% 作图
% figure(1)
% plot(app.UIAxes, traj(14, :), 'b');
% hold(app.UIAxes)
% plot(app.UIAxes, traj(13, :), 'r');
% legend(app.UIAxes, name_elements(14), name_elements(13))
% hold(app.UIAxes)
% plot(app.UIAxes, traj(9, :), 'g');
p1 = plot(tmesh, traj(13, :), 'r');
hold on;
p2 = plot(tmesh, traj(14, :), 'b');
hold on;
p3 = plot(tmesh, traj(9, :), 'g');
xlabel('Time, s')
ylabel('Substance amount, M')
legend(name_elements(13), name_elements(14), name_elements(9), 'location', 'northwest');


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
        for i=1:size(Q)
            if Q(i) < 0
                Q(i) = +0.0; % since propensity functions include square items, when the system lacks something
                % the propensity functions may be positive and bigger and bigger as iteration goes on
            end
        end
        qq = Q(Q > 0);
        Qs = sum(Q);                        %total propensity
        dt = -log(rand())/Qs;               %generating time to the next reaction
        R = sum(rand >= cumsum([0, Q])/Qs); %selecting reaction
        state = state + stoch(:, R);        %updating state
        if ~isempty(find(state < 0, 1))     %in case there are substrates of negative amount
            state = state - stoch(:, R);
            state = state + stoch(:, randperm(numel(qq), 1));
        end
        t = t + dt;                         %updating time
    
        %writing the output
        while cindx <= N && t > tmesh(cindx)
            trajectory(:, cindx) = state;
            cindx = cindx+1;
        end
%         pause(1e-3);
        b(1, [], []);
        if t >= 21600 && t <= 43200 %THC takein time interval
            par.c29 = 0.5; %parameter modified
        end
        if t > 43200
            par.c29 = 150;
        end
    end
%     b.release();
end


function trajectory = SSAv(tmesh, par, prop, stoch, init, nSim, num_reactions)
   %tmesh - time mesh on which solution should be returned
   %par - parameters of the pathway
   %prop - definition of propensity functions
   %stoch - stochiometric matrix
   %init - initial condition for the pathway
   %nSim - number of simulations to perform
 
   tmesh = tmesh(:);   %reshaping mesh to be vertical vector
   t = zeros(nSim, 1);  %current time for each simulation
   state = repmat(init(:)', nSim, 1); %current state variable for each simulation
   trajectory = zeros(nSim, length(init), length(tmesh)); %preparing output trajectory
   trajectory(:, :, 1) = state;%setting initial value as the first element in trajectory
   cindx = 2*ones(nSim, 1);%current trajectory indices
   N = length(tmesh); %number of time points
   aux = 1:nSim; %
 
   while ~isempty(t)
      Q = feval(prop, state, par);         %calculating propensities of the reactions
      Qs = sum(Q, 2);                     %total propensities
      dt = -log(rand(size(Qs, 1), 1))./Qs; %generating time to the next reaction
      P = bsxfun(@rdivide, cumsum([zeros(size(Qs, 1), 1) Q], 2), Qs); %probabilities for each reaction
      R = rem(sum(bsxfun(@ge,rand(size(Qs, 1), 1), P), 2), num_reactions)+1;                %selecting reaction
      state = state + stoch(:, R)';       %updating state
      t = t + dt;                        %updating time
 
     %writing the output
     update = t > tmesh(cindx);
     while any(update)
        %updating state
        iupdt = find(update);
        for i = 1:length(iupdt)
           trajectory(aux(iupdt(i)), :, cindx(iupdt(i))) = state(iupdt(i), :);
        end
        cindx = cindx+update;
 
        %removing finished simulations from the variables
        indx = cindx > N;
        if any(indx)
           cindx(indx) = [];
           t(indx) = [];
           aux(indx) = [];
           state(indx, :) = [];
           if isempty(cindx)
              break;
           end
        end
        update = t > tmesh(cindx);
     end
   end
end



