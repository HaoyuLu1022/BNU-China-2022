%% 反应速率：propensity function中的常数
par.c1 = 0.0751;
par.k1 = 9.23e6;
par.k2 = 6.76e15;
par.k3 = 9.23e6;
par.c2 = 1e-4;
par.c3 = 0.1;
par.c4 = 0.004;
par.c5 = 6.42e-5;
par.c6 = 6.42e-5;
par.c7 = 6.42e-5;
par.c8 = 1e-4;
par.c9 = 0.3; % 1e-41太小
par.c10 = 0.004;
par.c11 = 0.1;
par.c12 = 0.5;
par.c13 = 6.42e-5;
par.c28 = 0.05; 
par.c29 = 50;  % PmrB~P2 --> PmrB2 随外界刺激变化 两篇论文分别是0.51, 50 equation-free那篇将刺激强度定义为c28/c29
par.c24 = 1e-4;
par.c25 = 1.6;
par.c26 = 1e-4;
par.c27 = 1.61;
par.c14 = 6.42e-5;
par.c15 = 6.42e-5;
par.c16 = 6.42e-5;
par.c17 = 0.083;
par.c18 = 0.5;
par.c19 = 0.3;
par.c20 = 0.01;
par.c21 = 0.1;
par.c22 = 0.0042;
par.c23 = 6.42e-5;
par.Ka = 6.2e7;
% 共32个参数

%% propensity function
% prop = @(x,par)([par.c1*par.k2*x(1)*x(2)/(1+par.k2*x(1)*x(2)+par.k1*x(1)+par.k3*x(2)) ...
%     par.c2*par.k1*x(1)/(1+par.k2*x(1)*x(2)+par.k1*x(1)+par.k3*x(2)) ...
%     par.c3*x(3) ...
%     par.c4*x(3) ...
%     par.c5*x(4) ...
%     par.c6*x(11) ...
%     par.c7*x(11) ...
%     par.c8 ...
%     par.c9*x(5) ...
%     par.c10*x(5) ...
%     par.c13*x(6) ...
%     par.c11*x(6) ...
%     par.c12*x(6) ...
%     par.c28*x(6) ...
%     par.c29*x(9) ...
%     par.c24*x(4)*x(9) ...
%     par.c25*x(10) ...
%     par.c26*x(11)*x(6) ...
%     par.c27*x(12) ...
%     par.c14*x(7) ...
%     par.c15*x(8) ...
%     par.c16*x(9) ...
%     par.c17*x(11) ...
%     par.c18*x(2) ...
%     par.c19*par.k2*x(1)*x(2)/(1+par.k2*x(1)*x(2)+par.k1*x(1)+par.k3*x(2)) ...
%     par.c20*par.k1*x(1)/(1+par.k2*x(1)*x(2)+par.k1*x(1)+par.k3*x(2)) ...
%     par.c21*x(13) ...
%     par.c22*x(13) ...
%     par.c23*x(14) ...
%     ]);

prop = @(x,par)([par.c1.*par.k2.*x(:,1).*x(:,2)/(1+par.k2.*x(:,1).*x(:,2)+par.k1.*x(:,1)+par.k3.*x(:,2)) ...
    par.c2.*par.k1.*x(:,1)/(1+par.k2.*x(:,1).*x(:,2)+par.k1.*x(:,1)+par.k3.*x(:,2)) ...
    par.c3.*x(:,3) ...
    par.c4.*x(:,3) ...
    par.c5.*x(:,4) ...
    par.c6.*x(:,11) ...
    par.c7.*x(:,11) ...
    par.c8.*ones(size(x, 1), 1) ...
    par.c9.*x(:,5) ...
    par.c10.*x(:,5) ...
    par.c13.*x(:,6) ...
    par.c11.*x(:,6).*(x(:,6)-1) ...
    par.c12.*x(:,6) ...
    par.c28.*x(:,6) ...
    par.c29.*x(:,9) ...
    par.c24.*x(:,4).*x(:,9) ...
    par.c25.*x(:,10) ...
    par.c26.*x(:,11).*x(:,6) ...
    par.c27.*x(:,12) ...
    par.c14.*x(:,7) ...
    par.c15.*x(:,8) ...
    par.c16.*x(:,9) ...
    par.c17.*x(:,11).*(x(:, 11)-1) ...
    par.c18.*x(:,2) ...
    par.c19.*par.k2.*x(:,1).*x(:,2)/(1+par.k2.*x(:,1).*x(:,2)+par.k1.*x(:,1)+par.k3.*x(:,2)) ...
    par.c20.*par.k1.*x(:,1)/(1+par.k2.*x(:,1).*x(:,2)+par.k1.*x(:,1)+par.k3.*x(:,2)) ...
    par.c21.*x(:,13) ...
    par.c22.*x(:,13) ...
    par.c23.*x(:,14) ...
]);


%% 定义反应初始条件
% [RNAP PmrA~P2 mRNAa PmrA mRNAb PmrB PmrB~P PmrB2 PmrB~P2 PmrA.PmrB~P PmrA~P PmrA~P.PmrB2 mRNArep reporter]
init = [30; 0; 0; 3800; 0; 25; 0; 0; 0; 0; 0; 0; 0; 0];
init = gpuArray(init);
name_elements = ["RNAP", "PmrA~P2", "mRNAa", "PmrA", "mRNAb", "PmrB", "PmrB~P", "PmrB2", "PmrB~P2", "PmrA·PmrB~P2", "PmrA~P", "PmrA~P·PmrB2", "mRNArep", "reporter"];
num_elements = size(init);
num_reactions = size(detect_stoch, 2);

%% 定义化学计量矩阵
% detect_stoch = [];
detect_stoch = gpuArray(detect_stoch);

% 定义反应通道的时间序列
tmesh = gpuArray(linspace(0, 200, 200));
% tmesh = linspace(0, 50);

%% 反应通道仿真
% traj = SSA(tmesh, par, prop, detect_stoch, init);
nSim = 2;
traj = SSAv(tmesh, par, prop, detect_stoch, init, nSim, num_reactions);

% 作图
figure(1)
% plot(tmesh, traj(13, :), 'r')
for i = 1:nSim
    tmp = traj(i, 13, :);
    plot(tmesh, tmp(:), 'r'); % mRNA
    hold on;
end
hold on;
% plot(tmesh, traj(14, :), 'b')
for i = 1:nSim
    tmp = traj(i, 14, :);
    plot(tmesh, tmp(:), 'b'); % Protein
    hold on;
end
xlabel('Time, s')
ylabel('#Substance amount, M')
% legend("mRNArep", "Reporter")

%{
figure(1)
for i = 1:num_elements
    plot(tmesh, traj(i, :))
%     leg{i} = name_elements[i];
    hold on;
end
xlabel('Time, s')
ylabel('#Substances, M')
% legend(["RNAP", "PmrA~P2", "mRNAa", "PmrA", "mRNAb", "PmrB", "PmrB~P", "PmrB2", "PmrB~P2", "PmrA·PmrB~P2", "PmrA~P", "PmrA~P·PmrB2", "mRNArep", "reporter"]);
%}



%% Gillespie算法：朴素版 & vector优化版
function trajectory = SSA(tmesh, par, prop, stoch, init)
   %tmesh - time mesh on which solution should be returned
   %per - parameters of the pathway
   %prop - definition of propensity functions
   %stoch - stochiometric matrix
   %init - initial condition for the pathway
 
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
      R = sum(rand >= cumsum([0 Q])/Qs);  %selecting reaction
      state = state + stoch(:, R);        %updating state
      t = t + dt;                         %updating time
 
      %writing the output
      while cindx <= N && t > tmesh(cindx)
         trajectory(:, cindx) = state;
         cindx = cindx+1;
     end
   end
end


function trajectory = SSAv(tmesh, par, prop, stoch, init, nSim, num_reactions)
   %tmesh - time mesh on which solution should be returned
   %par - parameters of the pathway
   %prop - definition of propensity functions
   %stoch - stochiometric matrix
   %init - initial condition for the pathway
   %nSim - number of simulations to perform
 
   tmesh = tmesh(:);   %reshaping mesh to be vertical vector
   t = gpuArray(zeros(nSim, 1));  %current time for each simulation
   state = repmat(init(:)', nSim, 1); %current state variable for each simulation
   trajectory = gpuArray(zeros(nSim, length(init), length(tmesh))); %preparing output trajectory
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





