%%DEFINING REACTION RATES
par.kr = 1/10;     %transcription rate
par.kp = 10;   %translation rate
par.gr = 1/150; %mRNA degradation rate
par.gp = 1/30;  %protein degradation rate

%{
%%DEFINING PROPENSITY FUNCTIONS AS A FUNCTION HANDLE
prop = @(x,par)([par.kr,...      %transcription, one mRNA molecule created
                 par.kp*x(1),... %translation, one protein created
                 par.gr*x(1),... %mRNA degradation, one mRNA molecule removed
                 par.gp*x(2)]);  %protein degradation, one protein molecule removed
%}

%DEFINING PROPENSITY FUNCTIONS AS A FUNCTION HANDLE IN VECTORIZED MANNER
prop = @(x,par)([par.kr*ones(size(x,1),1),...
                 par.kp*x(:,1),...
                 par.gr*x(:,1),...
                 par.gp*x(:,2)]);
 
%%DEFINING INITIAL CONDITION, order [mRNA, Protein]
init = [0;1];
 
%%DEFINING STOICHIOMETRIC MATRIX
%column corresponds to the reaction, row corresponds to the molecule
%order as in prop and init variables
stoch = [1 0 -1 0; 0 1 0 -1];
 
%DEFINING TIME MESH FOR THE OUTPUT TRAJECTORY
tmesh = linspace(0, 1000, 100);

%simulating
nSim = 100;
traj = SSAv(tmesh, par, prop, stoch, init, nSim);
 
%plotting
figure(1)
for i=1:nSim
    tmp = traj(i, 1, :);
    plot(tmesh, tmp(:))
    hold on;
end
xlabel('Time, min')
ylabel('#mRNAs')
 
figure(2)
for i=1:nSim
    tmp = traj(i, 2, :);
    plot(tmesh, tmp(:))
    hold on;
end
xlabel('Time, min')
ylabel('#Proteins')


function trajectory = SSAv(tmesh, par, prop, stoch, init, nSim)
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
      R = sum(bsxfun(@ge,rand(size(Qs, 1), 1), P), 2);                %selecting reaction
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