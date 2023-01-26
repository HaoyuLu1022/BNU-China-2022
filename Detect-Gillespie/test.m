%%DEFINING REACTION RATES
par.kr = 1/10;     %transcription rate
par.kp = 10;   %translation rate
par.gr = 1/150; %mRNA degradation rate
par.gp = 1/30;  %protein degradation rate
 
%%DEFINING PROPENSITY FUNCTIONS AS A FUNCTION HANDLE
prop = @(x,par)([par.kr,...      %transcription, one mRNA molecule created
                 par.kp*x(1),... %translation, one protein created
                 par.gr*x(1),... %mRNA degradation, one mRNA molecule removed
                 par.gp*x(2)]);  %protein degradation, one protein molecule removed
 
%%DEFINING INITIAL CONDITION, order [mRNA, Protein]
init = [0;1];
 
%%DEFINING STOICHIOMETRIC MATRIX
%column corresponds to the reaction, row corresponds to the molecule
%order as in prop and init variables
stoch = [1 0 -1 0; 0 1 0 -1];
 
%DEFINING TIME MESH FOR THE OUTPUT TRAJECTORY
tmesh = linspace(0, 1000, 100);

%simulating
traj = SSA(tmesh, par,prop,stoch, init );
 
%plotting
figure(1)
plot(tmesh, traj(1, :))
xlabel('Time, min')
ylabel('#mRNAs')
 
figure(2)
plot(tmesh, traj(2, :), 'r')
xlabel('Time, min')
ylabel('#Proteins')


function trajectory = SSA(tmesh, par, prop, stoch, init )
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

