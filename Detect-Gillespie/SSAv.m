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