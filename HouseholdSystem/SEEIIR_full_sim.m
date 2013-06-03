% Simulate the full SEEIIR model in a finite population of households.

% -- want to output the infectious time series to compare to path int
%       results.
% -- modify Josh's SIR code to do this.
% -- alpha is now just scalar, not a vector.
% -- initially stick to one hh size of 4.


%Code for simulating SIR household disease dynamics

%h = pmf of household sizes 1, ..., length(h)

beta = 2; % within-household mixing rate.
sigma = 1;
gamma = 1; % within-household recovery rate.

alpha = 1; % between household mixing rate.

m = 1000; %number of households
%ini = initial number of households with 1 person infected
t_MAX = 150; %(maximum) simulation (epidemic) time
nsims = 1; %number of repeated simulations

%Current State matrix -- each row each household, columns are S,E1,E2,I1,I2
CurrState = zeros(m,5);

max_events = 4*5*m; % max number of possible events.

inf_out = zeros(max_events,1);
time_out = zeros(max_events,1);


%for nind = 1:nsims
    
%     %vector of size of each household
%     HouseholdSizes = mnrnd(1,h,m)';
%     HouseholdSizes = find(HouseholdSizes==1) - length(h)*[0:m-1]';
%     
%     %initial condition
%     CurrState(1:ini,1:2) = [HouseholdSizes(1:ini)-1,ones(ini,1)];
%     CurrState(ini+1:end,1) = HouseholdSizes(ini+1:end);
%     StoreState(nind,1,:) = CurrState(:,2);
%     
%     %size-mixing-biased distribution
%     pi = alpha(HouseholdSizes)'.*HouseholdSizes./sum(alpha(HouseholdSizes)'.*HouseholdSizes);
    
    HouseholdSizes = 4*ones(m,1); % all houses are size 4. 

    total_pop = 4*m;
    
    CurrState(:,1) = 4;
    CurrState(1,1) = 3;
    CurrState(1,2) = 1; % introduce one exposed.

    RateVector = zeros(5*m,1); % introduced this to make the algorthum a bit more efficent.

    %current time
    t = 0;
    count = 1;
    
    % sum(reshape(CurrState(:,2:1:5),4*m,1)) -- total E+I.
    
    %simulate stochastic dynamics
    while t < t_MAX && sum(reshape(CurrState(:,2:1:5),4*m,1)) > 0
        
        %current total infection other households        
%         tif(1) = pi(2:end)'*(CurrState(2:end,2)./HouseholdSizes(2:end));
%         for i = 2:m
%             tif(i) = pi([1:i-1,i+1:end])'*(CurrState([1:i-1,i+1:end],2)./HouseholdSizes([1:i-1,i+1:end]));
%         end

        % infectious prevalence for the population.
        
        tif = sum(sum(CurrState(:,[4,5])))/total_pop;
        
        % count the number of infected HOUSEHOLDS not just individuals.
        
        inf_out(count) = tif*total_pop; %sum(sum(CurrState(:,[4,5]),2)>0);
        time_out(count) = t;
        count = count +1;
        
       
        
        RateVector(1:m) = beta*CurrState(:,1).*sum(CurrState(:,[4,5]),2)./(HouseholdSizes-1) ...
                            + tif*alpha*CurrState(:,1);
                        
        RateVector(m+1:2*m) = 2*sigma*CurrState(:,2);
        RateVector(2*m+1:3*m) = 2*sigma*CurrState(:,3);
        RateVector(3*m+1:4*m) = 2*gamma*CurrState(:,4);
        RateVector(4*m+1:5*m) = 2*gamma*CurrState(:,5);
        
       
        %total rate
        TotRate = sum(RateVector);
        
        %time and event -- time
        t = t + exprnd(1/TotRate);
        

        
        %time and event -- event
        prob = rand;
        cprob = cumsum(RateVector)./TotRate;
        j = 1;
        while prob > cprob(j)
            j = j+1;
        end
        
        
        if j<m+1
            % infection.
            CurrState(j,1:2) = CurrState(j,1:2) + [-1,1];
            
        elseif j< 2*m+1
            % latent prog 1.
            CurrState(j-m,2:3) = CurrState(j-m,2:3) + [-1,+1];
            
        elseif j<3*m+1
            % latent prog 2. Becomes infectious.
            CurrState(j-2*m,3:4) = CurrState(j-2*m,3:4) + [-1,+1];
            
        elseif j<4*m+1
            % infectious prog 1.
            CurrState(j-3*m,4:5) = CurrState(j-3*m,4:5) + [-1,+1];
            
        else
            % recovery
            CurrState(j-4*m,5) = CurrState(j-4*m,5) -1;
          
        end
        
        
    end
    
    % record last event
    inf_out(count) = sum(sum(CurrState(:,[4,5]))); 
    time_out(count) = t;
 
    % plot the number of infectives as a function of time.
    plot(time_out(1:count),inf_out(1:count))

