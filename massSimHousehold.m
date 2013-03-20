function massSimHousehold(beta,gamma,k)
%gamma=2; beta=5; k=10;

maxIter=400;
timeCount=cell(maxIter,1);
infCount=cell(maxIter,1);
t=0:0.001:5;
infectedMatrix=zeros(1,length(t));
%infectedMatrix(1)=maxIter; %Always start with 1 infected.
for i=1:maxIter %Simulate
    [infNo times]=simHousehold(gamma,beta,k);
    if times(end) > max(t) %Extend max time if necessary
        infectedMatrix=[infectedMatrix zeros(1,length(t(end):0.01:times(end)))];
        t=[t t(end):0.01:times(end)];
    end
    
    %Want to know the number of infected at time t
    for ic=1:length(t)
        increase=find(times<=t(ic));
        increase=increase(end); %Finds the most recent event for the current point t(ic)
        
        infectedMatrix(ic)=infectedMatrix(ic)+infNo(increase);
        
    end
end
%Normalise
propInfected=infectedMatrix./maxIter;
plot(t,propInfected);
end
