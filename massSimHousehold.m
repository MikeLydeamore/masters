gamma=2; beta=500; k=10;

maxIter=100;
timeCount=cell(maxIter,1);
infCount=cell(maxIter,1);
t=0:0.01:5;
infectedMatrix=zeros(1,length(t));
infectedMatrix(1)=maxIter; %Always start with 1 infected.
for i=1:maxIter %Simulate
    [infNo times]=simHousehold(gamma,beta,k);
    if times(end) > max(t) %Extend max time if necessary
        disp(times);
        infectedMatrix=[infectedMatrix zeros(1,length(t(end):0.01:times(end)))];
        t=[t t(end):0.01:times(end)];
    end
    
   %Want to know the number of infected at time t
   for i=2:length(t)
       transitionAt=find(t(i)>times);
       transitionAt=transitionAt(end);
       infectedMatrix(transitionAt)=infectedMatrix(transitionAt)+infNo(transitionAt);
   end
end
%Normalise
propInfected=infectedMatrix./maxIter;
plot(t,propInfected);
        