beta=2;
gamma=2;
alpha=1;

N=10;
houseSize=4;

popStatus=zeros(N,3); %S-I-R
popStatus(2:end,1)=N*ones(N-1,1);
popStatus(1,1)=N-1; %1 infected
popStatus(1,2)=1;

max_events=3*houseSize*N;

noInfected=zeros(1,max_events);
evenTime=zeros(1,max_events);