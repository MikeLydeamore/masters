function [Q stateList]=genQ(beta,gamma,sigma,ro,tau,zeta,kappa,ita,k,absorbingOnIntro)
%Generates the Q matrix for a single household of size k including
%antiviral introduction with random exponential delay.

%The idea is to have a 'flip' state which will take us to another section
%of the Markov chain where antiviral rates are in effect. This should give
%our Q matrix the form
%     [ Q1 | I  |  0  ]
%     [  0 | Q2 |  KI ]
%     [  0 | 0  |  Q1 ]
%Arbitrarily give the delay a parameter 1
%
%Q1 should be exactly the same as our regular Q matrix from the model
%without antivirals
%
%Q2 will have the modified rates
%
%We essentially then have 2 'copies' of our state list, (S,E,I), one with
%antivirals and one without.
%This may make calculating the number of infected hard, will have to think
%more about this.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Upper left block                               %
%                          Centre block                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Our 2 blocks share exactly the same structure, but have different rates,
%so we might as well generate them simulateously.

%Generate state list
counter=1;
for s=1:k+1
    for e=1:k+1
        for i=1:k+1
            if (s-1)+(e-1)+(i-1)<=k
                stateList(:,counter)=[s-1;e-1;i-1];
                counter=counter+1;
            end
        end
    end
end
[m, noStates]=size(stateList);
Q1=sparse(noStates,noStates);
Q2=sparse(noStates,noStates);

%Infection events (S>0, I>0)
infectionEvents=intersect(find(stateList(1,:)>0),find(stateList(3,:)>0));

%Move to state (S-1, E+1, I) at rate beta*S*I/(k-1)
%For Q2, the rate is (1-ro)(1-tau)*beta*S*I/(k-1)
for i=1:length(infectionEvents)
    currState=infectionEvents(i);
    sminusone = find(stateList(1,:)==stateList(1,currState)-1);
    eplusone = find(stateList(2,:)==stateList(2,currState)+1);
    currenti = find(stateList(3,:)==stateList(3,currState));
    moveTo=intersect(intersect(sminusone,eplusone),currenti);
    Q1(currState,moveTo)=beta*stateList(1,currState)*stateList(3,currState)/(k-1);
    Q2(currState,moveTo)=(1-ro)*(1-tau)*beta*stateList(1,currState)*stateList(3,currState)/(k-1);
end

%Shedding events (E>0)
sheddingEvents=find(stateList(2,:)>0);

%Move to state (S, E-1, I+1) at rate sigma*E for both Q1 and Q2
for i=1:length(sheddingEvents)
    currState=sheddingEvents(i);
    currents = find(stateList(1,:)==stateList(1,currState));
    eminusone = find(stateList(2,:)==stateList(2,currState)-1);
    iplusone = find(stateList(3,:)==stateList(3,currState)+1);
    moveTo=intersect(intersect(currents,eminusone),iplusone);
    Q1(currState,moveTo)=sigma*stateList(2,currState);
    Q2(currState,moveTo)=sigma*stateList(2,currState);
end

%Recovery events (I>0)
recoveryEvents=find(stateList(3,:)>0);

%Move to state (S, E, I-1) at rate gamma*I for Q1 and rate (1+ita)*gamma*I
%for Q2.
for i=1:length(recoveryEvents)
    currState=recoveryEvents(i);
    currents = find(stateList(1,:)==stateList(1,currState));
    currente = find(stateList(2,:)==stateList(2,currState));
    iminusone = find(stateList(3,:)==stateList(3,currState)-1);
    moveTo=intersect(intersect(currents,currente),iminusone);
    Q1(currState,moveTo)=gamma*stateList(3,currState);
    Q2(currState,moveTo)=(1+ita).*gamma.*stateList(3,currState);
end

%Leave the diagonals until the last step.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Switch Block                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assume we always start with 1 infected (otherwise there is no epidemic) so
%initiate the arrival of antivirals.
%Need to remove the initial state of only having exposed (k-1, 1, 0, 0)
switchblock = zeta*eye(noStates);
initialState=find(stateList(1,:)==(k-1)&stateList(2,:)==1&stateList(3,:)==0);
switchblock(initialState,initialState)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Centre-Right Block                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switchbackblock=kappa*eye(noStates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Generate Q                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zeroMat=zeros(size(Q1));
if absorbingOnIntro==0
    Q=[Q1, switchblock zeroMat ; zeroMat, Q2 switchbackblock ; zeroMat, zeroMat, Q1];
elseif absorbingOnIntro==-1
    Q=[zeroMat, zeroMat, zeroMat; zeroMat, Q2, switchbackblock ; zeroMat, zeroMat, zeroMat];
else
    Q=[Q1, switchblock, zeroMat ; zeroMat, zeroMat, zeroMat ; zeroMat, zeroMat, zeroMat];
end
%Add in diagonals
rowSum=sum(Q,2);
[m, noStates]=size(Q);
for i=1:noStates
    Q(i,i)=-rowSum(i);
end

stateList=[stateList stateList stateList; zeros(1,length(stateList)), ones(1,length(stateList)), 2*ones(1,length(stateList)) ];

end