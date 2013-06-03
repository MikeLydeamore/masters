Last update: 8/3

Current files:
---expecGamma---
Solves the system of linear equations to give E[Gamma|X(0)=i] for any given Q and stateList and outer-household infection rate alpha


--genQ---
Generates the Q matrix for the infection and recovery of an SIR model.
Should be extendable for the SEIR model and also for more intricate models.
Is not super naive, but there is probably still room for efficiency.