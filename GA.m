% Programmer : Mj.Yazdani
% Email : Mj.Yazdani1988@gmail.com
% Date: January 29, 2016 

function []  =GA()
VarMin=-5.12; % Minimum
VarMax=5.12; % Maximum
PreferedBestCost=0;

nVar=5;  % number of genes
nPop=10; % number of initiate poplulation
MaxIt=600; % number of iteration
pc=0.8; % crossover rate
nc=2*round(pc*nPop/2); % number of crossover
gamma=0.05; % rate of genes that should be crossovered
pm=3; % mutation rate
nm=round(pm*nPop); % number of mutation
mu=0.02; % mutation rate in every chromosome

beta=8; % selection pressure

% initialization
pop=iCreatePopulation(VarMin, VarMax, nVar,nPop);
tmp=zeros(1,nPop);
for i=1: nPop
    tmp(i)=iFittness(pop(i,:));
end
pop=[pop,tmp'];

% sort Population
pop= sortrows(pop,nVar+1);

% Store Best and Worst Cost
%BestCost=pop(1,nVar+1);
WorstCost=pop(nPop,nVar+1);

mm=0;
CostAvgInEachStep=zeros(1,MaxIt);

% Main Loop
while (mm<=MaxIt)
    size(pop)
    % Calculate selection propbabilities
    P=exp(-beta*pop(:,nVar+1)/WorstCost);
    P=P/sum(P);

    % Crossover
    %popc=zeros(nVar,nPop);
    for k=1:nc/2
        i1=iRouleteWheelSelection(P);
        i2=iRouleteWheelSelection(P);

        % Select Parents
        p1=pop(i1,:);
        p2=pop(i2,:);

        % Apply  Crossover
        if k==1
            popc=[ iCrossover(p1,p2,gamma,VarMin,VarMax)];
        else
            popc=[popc ; iCrossover(p1,p2,gamma,VarMin,VarMax)];
        end
    end
    % evaluate cost function for crossover population
    tmp=zeros(1,nc);
    for i= 1 : nc
        tmp(i)=iFittness(popc(i,:));
    end
    %%%%%%popc=[popc,tmp'];
    % Mutation
    for k=1:nm
        i=randi([1 nPop]);
        if k==1
            popm=[ iMutation(pop(i,:),mu,VarMin,VarMax)];
        else
            popm=[popm ; iMutation(pop(i,:),mu,VarMin,VarMax)];
        end
    end

    %evaluate cost function for mutation population
    tmp=zeros(1,pm);
    for i= 1 : pm
        tmp(i)=iFittness(popm(i,:));
    end
    %popm=[popm,tmp'];

    %create merge population
    pop=[pop ; popc ; popm];

    %sort merge population
    pop= sortrows(pop,nVar+1);

    %truncate population
    pop=pop(1:nPop,:);

    %Best solution
    BestSol=pop(1,:);
    BestCost=pop(1,nVar+1);
    WorstCost=pop(nPop,nVar+1);

    TempMean=mean(pop);

    if mm==0
        CostAvgInEachStep= TempMean(end);
    else
        CostAvgInEachStep=[CostAvgInEachStep ;TempMean(end)];
    end

    mm=mm+1;
    if PreferedBestCost==BestSol(end)
        mm=MaxIt+1;
    end
end
plot(CostAvgInEachStep);
end

%%%%%%%%%%%%%%%%
function [out]=iCreatePopulation(VarMin, VarMax, nVar,nPop)
out= (VarMax-VarMin).*rand(nPop,nVar) + VarMin;
end
%%%%%%%%%%%%%%%%%%
function [out]=iCrossover(in1,in2,gamma,VarMin,VarMax)
alpha=unifrnd(-gamma,1+gamma,size(in1));

out1=alpha.*in1+(1-alpha).*in2;
out2=alpha.*in2+(1-alpha).*in1;

% out1=max(out1,VarMin);
% out1=min(out1,VarMax);
%
% out2=max(out2,VarMin);
% out2=min(out2,VarMax);

out=[out1 ; out2];
end
%%%%%%%%%%%%%%%%%%
function [out]=iFittness(in)
nCols = size(in, 2);
out=0;
for i= 1 : nCols
    out = out + in(i)*in(i) - 10*cos(2*pi*in(i));
end
out=10*nCols+out;
end
%%%%%%%%%%%%%%%%%%
function [out]=iMutation(in,mu,VarMin,VarMax)
pu=2*round(mu*(size(in,2)-1)/2);
for i=0:pu
    r=randi([1 (size(in,2)-1)]);
    in(1,r)= (VarMax-VarMin).*rand(1,1) + VarMin;

end
out=in;
end
%%%%%%%%%%%%%%%%%%
function [out]=iRouleteWheelSelection(in)
r=rand;
c=cumsum(in);
out=find(r<=c,1,'first');
end
%%%%%%%% The END %%%%%%%%%%
