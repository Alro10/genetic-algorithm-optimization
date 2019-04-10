%% Continuous Genetic Algorithm
%
% minimizes the objective function designated in ff
%
% Before beginning, set all the parameters in parts I, II, and III
%% Setup the GA
clear all; clc;
ff='testfunction'; % objective function
npar=2; % number of optimization variables
varhi=2; varlo=-1; % variable limits
%% II Stopping criteria
maxit=100; % max number of iterations
mincost=-999999; % minimum cost
%% III GA parameters
popsize=100; % set population size
mutrate=.01; % set mutation rate
selection=0.8; % fraction of population kept
Nt=npar; % continuous parameter GA Nt=#variables
keep=floor(selection*popsize); % #population members that survive
nmut=ceil((popsize-1)*Nt*mutrate); % total number of mutations
M=ceil((popsize-keep)/2); % number of matings // CEIL   Round towards plus infinity.
%% Create the initial population
iga=0; % generation counter initialized
par=(varhi-varlo)*rand(popsize,npar)+varlo; % random
Coords{1}=par;
cost=feval(ff,par); % calculates population cost using ff
[cost,ind]=sort(cost,'descend'); % min cost in element 1// SORT in ascending or descending order.
par=par(ind,:); % sort continuous
minc(1)=max(cost); % minc contains max of
meanc(1)=mean(cost); % meanc contains mean of population


%% Iterate through generations (Main Loop)
while iga<maxit
iga=iga+1; % increments generation counter
%_______________________________________________________
% Pair and mate
M=ceil((popsize-keep)/2); % number of matings
prob=flipud([1:keep]'/sum([1:keep])); % weights chromosomes
odds=[0 cumsum(prob(1:keep))']; % probability distribution function
pick1=rand(1,M); % mate #1 (vector of length M with random #s between 0 and 1)
pick2=rand(1,M); % mate #2
% ma and pa contain the indices of the chromosomes that will mate
% Choosing integer k with probability p(k)
%
ic=1;
while ic<=M
for id=2:keep+1
if pick1(ic)<=odds(id) && pick1(ic)>odds(id-1)
ma(ic)=id-1;
end
if pick2(ic)<=odds(id) && pick2(ic)>odds(id-1)
pa(ic)=id-1;
end
end
ic=ic+1;
end

%_______________________________________________________
% Performs mating using single point crossover
ix=1:2:keep; % index of mate #1
xp=ceil(rand(1,M)*Nt); % crossover point
r=rand(1,M); % mixing parameter
for ic=1:M
xy=par(ma(ic),xp(ic))-par(pa(ic),xp(ic)); % ma and pa mate
par(keep+ix(ic),:)=par(ma(ic),:); % 1st offspring
par(keep+ix(ic)+1,:)=par(pa(ic),:); % 2nd offspring
par(keep+ix(ic),xp(ic))=par(ma(ic),xp(ic))-r(ic).*xy; % 1st
par(keep+ix(ic)+1,xp(ic))=par(pa(ic),xp(ic))+r(ic).*xy; % 2nd
if xp(ic)<npar % crossover when last variable not selected
par(keep+ix(ic),:)=[par(keep+ix(ic),1:xp(ic))
par(keep+ix(ic)+1,xp(ic)+1:npar)];
par(keep+ix(ic)+1,:)=[par(keep+ix(ic)+1,1:xp(ic))
par(keep+ix(ic),xp(ic)+1:npar)];
end % if
end

%_______________________________________________________
% Mutate the population
mrow=sort(ceil(rand(1,nmut)*(popsize-1))+1);
mcol=ceil(rand(1,nmut)*Nt);
for ii=1:nmut
par(mrow(ii),mcol(ii))=(varhi-varlo)*rand+varlo;
% mutation
end % ii
%_______________________________________________________
% The new offspring and mutated chromosomes are
% evaluated
cost=feval(ff,par);
%_______________________________________________________
% Sort the costs and associated parameters
[cost,ind]=sort(cost,'descend');
par=par(ind,:);
Coords{iga+1}=par;
%_______________________________________________________
% Plot function  26_11_16
figure (2)
  [X,Y] = meshgrid(-1:.02:2, -1:.02:2);
        Z =sin(4*pi*X).*X-sin(4*pi*Y+pi).*Y+1;
        %hold on;
        % sin(4*pi*xx).*xx-sin(4*pi*yy+pi).*yy+1


        surf(X,Y,Z)
        xlabel('x')
        ylabel('y')
        zlabel('f(x,y)')

       hold on;
  pcolor(X,Y,Z);              % is really a SURF with its view set to directly above
 shading interp
%_______________________________________________________
% Plot offspring population at each generation
  plot3(Coords{iga+1}(:,1)',Coords{iga+1}(:,2)',cost, 'b*');
  hold off;
  %plot3(par(:,1)',par(:,2)',cost, 'b*');
  % Coords{iga+1}(:,1)'.*sin(4*pi*Coords{iga+1}(:,1)') - Coords{iga+1}(:,2)'.*sin(4*Coords{iga+1}(:,2)'+pi)+1
    %axis([-1 2 -1 2]);
    pause(0.1)
%_______________________________________________________
% Do statistics for a single nonaveraging run
minc(iga+1)=max(cost);
meanc(iga+1)=mean(cost);
%_______________________________________________________
% Stopping criteria
if iga>maxit || cost(1)<mincost
break
end
[iga cost(1)];
end %iga
%% Displays the output
day=clock;
disp(datestr(datenum(day(1),day(2),day(3),day(4),day(5),day(6)),0))
disp(['optimized function is ' ff])
format short g
disp(['popsize=' num2str(popsize) ' mutrate=' num2str(mutrate) ' # par=' num2str(npar)])
disp(['#generations=' num2str(iga) ' best cost=' num2str(cost(1))])
disp('best solution')
disp(num2str(par(1,:)))
disp('continuous genetic algorithm')
figure(1)
iters=0:length(minc)-1;

plot(iters,minc,iters,meanc,'r');
xlabel('generation');ylabel('fitness');
title('Fitness function')
legend('Best individual','Mean of population','Location','east')
