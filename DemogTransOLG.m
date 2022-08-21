%% Solves a deterministic OLG, including transition path: Demographic Transition
% Model has a pension system and progressive taxation.
%
% We consider a 'demographic transition'. This consists of a path for the
% conditional survival probabilities 'sj' and a path for the growth rate of
% the population of (model) age 1 'n'. Because sj is an age-dependent
% parameter the ParamPath.sj will be J-by-T. While ParamPath.n will be 1-by-T.
% Note that n is not the population growth rate, it is the growth rate of
% the age-1 population (where age 1 is understood to be model age 1, not
% actual age 1). We also use an initial distribution over age, called the
% 'age weights' mewj.
%
% Based on these we use VFI Toolkit command 'UpdateAgeWeights()' to create
% the correspondings ParamPath.mewj for the age weights.
%
% We can then just pass ParamPath to the transition command and it does the
% rest, automatically handling that some of these parameters are N_j-by-1
% (so their transition path is N_j-by-T) and realising that one of the
% ParamPath corresponds to the age weights.
%
% You obviously need to allow more points on the transition than the number of periods during 
% which your demographics change as the model needs to settle into the
% final equilibrium after the changing of demographics.
%
% We do the simplest thing and just have the pension system run a balanced
% budget every period. The pension tax rate is fixed, and the pension
% benefit amount will be adjusted to balanced the budget. Similarly the
% government has a progressive tax, and then government spending will be
% adjusted to run a balanced budget every period.
%
% One fun thing we can do now is plot a 'demographic pyramid'. In fact I do
% two, one for the initial agent distribution and one for the final agent
% distribution.
%
% This model is a deterministic OLG (no idiosyncratic shocks). You could
% easily add idiosyncratic shocks by changing z (that is changing n_z, z_grid and pi_z; 
% you would need to modify the return function and all the FnsToEvalaute to depend on 
% it where appropriate)
%
% Note that the initial guess for the general equilibrium is miles away
% from the actual solution. This makes the runtime much larger than need be
% if you used a decent inital guess. I leave it like this to demonstrate
% that you can typically solve without needing to know much about what the
% general eqm will look like.
%
% This code uses the actual demographic information for the USA from 1997
% to 2008. It is imported by 'Import_USdeathprobabilities1997to2018' and 
% 'Import_USpopulationAge20_1990to2019' for sj, and for n and initial mewj, 
% respectively. If you look inside these two scripts then it will 
% explain exactly where the original data was downloaded from.
%
% One decision variable: labour hours worked
% One endogenous state variable: assets
% No stochastic exogenous state variables
% Age


% Lets model agents from age 20 to age 100, so 81 periods
Params.agejshifter=19; % Age 20 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=100-Params.agejshifter; % =81, Number of period in life-cycle

% Grid sizes to use
n_d=31; % Endogenous labour choice (fraction of time worked)
n_a=401;
n_z=1; % This is how the VFI Toolkit thinks about deterministic models
N_j=Params.J; % Number of periods in finite horizon

transpathoptions.fastOLG=0; % If you set this to 1 then it is much faster, but requires substantially more gpu memory

figure_c=0; % I like to use a counter for the figures. Makes it easier to keep track of them when editing.

%% Parameters

% Discount rate
Params.beta = 0.99;
% Preferences
Params.sigma = 2; % Coeff of relative risk aversion (curvature of consumption)
Params.eta = 1.5; % Curvature of leisure (This will end up being 1/Frisch elasty)
Params.psi = 10; % Weight on leisure

Params.A=1; % Aggregate TFP. Not actually used anywhere.
% Production function
Params.alpha = 0.3; % Share of capital
Params.delta = 0.1; % Depreciation rate of capital

% Warm-glow of bequest
Params.warmglowparam1=1;
Params.warmglowparam2=2;

% Demographics
Params.Jr=67-Params.agejshifter; % Retirement age is 67 (remember j=1 is age 20) (You work when younger, not working from Jr on)
Params.agej=(1:1:Params.J)'; % Current 'j' age, so can check when you retire.
% Population growth rate
% E.g., Params.n=0.02; % percentage rate (expressed as fraction) at which population growths
% Note: n here is not the growth rate of the population (as not the whole
% population is anyway modelled). It is the growth rate of the population
% of age j=1 (here age 20).
Import_USpopulationAge20_1990to2019 % See this script for exact source of data
% We want 1997 to 2018 to match the sj
npath=USpopulationgrowthrate_Age20_1991to2019(7:end-1);
% I want to avoid negative population growth rates (Just because I have
% been too lazy to consider if they are possible, or will cause errors; I expect they would be fine but haven't tested. They may be a problem for the final period, but I doubt it, haven't thought it through properly though.)
npath(npath<=0)=0;
Params.n_initial=npath(1);
Params.n_final=npath(end);
% We begin by using n_initial
Params.n=Params.n_initial;

% The same code 'Import_USpopulationAge20_1990to2019' has also collected the age-distribution
Params.mewj_initial=USpopulation_Age20to100_1997/sum(USpopulation_Age20to100_1997);

% Conditional survival probabilities: sj is the probability of surviving to be age j+1, given alive at age j
% Most countries have calculations of these (as they are used by the government departments that oversee pensions)
% In fact I will here get data on the conditional death probabilities, and then survival is just 1-death.
% Here I just use them for the US, taken from "National Vital Statistics Report" "United States Life Tables" from 1997 to 2018 (https://www.cdc.gov/nchs/products/nvsr.htm note: there are spreadsheet versions via: https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/ open links and some contain a Table1.xls which are the relevant numbers)
% I took them from first column (Probability of dying between ages x and x+1; qx) of Table 1 (Total Population)
% Same documents provide the probabilities seperately for male and female, and by race (I don't use them, just mentioning)
% I copy-pasted them all into a spreadsheet and then saved it to csv, so that it can be easily imported here.
Import_USdeathprobabilities1997to2018
% Creates USdeathprobabilities1997to2018table and USdeathprobabilities1997to2018
% First, import the conditional death probabilities for 1997 to use as initial year
% Conditional death probabilities (dj covers Ages 0 to 100)
Params.dj_initial=USdeathprobabilities1997to2018table.year1997; 
Params.sj_initial=1-Params.dj_initial(21:101); % Conditional survival probabilities
Params.sj_initial(end)=0; % This is actually unnecessary, but looks nicer if we graph them
% Second, import the conditional death probabilities for 2018 to use as final year
Params.dj_final=USdeathprobabilities1997to2018table.year2018; 
Params.sj_final=1-Params.dj_final(21:101); % Conditional survival probabilities
Params.sj_final(end)=0; % This is actually unnecessary, but looks nicer if we graph them
% We begin by using sj_initial
Params.sj=Params.sj_initial;

% Labor efficiency units depend on age
Params.kappa_j=[linspace(1,3,50-Params.agejshifter), linspace(3,2,(Params.Jr-1)-(50-Params.agejshifter)),zeros(1,Params.J-Params.Jr+1)];
    % These are not done seriously, really they should be set to something like 'average hourly wage conditional on age' in the data.
    % I have made them increase until age 50 (j=31), then decrease, and then be zero from retirement at age 67.
    % For parameters that depend on age, we just make them a vector with a
    % length the same as the number of periods, VFI Toolkit then handles
    % them automatically.

% Taxes
Params.tau=0.15;
% In addition to payroll tax rate tau, which funds the pension system we will add a progressive 
% income tax which funds government spending.
% The progressive income tax takes the functional form:
% IncomeTax=eta1+eta2*log(Income)*Income; % This functional form is empirically a decent fit for the US tax system
% And is determined by the two parameters
Params.eta1=0.09; % eta1 will be determined in equilibrium to balance gov budget constraint
Params.eta2=0.053;

% Government spending
Params.GdivYtarget = 0.15; % Government spending as a fraction of GDP (this is essentially just used as a target to define a general equilibrium condition)

% Government debt target
Params.BdivYtarget=0.2;

% Fiscal response function: G=G_final+phi1*(taxrevenue-total gov spending)+phi2*(BdivYtarget-DdivY)
Params.phi1=0.6;
Params.phi2=0.5;
% Note: G is included in total gov spending, so can do rearrangement to get
%  G=(1/(1+phi1))*(G_final+phi1*(taxrevenue-total gov spending excluding G)+phi2*(BdivYtarget-DdivY))

%% Some initial values/guesses for variables that will be determined in general eqm
Params.r=0.11; % Interest rate
Params.pension=0.6; % Initial guess (this will be determined in general eqm)
Params.AccidentBeq=0.04; % Accidental bequests (this is the lump sum transfer)
Params.G=0.16; % Government expenditure
Params.eta1=0.18;

%% Grids
Params.kmax=5;
a_grid=linspace(0,Params.kmax,n_a)'; % Note: Implicitly, we are imposing borrowing constraint a'>=0
z_grid=1;
pi_z=1;

% Grid for labour choice
h_grid=linspace(0,1,n_d)';
% Switch into toolkit notation
d_grid=h_grid;

%% Now, create the return function 
DiscountFactorParamNames={'beta','sj'};

ReturnFn=@(h,aprime,a,z,agej,r,A,delta,alpha,sigma,psi,eta,Jr,pension,tau,kappa_j,J,warmglowparam1,warmglowparam2,AccidentBeq, eta1,eta2) DemogTransOLG_ReturnFn(h,aprime,a,z,agej,r,A,delta,alpha,sigma,psi,eta,Jr,pension,tau,kappa_j,J,warmglowparam1,warmglowparam2,AccidentBeq, eta1,eta2)
ReturnFnParamNames={'agej','r','A','delta','alpha','sigma','psi','eta','Jr','pension','tau','kappa_j','J','warmglowparam1','warmglowparam2','AccidentBeq','eta1','eta2'}; %It is important that these are in same order as they appear in 'FiscalOLG_ReturnFn'

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium

disp('Test ValueFnIter')
vfoptions=struct(); % Just using the defaults.
tic;
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
toc

%% Initial distribution of agents at birth (j=1)
jequaloneDist=zeros(n_a,n_z); % n_a by n_z.
jequaloneDist(1,1)=1; % Everyone is born with zero assets

%% Agents age distribution
% Normally you might calculate mewj using the code in the following lines
% that is commented out. But in the present example the 1997 value for
% parameter n is unusually high (relative to average of n) and so you would
% get somewhat 'extreme' result. We do what anyway seems the better option
% and just imported it from the data.
Params.mewj=Params.mewj_initial;

% % Many OLG models include some kind of population growth, and perhaps
% % some other things that create a weighting of different ages that needs to
% % be used to calculate the stationary distribution and aggregate variable.
% Params.mewj=ones(1,Params.J); % Marginal distribution of households over age
% for jj=2:length(Params.mewj)
%     Params.mewj(jj)=Params.sj(jj-1)*Params.mewj(jj-1)/(1+Params.n);
% end
% Params.mewj=Params.mewj./sum(Params.mewj); % Normalize to one

AgeWeightsParamNames={'mewj'}; % Many finite horizon models apply different weights to different 'ages'; eg., due to survival rates or population growth rates.

%% Test
disp('Test StationaryDist')
simoptions=struct(); % Just use the defaults
tic;
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
toc

%% General eqm variables
GEPriceParamNames={'r','pension','AccidentBeq','G','eta1'};

%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)

% Stationary Distribution Aggregates (important that ordering of Names and Functions is the same)
FnsToEvaluate.H = @(h,aprime,a,z) h; % Aggregate labour supply (in efficiency units)
FnsToEvaluate.L = @(h,aprime,a,z,kappa_j) kappa_j*h; % Aggregate  physical capital
FnsToEvaluate.K = @(h,aprime,a,z) a; % Aggregate labour supply (in efficiency units)
FnsToEvaluate.PensionSpending = @(h,aprime,a,z,pension,agej,Jr) (agej>=Jr)*pension; % Total spending on pensions
FnsToEvaluate.AccidentalBeqLeft = @(h,aprime,a,z,sj) aprime*(1-sj); % Accidental bequests left by people who die
FnsToEvaluate.IncomeTaxRevenue = @(h,aprime,a,z,eta1,eta2,kappa_j,r,delta,alpha,A) DemogTransOLG_ProgressiveIncomeTaxFn(h,aprime,a,z,eta1,eta2,kappa_j,r,delta,alpha,A); % Revenue raised by the progressive income tax (needed own function to avoid log(0) causing problems)

% General Equilibrium conditions (these should evaluate to zero in general equilbrium)
GeneralEqmEqns.capitalmarket = @(r,K,L,alpha,delta,A) r-alpha*A*(K^(alpha-1))*(L^(1-alpha)); % interest rate equals marginal product of capital net of depreciation
GeneralEqmEqns.pensions = @(PensionSpending,tau,L,r,A,alpha,delta) PensionSpending-tau*(A*(1-alpha)*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*L; % Retirement benefits equal Payroll tax revenue: pension*fractionretired-tau*w*H
GeneralEqmEqns.bequests = @(AccidentalBeqLeft,AccidentBeq,n) AccidentalBeqLeft/(1+n)-AccidentBeq; % Accidental bequests received equal accidental bequests left
GeneralEqmEqns.Gtarget = @(G,GdivYtarget,A,K,L,alpha) G-GdivYtarget*(A*K^(alpha)*(L^(1-alpha))); % G is equal to the target, GdivYtarget*Y
GeneralEqmEqns.govbudget = @(G,IncomeTaxRevenue) G-IncomeTaxRevenue; % Government budget balances (note that pensions are a seperate budget)
% Notice that the government budget now includes the interest on government debt (in stationary equilibrium debt will be constant)

%% Test
disp('Test AggVars')
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid);

%% Solve for the General Equilibrium

heteroagentoptions.verbose=1;
[p_eqm_initial,~,GeneralEqmConditions_initial]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, [], [], GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
% p_eqm_initial contains the general equilibrium parameter values
% Put this into Params so we can calculate things about the initial equilibrium
Params.r=p_eqm_initial.r;
Params.pension=p_eqm_initial.pension;
Params.AccidentBeq=p_eqm_initial.AccidentBeq;
Params.G=p_eqm_initial.G;
Params.eta1=p_eqm_initial.eta1;

% Calculate a few things related to the initial equilibrium.
[V_initial, Policy_initial]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
StationaryDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy_initial,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
% Can just use the same FnsToEvaluate as before.
AgeConditionalStats_initial=LifeCycleProfiles_FHorz_Case1(StationaryDist_initial,Policy_initial,FnsToEvaluate,[],Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid);
% Just to compare initial and final we will also do
AggVars_initial=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_initial, Policy_initial, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid);


% Store for use later
Params_initial=Params;

%%
% save ./SavedOutput/DemogTransOLG_initial.mat p_eqm_initial GeneralEqmConditions_initial StationaryDist_initial Params_initial AggVars_initial


%% Before we can calculate the final equilibrium we need the 'age weights' that will be used.
% For this reason we first set up the 'demographic transition' that we will
% later use when calculating the transition path. We can then just take the
% final period distribution over age and use this for the final equilbrium.
T=100 % Number of time periods
% First time I used T=150 (just a guess, deliberately large), but appeared clear that everything was settled
% easily within first 100 periods, so now using this.

%% Set up the demographic transition to an 'aging socitiey'
% Paths on n (growth rate of the population of model-age 1), sj (the
% conditional survival probability).
% From these we first use UpdateAgeWeights() to create the path on the
% distribution over ages, mewj.
% We can then simply include these three as the ParamPath for which we want
% to calculate the transition.

% Recall, the death probabilities from 1997 to 2018 were imported as a
% table but also as the following matrix, so we can just use
ParamPath.sj=1-USdeathprobabilities1997to2018(Params.agejshifter+2:Params.J+Params.agejshifter+1,:); % Matrix starts from age zero, so had to have an extra one to each of first and last ages
% Note that this just covers the first periods, to get the model to settle
% to the final stationary eqm we want to repeat the final column.
ParamPath.sj=[ParamPath.sj,ParamPath.sj(:,end).*ones(Params.J,T-size(ParamPath.sj,2))];

% Note: n here is not the growth rate of the population (as not the whole
% population is anyway modelled). It is the growth rate of the population
% of age j=1 (here age 20).
Import_USpopulationAge20_1990to2019
% We want 1997 to 2018 to match the sj
npath=USpopulationgrowthrate_Age20_1991to2019(7:end-1);
ParamPath.n=[npath,npath(end)*ones(1,T-length(npath))]; % Note: the 2018 growth rate of population of age 20 was essentially 0


% To handle aging we also need to update the 'marginal distribution of the population over age', 
% which we have denoted mewj. We just set this up beforehand and then
% input it as a ParamPath.mewj (because it is all exogenous we can
% calculate it beforehand). VFI Toolkit has a function to make this easy:
ParamPath.mewj=UpdateAgeWeights(Params.mewj_initial,ParamPath.sj,ParamPath.n);
% inputs: initial age distribution (1-by_N_j), the path on age conditional survival
% probability (T-by-N_j), the path on the the growth rate of the (model period)
% age 1 population (T-by-1).

Params.mewj_final=ParamPath.mewj(:,end); % Just to make it easy to call later

%% Solve for the final equilibrium
Params.sj=Params.sj_final;
Params.n=Params.n_final;
% Note that changing sj necessitates changing mewj as well
Params.mewj=Params.mewj_final;

% I need a better initial guess for bequests as was otherwise having difficulty converging
Params.AccidentBeq=0.05;

% We need to fix the taxation (we want to keep it at initial value)
% So drop eta1 from the general eqm prices
GEPriceParamNames={'r','pension','AccidentBeq','G'};
% As a result, we cannot target the size of government, so drop the target for G as fraction of Y
GeneralEqmEqns=rmfield(GeneralEqmEqns,'Gtarget');
% Note: The alternative to this would be to keep eta1 and Gtarget, and then in the transition to add ParamsPath.eta1.
% Here we are doing that G adjusts to keep eta1 unchanged (the tax rate is unchanged, but the tax revenue will change). 
% This alternative would instead have eta1 adjust to keep G at GdivYtarget (in final period, not whole transition path).
% Which of this version and the alternative is preferred depends on exactly what kind of reform you want to model.

[p_eqm_final,~,GeneralEqmConditions_final]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, [], [], GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
% p_eqm_final contains the general equilibrium parameter values
% Put this into Params so we can calculate things about the initial equilibrium
Params.r=p_eqm_final.r;
Params.pension=p_eqm_final.pension;
Params.AccidentBeq=p_eqm_final.AccidentBeq;
Params.G=p_eqm_final.G;
% Params.eta1=p_eqm_final.eta1;

% Calculate a few things related to the initial equilibrium.
[V_final, Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
StationaryDist_final=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy_final,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
% Can just use the same FnsToEvaluate as before.
AgeConditionalStats_final=LifeCycleProfiles_FHorz_Case1(StationaryDist_final,Policy_final,FnsToEvaluate,[],Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid);
% Just to compare initial and final we will also do
AggVars_final=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_final, Policy_final, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid);


Params_final=Params;

%%
% save ./SavedOutput/DemogTransOLG_final.mat p_eqm_final GeneralEqmConditions_final  V_final Params_final AggVars_final

%% This is just so that it is easy to restart from the transition path itself
% save ./SavedOutput/DemogTransOLG_more.mat
% load ./SavedOutput/DemogTransOLG_more.mat 

%% Plot the life cycle profiles of capital and labour for the inital and final eqm.

figure_c=figure_c+1;
figure(figure_c)
subplot(3,1,1); plot(1:1:Params.J,AgeConditionalStats_initial.H.Mean,1:1:Params.J,AgeConditionalStats_final.H.Mean)
title('Life Cycle Profile: Hours Worked')
legend({'Initial','Final'})
subplot(3,1,2); plot(1:1:Params.J,AgeConditionalStats_initial.L.Mean,1:1:Params.J,AgeConditionalStats_final.L.Mean)
title('Life Cycle Profile: Labour Supply')
legend({'Initial','Final'})
subplot(3,1,3); plot(1:1:Params.J,AgeConditionalStats_initial.K.Mean,1:1:Params.J,AgeConditionalStats_final.K.Mean)
title('Life Cycle Profile: Assets')
legend({'Initial','Final'})
% saveas(figure_c,'./SavedOutput/Graphs/DemogTransOLG_LifeCycleProfiles','pdf')

%% Compute the Transition
% T=100 % Number of time periods
% Already declare T above as needed to create path on mewj, the population age-weights.

% We consider a 'demographic transition'. This consists of a path for the
% conditional survival probabilities 'sj' and a path for the growth rate of
% the population of (model) age 1 'n'. Because sj is an age-dependent
% parameter the ParamPath.sj will be J-by-T. The path for n will be 1-by-T.
%
% We then create the implied ParamPath.mewj of the age-distribution of
% agents using UpdateAgeWeights().
%
% The transition path commands automatically recognise that some of our
% ParamPath are J-by-T and interpret these as transitions in age-dependent
% parameters.

%% Setting up the things that will be determined in the transition paths
% GEPriceParamNames={'r','pension','AccidentBeq','G'};
% You do not need to tell the transition path what the general eqm prices
% are, it knows that they are the things in PricePath0.

% We need to give an initial guess for the price path on the Gen Eqm variables:
PricePath0.r=[linspace(p_eqm_initial.r, p_eqm_final.r, floor(T/3))'; p_eqm_final.r*ones(T-floor(T/3),1)]; % For each price/parameter that will be determined in general eqm over transition path PricePath0 is matrix of size T-by-1
PricePath0.pension=[linspace(p_eqm_initial.pension, p_eqm_final.pension, floor(T/3))'; p_eqm_final.pension*ones(T-floor(T/3),1)]; % For each price/parameter that will be determined in general eqm over transition path PricePath0 is matrix of size T-by-1
PricePath0.AccidentBeq=[linspace(p_eqm_initial.AccidentBeq, p_eqm_final.AccidentBeq, floor(T/3))'; p_eqm_final.AccidentBeq*ones(T-floor(T/3),1)]; % For each price/parameter that will be determined in general eqm over transition path PricePath0 is matrix of size T-by-1
PricePath0.G=[linspace(p_eqm_initial.G, p_eqm_final.G, floor(T/3))'; p_eqm_final.G*ones(T-floor(T/3),1)]; % For each price/parameter that will be determined in general eqm over transition path PricePath0 is matrix of size T-by-1

TransPathGeneralEqmEqns.capitalmarket = @(r,K,L,alpha,delta) r-alpha*(K^(alpha-1))*(L^(1-alpha)); % interest rate equals marginal product of capital net of depreciation
TransPathGeneralEqmEqns.pensionbalance = @(PensionSpending,tau,L,r,A,alpha,delta) PensionSpending-tau*(A*(1-alpha)*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*L; % Retirement benefits equal Payroll tax revenue: pension*fractionretired-tau*w*H
TransPathGeneralEqmEqns.bequests = @(AccidentalBeqLeft_tminus1,AccidentBeq,n) AccidentalBeqLeft_tminus1/(1+n)-AccidentBeq; % Accidental bequests received equal accidental bequests left
TransPathGeneralEqmEqns.govbudgetbalance = @(G,IncomeTaxRevenue) G-IncomeTaxRevenue; % Government budget balances (note that pensions are a seperate budget), has a -0.1(B-BdivYtarget*Y) term that represents the government adjusting government spending to get debt toward target
% Note: require the pension system to run a balanced budget every period.
% Require the government to run a balanced budget every period.

% To be able to use a '_tminus1' as input, we have to provide the initial value
transpathoptions.initialvalues.AccidentalBeqLeft=Params_initial.AccidentBeq/(1+Params_initial.n);

transpathoptions.GEnewprice=3;
% Need to explain to transpathoptions how to use the GeneralEqmEqns to
% update the general eqm transition prices (in PricePath).
transpathoptions.GEnewprice3.howtoupdate=... % a row is: GEcondn, price, add, factor
    {'capitalmarket','r',0,0.5;... % capitalmarket GE condition will be positive if r is too big, so subtract
    'pensionbalance','pension',0,0.5;... % pensionbalance GE condition will be positive if pension is too big, so subtract
    'bequests','AccidentBeq',1,0.5;... % bequests GE condition will be positive if AccidentBeq is too small, so add
    'govbudgetbalance','G',0,0.5}; % govbudget GE condition will be positive if G is too big, so subtract
% Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
% Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0


% Now just run the TransitionPath_Case1 command (all of the other inputs are things we had already had to 
% define to be able to solve for the initial and final equilibria)
transpathoptions.verbose=1;
transpathoptions.graphpricepath=1;
transpathoptions.graphaggvarspath=1;
transpathoptions.tolerance=5*10^(4); % default would be 10^(-4)
PricePath=TransitionPath_Case1_FHorz(PricePath0, ParamPath, T, V_final, StationaryDist_initial, n_d, n_a, n_z, N_j, pi_z, d_grid,a_grid,z_grid, ReturnFn,  FnsToEvaluate, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightsParamNames, transpathoptions);
% The path for the stock variables is simply included as part of PricePath. While the difference between a stock variable and a price variable
% is important when computing the general equilibrium transition path itself it ceases to be important from then on (similar to how the
% distinction between a stock variable and a general eqm price was unimportant when we computed the initial and final stationary equilibria).

%%
% save ./SavedOutput/DemogTransOLG_TransPath.mat ParamPath PricePath T V_final StationaryDist_initial

%%
% You can calculate the value and policy functions for the transition path
[VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions);

% You can then use these to calculate the agent distribution for the transition path
AgentDistPath=AgentDistOnTransPath_Case1_FHorz(StationaryDist_initial,PricePath, ParamPath, PolicyPath, AgeWeightsParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params);

%% Graph some things from the transition

FnsToEvaluate.K = @(h, aprime,a,z) a; % Aggregate physical capital
FnsToEvaluate.H = @(h, aprime,a,z) h; % Aggregate hours worked
FnsToEvaluate.L = @(h, aprime,a,z,kappa_j) kappa_j*h; % Aggregate labour supply (in efficiency units)
FnsToEvaluate.C = @(h, aprime,a,z,agej,Jr,r,pension,tau,kappa_j,alpha,delta,A,eta1,eta2) DemogTransOLG_ConsumptionFn(h,aprime,a,z,agej,Jr,r,pension,tau,kappa_j,alpha,delta,A,eta1,eta2);
FnsToEvaluate.I = @(h, aprime,a,z, delta) aprime-(1-delta)*a; % Note that this will be zero by definition here as capital is exogenous.
FnsToEvaluate.IncomeTaxRevenue = @(h,aprime,a,z,eta1,eta2,kappa_j,r,delta,alpha,A) DemogTransOLG_ProgressiveIncomeTaxFn(h,aprime,a,z,eta1,eta2,kappa_j,r,delta,alpha,A); % Revenue raised by the progressive income tax (needed own function to avoid log(0) causing problems)

AggVars_initial=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_initial, Policy_initial, FnsToEvaluate, Params_initial, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid);
AggVarsPath=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate,AgentDistPath,PolicyPath, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, transpathoptions);

Y_initial=Params.A*(AggVars_initial.K.Mean.^Params.alpha).*(AggVars_initial.L.Mean.^(1-Params.alpha));
KdivL_initial=(Params_initial.r/(Params.alpha*Params.A)).^(1/(Params.alpha-1));
w_initial=Params.A*(1-Params.alpha)*(KdivL_initial.^Params.alpha);

Path_Y=Params.A*(AggVarsPath.K.Mean.^Params.alpha).*(AggVarsPath.L.Mean.^(1-Params.alpha));
Path_KdivL=(PricePath.r/(Params.alpha*Params.A)).^(1/(Params.alpha-1));
Path_w=Params.A*(1-Params.alpha)*(Path_KdivL.^Params.alpha);

figure_c=figure_c+1;
figure(figure_c)
subplot(4,2,1); plot([Y_initial;Path_Y])
axis tight
title('Output (Y)')
subplot(4,2,2); plot([AggVars_initial.C.Mean;AggVarsPath.C.Mean])
axis tight
title('Consumption (C)')
subplot(4,2,3); plot([AggVars_initial.K.Mean;AggVarsPath.K.Mean])
axis tight
title('Physical Capital (K)')
subplot(4,2,4); plot([AggVars_initial.L.Mean;AggVarsPath.L.Mean])
axis tight
title('Labor Supply (L)')
subplot(4,2,5); plot([w_initial;Path_w])
axis tight
title('Wages (w)')
subplot(4,2,6); plot([Params_initial.pension;PricePath.pension])
axis tight
title('Pensions')
subplot(4,2,7); plot([AggVars_initial.I.Mean;AggVarsPath.I.Mean])
axis tight
title('Investment (I)')
subplot(4,2,8); plot([Params_initial.r;PricePath.r] - Params.delta)
axis tight
title('Interest Rate (r)')
% saveas(figure_c,'./SavedOutput/Graphs/DemogTransOLG_TransPathFig1','pdf')

figure_c=figure_c+1;
figure(figure_c)
subplot(2,1,1); plot([Params_initial.G;PricePath.G])
axis tight
title('Government Consumption (G)')
subplot(2,1,2); plot([AggVars_initial.IncomeTaxRevenue.Mean;AggVarsPath.IncomeTaxRevenue.Mean])
axis tight
title('Income Tax Revenue')
% saveas(figure_c,'./SavedOutput/Graphs/DemogTransOLG_TransPathFig2','pdf')

%% Plot the initial and final 'demographic pyramids'
InitialAgeDist=shiftdim(sum(sum(AgentDistPath(:,:,:,1),1),2),2); % Just the marginal distribution over ages (for "0-th" period of transition)
FinalAgeDist=shiftdim(sum(sum(AgentDistPath(:,:,:,end),1),2),2); % Just the marginal distribution over ages (for final period of transition)

figure_c=figure_c+1;
figure(figure_c)
hold on
pyramidright = barh((1:1:Params.J)+Params.agejshifter,InitialAgeDist/2,'hist'); % I need the divided by two to get the symmetry (half on each 'side')
pyramidleft = barh((1:1:Params.J)+Params.agejshifter,-InitialAgeDist/2,'hist'); % 'minus', so it goes to other 'side'
set(pyramidright,'FaceColor','b')
set(pyramidleft,'FaceColor','b')
hold off
ylabel('Age')
xlabel('Fraction of population')
title('Initial Demographic Pyramid')
% saveas(figure_c,'./SavedOutput/Graphs/DemogTransOLG_DemographicPyramid1','pdf')

figure_c=figure_c+1;
figure(figure_c)
hold on
pyramidright = barh((1:1:Params.J)+Params.agejshifter,FinalAgeDist/2,'hist'); % I need the divided by two to get the symmetry (half on each 'side')
pyramidleft = barh((1:1:Params.J)+Params.agejshifter,-FinalAgeDist/2,'hist'); % 'minus', so it goes to other 'side'
set(pyramidright,'FaceColor','b')
set(pyramidleft,'FaceColor','b')
hold off
ylabel('Age')
xlabel('Fraction of population')
title('Final Demographic Pyramid')
% saveas(figure_c,'./SavedOutput/Graphs/DemogTransOLG_DemographicPyramid2','pdf')


%% Graph the agent distribution over ages to check that the demographic transition is working correctly
% Obviously you don't normally need to do this. Is simply to illustrate what is going on.
% VFI Toolkit works by taking the 'age weights' mewj and imposing them on the agent distribution.

% To check that the demographic transition is working, let's plot some of
% the ParamPath.mewj against the age-distribution of the AgentDistPath.
figure_c=figure_c+1;
figure(figure_c)
subplot(2,1,1); plot(1:1:Params.J,ParamPath.mewj(:,1))
hold on
subplot(2,1,1); plot(1:1:Params.J,ParamPath.mewj(:,6))
subplot(2,1,1); plot(1:1:Params.J,ParamPath.mewj(:,11))
subplot(2,1,1); plot(1:1:Params.J,ParamPath.mewj(:,16))
subplot(2,1,1); plot(1:1:Params.J,ParamPath.mewj(:,21))
hold off
title('Age Distribution over Transition: mewj')
AgentDistPath_age=shiftdim(sum(sum(AgentDistPath,1),2),2); % Just the distribution over ages (and transition periods)
subplot(2,1,2); plot(1:1:Params.J,AgentDistPath_age(:,1))
hold on
subplot(2,1,2); plot(1:1:Params.J,AgentDistPath_age(:,6))
subplot(2,1,2); plot(1:1:Params.J,AgentDistPath_age(:,11))
subplot(2,1,2); plot(1:1:Params.J,AgentDistPath_age(:,16))
subplot(2,1,2); plot(1:1:Params.J,AgentDistPath_age(:,21))
hold off
legend('0th period','5th period','10th period','15th period','20th period')
title('Age Distribution over Transition: AgentDist')
% saveas(figure_c,'./SavedOutput/Graphs/DemogTransOLG_TransPathFig3','pdf')
% As is clear from the graphs, the Agent Distribution now contains the age weights.



