% SENSPRCCREG - Example file by Jaimit Parikh to run the PRCC and Linear
% Code Edited and generalized by Patrick Hanafin
% Regression based sensitivity analysis on the QSP AD model


% 1. Get parameters of the AD model
rng(27);
parameters = antiPCSK9modelPars();
P = struct2cell(parameters);


% 2. Sample using Latin hypercube sampling to generate parameter values 
% for running simulations
parsName = string(fieldnames(parameters))'; % extract parameter names into a string
nParameters = length(parsName); % number of parameters extracted from length of the parameter names
problem.lowerBound = (horzcat(cell2mat(struct2cell(parameters)))*0.9)'; % Set lowerbound for LHS
problem.upperBound = (horzcat(cell2mat(struct2cell(parameters)))*1.1)'; % Set upperbound for LHS
problem.num_vars = nParameters; % holds number of parameters in problem structure
problem.names = parsName; % Set names of parameters in structure
N = 100*nParameters; % number of samples = 100 * Number of Parameters
samples = latinHypercubeSampling(problem,  N); % Sample parameters using LHS

% Simulating AD model for the desired parameter set and extracting state
% variables
%de = adDrugEffects(); % Extract drug effect for ODEs
IC = getInitialConditions(); % Extract initial conditions of ODEs
%EASI0=72*(2*IC(2)+2*(1-IC(1)))/4; % Baseline EASI score. To be used for aditional model output
t0 = 0; % Initial Time, set to 0
tfinal = 50; % Time used for PRCC analysis 24 weeks is used for EASI Score in AD Assessment

% Create empty arrays for user-chosen outputs of interest
% Outputs will be used for PRCC and plotting
% MUST be defined by the USER
LDLc = zeros(1, N); aPCSK9 =zeros(1,N); CP = zeros(1, N);

% Simulate model for each parameter set, N, and for each user-chosen output
pN = cell(1,N);
for ii = 1:N
    pN{ii} = updatePars(parameters, parsName, samples(ii, :));
    [T,Y] = ode15s(@(t, y)odefun(t, y, pN{ii}),...
        linspace(t0, tfinal, 50), IC);
    LDLc(ii) = Y(end, 2); %Assign ODE output of interest to USER-generated empty array
    aPCSK9(ii) = Y(end, 5); %Assign ODE output of interest to USER-generated empty array
    CP(ii) = Y(end, 3); %EASI Score Improvement
    %PIE24 is an output of interest and is assigned to USER-generated empty array
end

% Generate names vector for Output and Input Selected Outputs into cell
% USER MUST insert names for outputs of interest and input them into cell
% named "OUTPUTS"
OUTPUTNames=["LDLc", "aPCSK9", "CP"];
OUTPUTS={LDLc,aPCSK9,CP};

% Loop through model parameters to produce a scatter plot of each output
% Based on LHS samples
% Currently commented out to increase script run speed
%for ii =1:nParameters
    %plotScatter(samples, OUTPUTS,OUTPUTNames, parsName,ii, strcat(parsName(ii),'_scatter.png'));
%end

% 3. Perform PRCC based sensitivity analysis
% Generate PRCC results for each output of interest and enter into struct
S_OUTPUTS=struct([]);
for ii = 1:length(OUTPUTS)
    S_OUTPUTS= [S_OUTPUTS,prccAnalysis(samples, cell2mat(OUTPUTS(ii))', "partialcorr")];
end

% Generate color grid/heat map of PRCC results for all outputs of interest
% and all parameters
plotSAresults(S_OUTPUTS,OUTPUTNames,parsName, './prcc.png', 'PRCC Analysis')

% Gernate bar plots for PRCC  for each Output of Interest
for ii =1:length(OUTPUTS)
    plotSABarResults(S_OUTPUTS(ii).prcc(1:end-1, end),parsName,OUTPUTNames(ii),'_prccBar.png')
end
    
% Generate ordered Bar plots for PRCC  for each Output of Interest
for ii =1:length(OUTPUTS)
    plotSABarOderedResults(S_OUTPUTS(ii).prcc(1:end-1, end),parsName,OUTPUTNames(ii),'_prccBar.png')
end

% Scatter plot function for each output of interest based on Latin
% Hypercube Sampling of a designated parameter. Dimensions of figure based
% on total Outputs of interest
function plotScatter(samples,OUTPUTS,OUTPUTNames, parsName,parsNum, fname)
f = figure('DefaultAxesFontSize', 14, 'Position', [40, 40, 1400, 800]);
DIM=ceil(sqrt(length(OUTPUTS))); % Generate Dimension of figures based on OUTPUTS length
for ii = 1:length(OUTPUTS) % Create plot for each Output of interest
    subplot(DIM, DIM, ii);
    plot(samples(:,parsNum), cell2mat(OUTPUTS(ii)), 'bo'); xlabel(parsName(parsNum)); ylabel(OUTPUTNames(ii));
end
sgtitle(parsName(parsNum)); 
exportgraphics(f, fname, 'resolution', 300);
end

% Update Model Parameters based on Latin Hypercube Sampling
function p = updatePars(p, parsName, parsValue)
 
for ii = 1:length(parsName)
    p.(parsName(ii)) = parsValue(ii);
end

end

% Function for PRCC Bar plots for an outputs of interest
function plotSABarResults(SAdat, parsName,tname,fname) 
f = figure('DefaultAxesFontSize', 14);
bar(abs(SAdat)); xticks(1:length(parsName));xticklabels(parsName); xtickangle(45); ylabel('|PRCC|');
title(tname);
exportgraphics(f,strcat(tname,fname), 'resolution', 300);
end

% Function for PRCC Bar plots for an ouputs of interest in descending order of
% correlation
function plotSABarOderedResults(SAdat, parsName,tname,fname) 
f = figure('DefaultAxesFontSize', 14);
SAdat=abs(SAdat);
[~,I]=sort(SAdat,'descend');
bar(SAdat(I)); xticks(1:length(parsName));xticklabels(parsName(I)); ylabel('|PRCC|');
title(tname);
exportgraphics(f,strcat(tname,fname), 'resolution', 300);
end

% Function for PRCC heatmap with all outputs/parameters of interest
function plotSAresults(S_OUTPUTS,OUTPUTNames ,parsName, fname, tname)
AreaDat=zeros(length(parsName),length(OUTPUTNames));
for ii =1:size(AreaDat,2) % Reformat PRCC result into array for plotting
    AreaDat(:,ii)=S_OUTPUTS(ii).prcc(1:end-1, end);
end
f = figure('DefaultAxesFontSize', 14);
heatmap(OUTPUTNames, parsName,...
abs(AreaDat),'Colormap',jet); % Plot the absolute value of the correlation matrix
title(tname);
exportgraphics(f,fname, 'resolution', 300);

end


% Model Code
% Parameters of the antiPCKS9 model
function p =  antiPCSK9modelPars()
p.LDLparticleCE  =  0.92;      % milligram/nanomole 
p.circ_volume  = 5;         % liter              
p.clearance_hepatic_fraction  = 0.8;                          
p.baseline_hepatic_cholesterol   = 6000; %      milligram   ??       
p.maxSREBP2level  = 2 ;                            % ??
%p.minSREBP2level = 0;                             % ??
p.LDLcClearanceRate =0.187; %1/day              
p.deciliter_to_liter = 10;       %deciliter/liter    
p.dilipidemic_index = 1;                            
p.PK_Ka  =  0.247; %     1/day              
p.PK_Kel_F  = 0.0443; %    1/day              
p. PK_V2_F  =  5640; %      milliliter         
p.PK_ComplexClearanceRate = 0.182; %     1/day              
p.PK_kon  = 0.84 ;     %1/nanomole/day     
p.PK_koff  = 1.11; %      1/day              
p.PK_V3  = 3990; %      milliliter         
p.PK_Q   = 615; %       milliliter/day     
p.Baselinepcsk9  = 281.94; %    nanogram/milliliter
p.LDLrIndClearanceRate = 0.02;       
p.AbsorptionFraction = 0.5;
p.CholesterolIntakeDiet = 300;
p.StatinEffectOnCholesterolSynthesis = 1;
p.BaselineCholesterolSynthesisRate = 800;
p.LDLparticleProdRate = 0.03;
p.LossRate = 0.23;
p.HDLcClearanceRate = 0.3; % / day
p.pcsk9SynthesisRate = 9.4488;
p.pcsk9_synthesis_Vm_up = 2;
p.pcsk9_synthesis_Km_up = 1.5;
p.pcsk9_synthesis_Vm_down = 0.7;
p.pcsk9_synthesis_Km_down = 0.5;
p.pcsk9ClearanceRate = 2.48;
p.LDLr0 = 1;
%p.gamma = 0;
p.LDLrSynthesis = 0.5;
p.LDLr_expression_Vm_up = 3;
p.LDLr_expression_Km_up = 1.5;
p.LDLr_expression_Vm_down = 0.7;
p.LDLr_expression_Km_down = 0.5;
p.LDLrClearance = 0.5;
p.pcsk9_on_LDLr = 0.75;
p.pcsk9_on_LDLr_range = 650;
p.SREBP2 = 1;
p.circ_pcsk9_ngperml = 281.94;
p.HDLch = 50;
end

function IC = getInitialConditions()
hepatic_cholesterol = 6000;
LDLc = 4000;
circ_pcsk9 = 3.81;
surface_LDLr = 1;
antipcsk9 = 0.0001;
antipcsk9_dose =400000;
complex = 0;
peripheral = 0;
IC = [hepatic_cholesterol, LDLc, circ_pcsk9, surface_LDLr, ...
    antipcsk9, antipcsk9_dose, complex, peripheral]; %...
end

function dydt = odefun(t,y,p)
dydt = zeros(length(y), 1);
hepatic_cholesterol = y(1);
LDLc = y(2);
circ_pcsk9 = y(3);
surface_LDLr = y (4);
antipcsk9 = y(5);
antipcsk9_dose = y(6);
complex = y(7) ;
peripheral = y(8);

% repeated assignment

SREBP2 = transform(hepatic_cholesterol,...
    cat(2, p.maxSREBP2level, 0), ...
    p.baseline_hepatic_cholesterol, 3);

circ_pcsk9_ngperml = circ_pcsk9 * 74;
    
% Fluxes
DietCholesterolAbsorption = p.AbsorptionFraction * p.CholesterolIntakeDiet;

CholesterolSynthesis = p.BaselineCholesterolSynthesisRate...
    * p.StatinEffectOnCholesterolSynthesis;

LDLFormation = p.LDLparticleProdRate * hepatic_cholesterol...
    * p.circ_volume * p.LDLparticleCE;

LDLcleranceToHepatic = p.LDLcClearanceRate * p.dilipidemic_index *...
    LDLc * surface_LDLr * p.clearance_hepatic_fraction;

CholesterolLost = p.LossRate * hepatic_cholesterol;

HDLclerance = p.HDLcClearanceRate * p.HDLch * p.circ_volume * p.deciliter_to_liter...
    * p.clearance_hepatic_fraction;

LDLcleranceToPeriphery = p.LDLcClearanceRate * p.dilipidemic_index * LDLc * (1- ...
    p.clearance_hepatic_fraction) *surface_LDLr;


pcsk9_synthesis = p.pcsk9SynthesisRate* SREBP2_reg(SREBP2,...
    p.pcsk9_synthesis_Vm_up,...
    p.pcsk9_synthesis_Km_up,...
    p.pcsk9_synthesis_Vm_down,...
    p.pcsk9_synthesis_Km_down);

pcsk9_clearance = p.pcsk9ClearanceRate * circ_pcsk9 ...
    * (surface_LDLr / p.LDLr0)^0;

LDLr_expression = p.LDLrSynthesis * SREBP2_reg(SREBP2, ...
    p.LDLr_expression_Vm_up, p.LDLr_expression_Km_up, ...
    p.LDLr_expression_Vm_down, p.LDLr_expression_Km_down);

LDLr_clearance = p.LDLrClearance * surface_LDLr * ...
    transform(circ_pcsk9_ngperml,...
    cat(2,(1 - p.pcsk9_on_LDLr), (1 + p.pcsk9_on_LDLr)), ...
    p.Baselinepcsk9, p.pcsk9_on_LDLr_range,'lin');


absorption = ((p.PK_Ka * antipcsk9_dose / p.PK_V2_F) * (1000 / 150) ) ;

binding = p.PK_kon * antipcsk9 * circ_pcsk9 - p.PK_koff * complex;

antipcsk9_clearance = p.PK_Kel_F * antipcsk9;

complex_clearance = p.PK_ComplexClearanceRate * complex;

dose_compartment_clearance = p.PK_Ka * antipcsk9_dose;

distribution = p.PK_Q / p.PK_V2_F * (antipcsk9 - peripheral);

redistribution = p.PK_Q / p.PK_V3 * (peripheral - antipcsk9);

LDLrIndependentCleranceToHepatic = p.LDLrIndClearanceRate * ...
    p.clearance_hepatic_fraction * LDLc;

LDLrIndependentCleranceToPeriphery = p.LDLrIndClearanceRate * ...
    (1 - p.clearance_hepatic_fraction) * LDLc;


% ODEs
dhepatic_cholesterol_dt = DietCholesterolAbsorption + ...
    CholesterolSynthesis - LDLFormation + LDLcleranceToHepatic - ...
    CholesterolLost + HDLclerance + LDLrIndependentCleranceToHepatic;

dLDLc_dt = LDLFormation - LDLcleranceToHepatic - ...
    LDLcleranceToPeriphery - LDLrIndependentCleranceToHepatic - ...
    LDLrIndependentCleranceToPeriphery;

dcirc_pcsk9_dt = pcsk9_synthesis - pcsk9_clearance - binding;

dsurface_LDLr_dt = (LDLr_expression - LDLr_clearance);

dantipcsk9_dt = absorption - binding - antipcsk9_clearance - distribution;

dantipcsk9_dose_dt = -dose_compartment_clearance;

dcomplex_dt = binding - complex_clearance;

dperipheral_dt = -redistribution;


dydt(1) = dhepatic_cholesterol_dt;
dydt(2) = dLDLc_dt;
dydt(3) = dcirc_pcsk9_dt;
dydt(4) = dsurface_LDLr_dt;
dydt(5) = dantipcsk9_dt;
dydt(6) = dantipcsk9_dose_dt;
dydt(7) = dcomplex_dt ;
dydt(8) = dperipheral_dt;

end

function [Yout] = SREBP2_reg(SREBP2,Vmax_up, Km_up, Vmax_down, Km_down)

if ( SREBP2 >= 1 )
    In = (SREBP2-1);
    Yout = 1 + (Vmax_up-1) * In / (In + (Km_up-1)); 
end

if ( SREBP2 < 1 )
    In = (1-SREBP2);
    Yout = 1 - Vmax_down * In / (In + Km_down); 
end
 
end

function [Yout] = transform(Xin,Yrange,Xic50,varargin)

nVarargs = length(varargin);
Xscale='log'; % default is log
Xrange=3;     % log-order if log; actual range if linear
if (nVarargs==1) 
    Xrange=varargin{1}; 
end
if (nVarargs==2)
    Xrange=varargin{1}; Xscale=varargin{2};
end

if strcmp(Xscale,'log')
    expT=3.97865*Xrange^-0.995; % This is fitted to give 99% output at the extreme of the range specified
    Yout=( Yrange(2) - Yrange(1) ) * (Xin^expT/(Xin^expT+Xic50^expT)) + Yrange(1);
end

 if strcmp(Xscale,'lin')
    Xmin=Xic50-Xrange/2;Xmax=Xic50+Xrange/2;
    expT=3.97865*2^-0.995; % optimized for range of -1 to 1;
    Xmod=10^( 2*(Xin-Xmin)/(Xmax-Xmin) -1 )  ; 
    Xmodic50=10^( 2*(Xic50-Xmin)/(Xmax-Xmin) -1 );  
    Yout=( Yrange(2) - Yrange(1) ) * (Xmod^expT/(Xmod^expT+Xmodic50^expT)) + Yrange(1);
 end
 
end