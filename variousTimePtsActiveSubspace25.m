%Comment Code and Put in Google Drive
rng(27);
%USER EDIT: Choose numbers of points to sample
M = 50;

%Gets the parameters of the PCSK9 model
parameters = antiPCSK9modelPars();
fields = fieldnames(parameters);
P = struct2cell(parameters);

ScaledPoints = ScaleCreate(M, P);

%Set up for complex step gradient approximation
h = 1e-14;
gradApproxs = cell(27,1);
G = cell(27,1);
U = cell(27,1);
S = cell(27,1);
V = cell(27,1);
DS = cell(27,1);
for i = 1:27
    gradApproxs{i} = zeros(length(P), M);
end

%Runs model for, iteratively using a complex step derviative approximation
%at each parameter for each row (set of input points)
    for i = 1:M
        for j = 1:length(P)
            P{j} = complex((P{j}*ScaledPoints(i,j)), h);
            parameters = cell2struct(P, fields);
            IC = getInitialConditions();
            t0 = 0; td1 = 28; td2 = 56; td3 = 84; td4 = 112; tfinal = 140;
            [T1,Y1] = ode15s(@(t, y)odefun(t, y, parameters),linspace(t0, td1, 280),IC);
            for iii = 1:8
                IC(iii) = Y1(end, iii);
            end
            IC(6) = 400000;
            [T2,Y2] = ode15s(@(t, y)odefun(t, y, parameters),linspace(td1, td2, 280),IC);
            for iii = 1:8
                IC(iii) = Y2(end, iii);
            end
            IC(6) = 400000;
            [T3,Y3] = ode15s(@(t, y)odefun(t, y, parameters),linspace(td2, td3, 280),IC);
            for iii = 1:8
                IC(iii) = Y3(end, iii);
            end
            IC(6) = 400000;
            [T4,Y4] = ode15s(@(t, y)odefun(t, y, parameters),linspace(td3, td4, 280),IC);
            for iii = 1:8
                IC(iii) = Y4(end, iii);
            end
            IC(6) = 400000;
            [T5,Y5] = ode15s(@(t, y)odefun(t, y, parameters),linspace(td4, tfinal, 280),IC);
            Y = [Y1; Y2; Y3; Y4; Y5];
            T = [T1; T2; T3; T4; T5];
            for q = 1:27
                gradApproxs{q}(j, i) = imag(Y((50*q),2)) / h; %USER EDIT: Desired output %FIX HERE
            end
            P{j} = real(P{j});
            P{j} = P{j} / ScaledPoints(i,j);
        end
    end

%Forms G as described in Active Subspace Paper
for i = 1:27
    G{i} = (1/sqrt(M)) * gradApproxs{i};
end

%Take SVD of G (Note: Singular values in
%S tell us which columns of U are most important)
for i = 1:27
    [U{i},S{i},V{i}] = svd(G{i});
end
%Number of active subspaces: Choose using S matrix
for i = 1:27
    DS{i} = diag(S{i});
end

FULLactivityScores = cell(1, 27);

for q = 1:27
    activityScores = zeros(length(P), 1);
    for i = 1:length(P)
        for j = 1:25
            activityScores(i) = activityScores(i) + DS{q}(j) * (U{q}(i,j))^2;
        end
    end
    FULLactivityScores{q} = activityScores;
end

plotfields = string(fieldnames(parameters))';

ltFAS = cell2mat(FULLactivityScores);
ltFAS = ltFAS + 1;
ltFAS = log(ltFAS);

figure;
pcolor(ltFAS);
colorbar;
xlabel('Time (5 Days)', 'FontSize', 22)
ylabel('Parameter', 'FontSize', 22)
title('Active Subspace Results at Different Time Points', 'FontSize', 26)
ax = gca;
properties(ax);
ax.LineWidth = 1.2;
set(ax, 'FontSize', 18)

% Parameters of the antiPCKS9 model
function p =  antiPCSK9modelPars()
p.LDLparticleCE  =  0.92;      % milligram/nanomole 
p.circ_volume  = 5;         % liter              
p.clearance_hepatic_fraction  = 0.8;                          
p.baseline_hepatic_cholesterol   = 6000; %      milligram   ??       
p.maxSREBP2level  = 2 ;                            % ??
p.minSREBP2level = 0;                             % ??
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
p.gamma = 0;
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
    cat(2, p.maxSREBP2level, p.minSREBP2level), ...
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
    * (surface_LDLr / p.LDLr0)^p.gamma;

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