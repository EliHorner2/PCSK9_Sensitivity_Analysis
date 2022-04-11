addpath('functions');
rng(27);
M = 100; %Should be about 10 * dimensions
SigPars = 20;
%[V,I] = maxk(FAS25, SigPars);
I = [24, 11, 7, 35, 25, 10, 40, 19, 41, 22, 5, 26, 3, 32, 27, 9, 1, 13, 34, 20];
parameters = antiPCSK9modelPars();
fields = fieldnames(parameters);
P = struct2cell(parameters);

%5 percent for unknowns, 0 for held constants
pmins = zeros(length(P), 1);
pmaxs = zeros(length(P), 1);
pmins(1) = 0.99;
pmaxs(1) = 1.01;
pmins(2) = 0.99;
pmaxs(2) = 1.01;
pmins(3) = 0.99;
pmaxs(3) = 1.01;
pmins(4) = 0.99;
pmaxs(4) = 1.01;
pmins(5) = 0.99;
pmaxs(5) = 1.01;
pmins(6) = 0.99;
pmaxs(6) = 1.01;
pmins(7) = 0.65;
pmaxs(7) = 1.35;
pmins(8) = 0.99;
pmaxs(8) = 1.01;
pmins(9) = 0.95;
pmaxs(9) = 1.05;
pmins(10) = 0.48;
pmaxs(10) = 1.52;
pmins(11) = 0.95;
pmaxs(11) = 1.05;
pmins(12) = 0.70;
pmaxs(12) = 1.30;
pmins(13) = 0.99;
pmaxs(13) = 1.01;
pmins(14) = 0.99;
pmaxs(14) = 1.01;
pmins(15) = 0.99;
pmaxs(15) = 1.01;
pmins(16) = 0.80;
pmaxs(16) = 1.20;
pmins(17) = 0.99;
pmaxs(17) = 1.01;
pmins(18) = 0.99;
pmaxs(18) = 1.01;
pmins(19) = 0.99;
pmaxs(19) = 1.01;
pmins(20) = 0.99;
pmaxs(20) = 1.01;
pmins(21) = 0.80;
pmaxs(21) = 1.20;
pmins(22) = 0.75;
pmaxs(22) = 1.25;
pmins(23) = 0.90;
pmaxs(23) = 1.10;
pmins(24) = 0.55;
pmaxs(24) = 1.45;
pmins(25) = 0.95;
pmaxs(25) = 1.05;
pmins(26) = 0.99;
pmaxs(26) = 1.01;
pmins(27) = 0.82;
pmaxs(27) = 1.18;
pmins(28) = 0.95;
pmaxs(28) = 1.05;
pmins(29) = 0.95;
pmaxs(29) = 1.05;
pmins(30) = 0.95;
pmaxs(30) = 1.05;
pmins(31) = 0.95;
pmaxs(31) = 1.05;
pmins(32) = 0.81;
pmaxs(32) = 1.19;
pmins(33) = 0.95;
pmaxs(33) = 1.05;
pmins(34) = 0.95;
pmaxs(34) = 1.05;
pmins(35) = 0.55;
pmaxs(35) = 1.45;
pmins(36) = 0.95;
pmaxs(36) = 1.05;
pmins(37) = 0.95;
pmaxs(37) = 1.05;
pmins(38) = 0.95;
pmaxs(38) = 1.05;
pmins(39) = 0.95;
pmaxs(39) = 1.05;
pmins(40) = 0.55;
pmaxs(40) = 1.45;
pmins(41) = 0.95;
pmaxs(41) = 1.05;
pmins(42) = 0.95;
pmaxs(42) = 1.05;
pmins(43) = 0.95;
pmaxs(43) = 1.05;
pmins(44) = 0.95;
pmaxs(44) = 1.05;
pmins(45) = 0.95;
pmaxs(45) = 1.05;

%Start Loop
Loops = 100;
RP = cell(Loops, 1);
for ii = 1:Loops
    LatinPoints = lhsdesign(M, SigPars);
    ScaledPoints = zeros(M, SigPars);
    Outs = cell(M, 1);

    for i = 1:M
        for j = 1:SigPars
            ScaledPoints(i, j) = (((pmaxs(I(j)) - pmins(I(j))) * LatinPoints(i, j)) + pmins(I(j)));
        end
    end

    for i = 1:M
        for j = 1:SigPars
            P{I(j)} = P{I(j)} * ScaledPoints(i, j);
        end
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
        Outs{i} = Y(550, 2);% Y(150,3) Y(150, 5)];
        for j = 1:SigPars
            P{I(j)} = P{I(j)} / ScaledPoints(i, j);
        end
    end

    Outs = cell2mat(Outs);
    %y = zeros(1,M);
    %for k = 1:M
        %y(k) = Outs{k};
    %end

    options.nugget_est=true;
    model = ppgasp(LatinPoints, Outs, options);
    RP{ii} = model.range_par;
end


RPMat = cell2mat(RP');

COPY = RPMat;

for i = 1:SigPars
    for j = 1:100
        if COPY(i,j) > 20
            COPY(i,j) = 20;
        end
    end
end

Medians = zeros(SigPars, 1);
for i = 1:SigPars
    Medians(i) = median(COPY(i,:));
end

Means = zeros(SigPars, 1);
for i = 1:SigPars
    Means(i) = mean(COPY(i,:));
end

[VV, II] = mink(Medians, SigPars);
[VVV, III] = mink(Means, SigPars);
% figure;
% tiledlayout(2,2);
% nexttile
% histogram(COPY(1,:), 50)
% hold on
% xline(Medians(1), '-', 'Median', 'Color', 'r')
% title(I(1))
% nexttile
% histogram(COPY(2,:), 50)
% hold on
% xline(Medians(2), '-', 'Median', 'Color', 'r')
% title(I(2))
% nexttile
% histogram(COPY(3,:), 50)
% hold on
% xline(Medians(3), '-', 'Median', 'Color', 'r')
% title(I(3))
% nexttile
% histogram(COPY(4,:), 50)
% hold on
% xline(Medians(4), '-', 'Median', 'Color', 'r')
% title(I(4))
% nexttile
% histogram(COPY(5,:), 50)
% hold on
% xline(Medians(5), '-', 'Median', 'Color', 'r')
% title(I(5))
% nexttile
% histogram(COPY(6,:), 50)
% hold on
% xline(Medians(6), '-', 'Median', 'Color', 'r')
% title(I(6))
% nexttile
% histogram(COPY(7,:), 50)
% hold on
% xline(Medians(7), '-', 'Median', 'Color', 'r')
% title(I(7))
% nexttile
% histogram(COPY(8,:), 50)
% hold on
% xline(Medians(8), '-', 'Median', 'Color', 'r')
% title(I(8))
% nexttile
% histogram(COPY(9,:), 50)
% hold on
% xline(Medians(9), '-', 'Median', 'Color', 'r')
% title(I(9))

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