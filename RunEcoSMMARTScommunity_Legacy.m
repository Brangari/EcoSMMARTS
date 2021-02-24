
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% % % % % % % % % EcoSMMARTS_Legacy()

% % % % EcoSMMARTS - Community version simulation (Research paper: Brangarí et al, 2021, SBB)
% % % % In this updated version of the model, a higher degree of microbial diversity is included (functional groups with distinct life strategies, and competition for resources)
% % % % Model is used to study how communities respond to soil rewetting and the dependence of these responses on the history of soil moisture and water stress
% % % % By Albert C Brangarí. Lund University, Department of Biology, 2020.

close  all;
fprintf('\n \n Starting EcoSMMARTS simulation...')

%% % % % % % % % % Parameter definition
    
    % % % % C compartments: POC (particulate organic matter), DOC (dissolved organic carbon), AC (active cells), EPS (extracellular polymeric substances), EZ (enzymes), DC (dormant cells), OSac (osmolytes in AC), OSdc (osmolytes in DC), CR (cell residues or "zombie cells")    
  	% % % % This simulation uses two microbial functional groups: fast growers (FG) vs drought resistors (DR)
    % % % % FG: with higher growth rates and CUE, slower production of osmolytes, larger decay and death rates under osmotic stress, and larger drought-legacy effects (‘memory’)
    % % % % DR: lower growth and death rates, less efficient use of C, more efficient osmoregulation, and shorter legacy effects
    % % % % Parameters for FG are defined first in vectors (1st position), DR are defined second (2nd position in vectors)
    
    vec = [0 1 2 0 1 2 2 2 1];     % Number of pools of each C compartment [POC DOC AC EPS EZ DC OSac OSdc CR]. 0 for deactivation (Note that all combinations possible are not implemented in this version)
    
    POC_in = 15;                % Initial concentration of POC [mg/cm^3] -> Brangarí et al. (2018)
    DOC_in = 0.05/0.0946;       % Initial mass of DOC [mg/cm^3] -> Lawrence et al. (2009) (concentration / VWC)
    EE_in = 6e-3;               % Initial concentration of EE [mg/cm^3] -> Brangarí et al. (2018)
    AC_in = [7e-1 7e-1]'/2;     % Initial concentration of AC [mg/cm^3] -> Lawrence et al. (2009) (distributed among functional groups)
    DC_in = AC_in;              % Initial concentration of DC [mg/cm^3] -> Blagodatskaya and Kuzyakov (2013)
    CR_in = 0;                  % Initial concentration of CR [mg/cm^3] -> Assumed
    OSac_in = [365.7 365.7]';   % Initial concentration of OSac [mg/cm^3] -> Estimated assuming equilibrium with initial water potential
    OSdc_in = OSac_in;          % Initial concentration of OSdc [mg/cm^3] -> Estimated assuming equilibrium with initial water potential
    ThetaS = 0.32;              % Saturated water content [-] -> Brangarí et al. (2020)
    ThetaR = 0;                 % Residual water content [-] -> Brangarí et al. (2020)
    a = 0.0172;                 % Empirical coefficient WRC [1/cm] -> Brangarí et al. (2020)
    n = 1.3824;                 % Empirical coefficient WRC [1/cm] -> Brangarí et al. (2020)
    Gamma = 1.8;                % Factor of tortuosity [-] -> Hamamoto et al. (2010)
    
    ADD = 0;                        % Litter input [mg/cm^3] -> Assumed at steady-state
    Lambda_EE = 0.01;               % Coefficient of carbon allocation diverted towards enzymes [-] -> Brangarí et al. (2018)
    Y_M = [0.7 0.5]';               % y: Yield coefficient of carbon uptake [-] -> Larger for FG, smaller for DR
    MuCA = [3.5e-3 1e-3]';          % Maximum specific uptake rate [1/min]
    MuCZ = MuCA(1)/10;              % Maximum specific mineralization rate by cell residues [1/min] -> From calibration
    NuPOC = 1.7e-5;                 % Maximum specific POC decomposition rate [1/min] -> From calibration
    Yos = 0.6;                      % Yield coefficient of OS elimination [-] -> Brangarí et al. (2020) (if =0: no assimilation -> all to respiration)
    Zre = 1e-2;                     % Yield coefficient of reactivation [-] -> Brangarí et al. (2020)
    multD = 3;                      % Coefficient of aggregate disruption [-] -> Calibrated
    K_C = 5;                        % Half-saturation constant of DOC [mg/cm^3] (never set to 0) -> Brangarí et al. (2020)
    K_EE = EE_in/10;                % Half-saturation constant of EZ [mg/cm^3] (never set to 0) -> Assumed
    KdAC = [1e-5 5e-6]';            % Constant rate of decay for AC [1/min]
    KdACs = [1e-4 1e-6]';           % Constant rate of decay for AC by stress [1/min]
    KdDC = KdAC/10;                 % Constant rate of decay for DC [1/min] -> Assumed
    KdEZ = 1e-6;                    % Constant rate of decay for EZ [1/min] -> Assumed
    KdCR = KdAC(1)*10;              % Constant rate of decay for CR [1/min] -> Assumed
    Lambda_r = 0.5;                 % Coefficient of recycling [-] -> Assumed
    Lambda_z = 0.8;                 % Coefficient of cell residues permanence [-] -> Assumed
    Tau_i = [2e-6 2e-6]';           % Maximum specific inactivation rate [1/min] -> Calibrated
    Tau_r = 10*Tau_i./[1.43 1.43]';	% Maximum specific reactivation rate [1/min] -> Konopka (2000)
    Tau_OS = [1e-3 2e-3]';          % Maximum specific rate of osmoregulation [1/min] -> Calibrated
    tMemS = 10*24*60;               % Coefficients of memory length for soil [min] -> Brangarí et al. (2020)
    tMemB = [2 1]'*24*60;           % Coefficients of memory length for microbes [min] -> Calibrated
    RelCinlet = 0;                  % Concentration of DOC and EZ in the inlet water (0: if deonized/rainfall water, 1: water carring compounds (same concentration as in control volume))
    
%% % % % % % % % % Run simulations
    
    DT = 30;   % Delta t results [min]
    fprintf('\n \n Exposing soils to the different soil moisture regimes...')
    fprintf('\n \n (drying and rewetting virtual samples requires some time)')
    
% % Initialization phase: optimal conditions
    TimeD = 45;         % Length of the first part of the pretreatment [Optimal WHC] [d]
    modeRW = 'C';       % Mode used to define moisture. C: no rewetting (constant initial water potential), RW-D: rewetting followed by drying (maximum saturation and decrease), RW-C: rewetting followed by contant conditions
    Whc0 = 50;          % Pretreatment water holding capacity [%]
    FirstRW = 10;       % Time of the first RW event [d]
    WhcRW = NaN;        % Water holding capacity after rewetting [%]
    Dt_RW = NaN;        % Interval between RW events [d]
    ChT_FC = NaN;       % Characteristic time of drainage (only relevant after rainfall events) [d]
    ChT_0 = NaN;        % Characteristic time of evaporation (only relevant after rainfall events) [d]
    T_hist = 0;         % Moisture history of soils (previous treatment): time [d]
    Whc_hist = Whc0;    % Moisture history of soils (previous treatment): WHC [%]
    ThetaFC_in = 0.157856735938127/2; % Initial field capacity -> From Brangarí et al. (2020)
    [POC_0,DOC_0,AC_0,DC_0,EZ_0,OSac_0,OSdc_0,CR_0,Psi_0,ThetaW_0,ThetaFC,RespAC_0,RespCR_0,RespOS_0,RespDC_0,Growth_0,MassTot_0,MassPart_0,TimeDay_0] = EcoSMMARTS_community(TimeD,DT,modeRW,Whc0,FirstRW,WhcRW,Dt_RW,ChT_FC,ChT_0,T_hist,Whc_hist,...
                                                                                                         ThetaS,ThetaR,a,n,Gamma,POC_in,DOC_in,EE_in,AC_in,DC_in,CR_in,OSac_in,OSdc_in,[0 0 0 0 ThetaFC_in*DOC_in],[0 0 0 0 ThetaFC_in*DOC_in],...
                                                                                                         ADD,Lambda_EE,Y_M,MuCA,MuCZ,NuPOC,Yos,Zre,multD,K_C,K_EE,KdAC,KdACs,...
                                                                                                         KdDC,KdEZ,KdCR,Lambda_r,Lambda_z,Tau_i,Tau_r,Tau_OS,tMemS,tMemB,RelCinlet,vec);
% % Moisture variation phase: Control (CONTROL) - Optimal WHC
    TimeD = 365-TimeD;	% Length of the second part [d]
    modeRW = 'C';
    Whc0 = ThetaW_0(end)*100/ThetaFC;  % Pretreatment water holding capacity [%]
    FirstRW = NaN;      % Time of the first RW event [d]
    WhcRW = NaN;        % Water holding capacity after rewetting [%]
    Dt_RW = NaN;        % Interval between RW events [d]
    ChT_FC = NaN;       % Characteristic time of drainage (only relevant after rainfall events) [d]
    ChT_0 = NaN;        % Characteristic time of evaporation (only relevant after rainfall events) [d]
    T_hist = TimeDay_0; % Moisture history of soils (previous treatment): time [d]
    Whc_hist = ThetaW_0*100/ThetaFC;  % Moisture history of soils (previous treatment): WHC [%]
    [POC_0C,DOC_0C,AC_0C,DC_0C,EZ_0C,OSac_0C,OSdc_0C,CR_0C,Psi_0C,ThetaW_0C,~,RespAC_0C,RespCR_0C,RespOS_0C,RespDC_0C,Growth_0C,MassTot_0C,MassPart_0C,TimeDay_0C] = EcoSMMARTS_community(TimeD,DT,modeRW,Whc0,FirstRW,WhcRW,Dt_RW,ChT_FC,ChT_0,T_hist,Whc_hist,...
                                                                                                         ThetaS,ThetaR,a,n,Gamma,POC_0(end),DOC_0(end),EZ_0(end),AC_0(:,end),DC_0(:,end),CR_0(end),OSac_0(:,end),OSdc_0(:,end),MassTot_0(:,end),MassPart_0(:,end),...
                                                                                                         ADD,Lambda_EE,Y_M,MuCA,MuCZ,NuPOC,Yos,Zre,multD,K_C,K_EE,KdAC,KdACs,...
                                                                                                         KdDC,KdEZ,KdCR,Lambda_r,Lambda_z,Tau_i,Tau_r,Tau_OS,tMemS,tMemB,RelCinlet,vec);
% % Moisture variation phase: Severe scenario (SEV)
    modeRW = 'RW-D';
    Whc0 = ThetaW_0(end)*100/ThetaFC;   % Pretreatment water holding capacity [%]
    FirstRW = 0;    % Time of the first RW event [d]
    WhcRW = 4;      % Water holding capacity after rewetting [%]
    Dt_RW = 40.001; % Interval between RW events [d]
    ChT_FC = 1;     % Characteristic time of drainage (only relevant after rainfall events) [d]
    ChT_0 = 20;     % Characteristic time of evaporation (only relevant after rainfall events) [d]
    T_hist = TimeDay_0;                 % Moisture history of soils (previous treatment): time [d]
    Whc_hist = ThetaW_0*100/ThetaFC;    % Moisture history of soils (previous treatment): WHC [%]
    [POC_0M,DOC_0M,AC_0M,DC_0M,EZ_0M,OSac_0M,OSdc_0M,CR_0M,Psi_0M,ThetaW_0M,~,RespAC_0M,RespCR_0M,RespOS_0M,RespDC_0M,Growth_0M,MassTot_0M,MassPart_0M,TimeDay_0M] = EcoSMMARTS_community(TimeD,DT,modeRW,Whc0,FirstRW,WhcRW,Dt_RW,ChT_FC,ChT_0,T_hist,Whc_hist,...
                                                                                                         ThetaS,ThetaR,a,n,Gamma,POC_0(end),DOC_0(end),EZ_0(end),AC_0(:,end),DC_0(:,end),CR_0(end),OSac_0(:,end),OSdc_0(:,end),MassTot_0(:,end),MassPart_0(:,end),...
                                                                                                         ADD,Lambda_EE,Y_M,MuCA,MuCZ,NuPOC,Yos,Zre,multD,K_C,K_EE,KdAC,KdACs,...
                                                                                                         KdDC,KdEZ,KdCR,Lambda_r,Lambda_z,Tau_i,Tau_r,Tau_OS,tMemS,tMemB,RelCinlet,vec);
% % Moisture variation phase: Moderate scenario (MOD)
    modeRW = 'RW-D';
    Whc0 = ThetaW_0(end)*100/ThetaFC;   % Pretreatment water holding capacity [%]
    FirstRW = 0;  	% Time of the first RW event [d]
    WhcRW = 8;      % Water holding capacity after rewetting [%]
    Dt_RW = 20.001; % Interval between RW events [d]
    ChT_FC = 1;     % Characteristic time of drainage (only relevant after rainfall events) [d]
    ChT_0 = 20;     % Characteristic time of evaporation (only relevant after rainfall events) [d]
    T_hist = TimeDay_0;                 % Moisture history of soils (previous treatment): time [d]
    Whc_hist = ThetaW_0*100/ThetaFC;    % Moisture history of soils (previous treatment): WHC [%]
    [POC_0H,DOC_0H,AC_0H,DC_0H,EZ_0H,OSac_0H,OSdc_0H,CR_0H,Psi_0H,ThetaW_0H,~,RespAC_0H,RespCR_0H,RespOS_0H,RespDC_0H,Growth_0H,MassTot_0H,MassPart_0H,TimeDay_0H] = EcoSMMARTS_community(TimeD,DT,modeRW,Whc0,FirstRW,WhcRW,Dt_RW,ChT_FC,ChT_0,T_hist,Whc_hist,...
                                                                                                         ThetaS,ThetaR,a,n,Gamma,POC_0(end),DOC_0(end),EZ_0(end),AC_0(:,end),DC_0(:,end),CR_0(end),OSac_0(:,end),OSdc_0(:,end),MassTot_0(:,end),MassPart_0(:,end),...
                                                                                                         ADD,Lambda_EE,Y_M,MuCA,MuCZ,NuPOC,Yos,Zre,multD,K_C,K_EE,KdAC,KdACs,...
                                                                                                         KdDC,KdEZ,KdCR,Lambda_r,Lambda_z,Tau_i,Tau_r,Tau_OS,tMemS,tMemB,RelCinlet,vec);

%% % % % % % % % % Drying of soils pre-experiment

    TimeD = 35;     % Length of the drying pre-experiment [d]
    
% % CONTROL
    modeRW = 'C';
    Whc0 = ThetaW_0C(end)*100/ThetaFC;  % Pretreatment water holding capacity [%]
    FirstRW = NaN;  % Time of the first RW event [d]
    WhcRW = NaN;    % Water holding capacity after rewetting [%]
    Dt_RW = NaN;    % Interval between RW events [d]
    ChT_FC = NaN;   % Characteristic time of drainage (only relevant after rainfall events) [d]
    ChT_0 = NaN;    % Characteristic time of evaporation (only relevant after rainfall events) [d]
    T_hist = [TimeDay_0 TimeDay_0(end)+TimeDay_0C];     % Moisture history of soils (previous treatment): time [d]
    Whc_hist = [ThetaW_0 ThetaW_0C]*100/ThetaFC;        % Moisture history of soils (previous treatment): WHC [%]
    [POC_dC,DOC_dC,AC_dC,DC_dC,EZ_dC,OSac_dC,OSdc_dC,CR_dC,Psi_dC,ThetaW_dC,~,RespAC_dC,RespCR_dC,RespOS_dC,RespDC_dC,Growth_dC,MassTot_dC,MassPart_dC,TimeDay_dC] = EcoSMMARTS_community(TimeD,DT,modeRW,Whc0,FirstRW,WhcRW,Dt_RW,ChT_FC,ChT_0,T_hist,Whc_hist,...
                                                                                                         ThetaS,ThetaR,a,n,Gamma,POC_0C(end),DOC_0C(end),EZ_0C(end),AC_0C(:,end),DC_0C(:,end),CR_0C(end),OSac_0C(:,end),OSdc_0C(:,end),MassTot_0C(:,end),MassPart_0C(:,end),...
                                                                                                         ADD,Lambda_EE,Y_M,MuCA,MuCZ,NuPOC,Yos,Zre,multD,K_C,K_EE,KdAC,KdACs,...
                                                                                                         KdDC,KdEZ,KdCR,Lambda_r,Lambda_z,Tau_i,Tau_r,Tau_OS,tMemS,tMemB,RelCinlet,vec);
% % SEV
    modeRW = 'RW-D';
    Whc0 = ThetaW_0M(end)*100/ThetaFC;      % Pretreatment water holding capacity [%]
    FirstRW = 0;    % Time of the first RW event [d]
    WhcRW = 4;      % Water holding capacity after rewetting [%]
    Dt_RW = NaN;    % Interval between RW events [d]
    ChT_FC = 1;     % Characteristic time of drainage (only relevant after rainfall events) [d]
    ChT_0 = 20;     % Characteristic time of evaporation (only relevant after rainfall events) [d]
    T_hist = [TimeDay_0 TimeDay_0(end)+TimeDay_0M];     % Moisture history of soils (previous treatment): time [d]
    Whc_hist = [ThetaW_0 ThetaW_0M]*100/ThetaFC;        % Moisture history of soils (previous treatment): WHC [%]
    [POC_dM,DOC_dM,AC_dM,DC_dM,EZ_dM,OSac_dM,OSdc_dM,CR_dM,Psi_dM,ThetaW_dM,~,RespAC_dM,RespCR_dM,RespOS_dM,RespDC_dM,Growth_dM,MassTot_dM,MassPart_dM,TimeDay_dM] = EcoSMMARTS_community(TimeD,DT,modeRW,Whc0,FirstRW,WhcRW,Dt_RW,ChT_FC,ChT_0,T_hist,Whc_hist,...
                                                                                                         ThetaS,ThetaR,a,n,Gamma,POC_0M(end),DOC_0M(end),EZ_0M(end),AC_0M(:,end),DC_0M(:,end),CR_0M(end),OSac_0M(:,end),OSdc_0M(:,end),MassTot_0M(:,end),MassPart_0M(:,end),...
                                                                                                         ADD,Lambda_EE,Y_M,MuCA,MuCZ,NuPOC,Yos,Zre,multD,K_C,K_EE,KdAC,KdACs,...
                                                                                                         KdDC,KdEZ,KdCR,Lambda_r,Lambda_z,Tau_i,Tau_r,Tau_OS,tMemS,tMemB,RelCinlet,vec);
% % MOD
    modeRW = 'RW-D';
    Whc0 = ThetaW_0H(end)*100/ThetaFC;      % Pretreatment water holding capacity [%]
    FirstRW = 0;    % Time of the first RW event [d]
    WhcRW = 4;      % Water holding capacity after rewetting [%]
    Dt_RW = NaN;    % Interval between RW events [d]
    ChT_FC = 1;     % Characteristic time of drainage (only relevant after rainfall events) [d]
    ChT_0 = 20;     % Characteristic time of evaporation (only relevant after rainfall events) [d]
    T_hist = [TimeDay_0 TimeDay_0(end)+TimeDay_0H];     % Moisture history of soils (previous treatment): time [d]
    Whc_hist = [ThetaW_0 ThetaW_0H]*100/ThetaFC;        % Moisture history of soils (previous treatment): WHC [%]
    [POC_dH,DOC_dH,AC_dH,DC_dH,EZ_dH,OSac_dH,OSdc_dH,CR_dH,Psi_dH,ThetaW_dH,~,RespAC_dH,RespCR_dH,RespOS_dH,RespDC_dH,Growth_dH,MassTot_dH,MassPart_dH,TimeDay_dH] = EcoSMMARTS_community(TimeD,DT,modeRW,Whc0,FirstRW,WhcRW,Dt_RW,ChT_FC,ChT_0,T_hist,Whc_hist,...
                                                                                                         ThetaS,ThetaR,a,n,Gamma,POC_0H(end),DOC_0H(end),EZ_0H(end),AC_0H(:,end),DC_0H(:,end),CR_0H(end),OSac_0H(:,end),OSdc_0H(:,end),MassTot_0H(:,end),MassPart_0H(:,end),...
                                                                                                         ADD,Lambda_EE,Y_M,MuCA,MuCZ,NuPOC,Yos,Zre,multD,K_C,K_EE,KdAC,KdACs,...
                                                                                                         KdDC,KdEZ,KdCR,Lambda_r,Lambda_z,Tau_i,Tau_r,Tau_OS,tMemS,tMemB,RelCinlet,vec);

% % Merging solution vectors of the pretreatment (and sum of multi-compounds)
    DOC_LsolC = [DOC_0 DOC_0C DOC_dC];
    DOC_LsolM = [DOC_0 DOC_0M DOC_dM];
    DOC_LsolH = [DOC_0 DOC_0H DOC_dH];
    POC_LsolC = [POC_0 POC_0C POC_dC];
    POC_LsolM = [POC_0 POC_0M POC_dM];
    POC_LsolH = [POC_0 POC_0H POC_dH];
    AC_LsolC = [AC_0 AC_0C AC_dC; sum([AC_0 AC_0C AC_dC],1)];
    AC_LsolM = [AC_0 AC_0M AC_dM; sum([AC_0 AC_0M AC_dM],1)];
    AC_LsolH = [AC_0 AC_0H AC_dH; sum([AC_0 AC_0H AC_dH],1)];
    DC_LsolC = [DC_0 DC_0C DC_dC; sum([DC_0 DC_0C DC_dC],1)];
    DC_LsolM = [DC_0 DC_0M DC_dM; sum([DC_0 DC_0M DC_dM],1)];
    DC_LsolH = [DC_0 DC_0H DC_dH; sum([DC_0 DC_0H DC_dH],1)];
    EZ_LsolC = [EZ_0 EZ_0C EZ_dC];
    EZ_LsolM = [EZ_0 EZ_0M EZ_dM];
    EZ_LsolH = [EZ_0 EZ_0H EZ_dH];
    OSac_LsolC = [OSac_0 OSac_0C OSac_dC; sum([OSac_0 OSac_0C OSac_dC].*[AC_0 AC_0C AC_dC],1)./sum([AC_0 AC_0C AC_dC],1)];
    OSac_LsolM = [OSac_0 OSac_0M OSac_dM; sum([OSac_0 OSac_0M OSac_dM].*[AC_0 AC_0M AC_dM],1)./sum([AC_0 AC_0M AC_dM],1)];
    OSac_LsolH = [OSac_0 OSac_0H OSac_dH; sum([OSac_0 OSac_0H OSac_dH].*[AC_0 AC_0H AC_dH],1)./sum([AC_0 AC_0H AC_dH],1)];
    OSdc_LsolC = [OSdc_0 OSdc_0C OSdc_dC; sum([OSdc_0 OSdc_0C OSdc_dC].*[DC_0 DC_0C DC_dC],1)./sum([DC_0 DC_0C DC_dC],1)];
    OSdc_LsolM = [OSdc_0 OSdc_0M OSdc_dM; sum([OSdc_0 OSdc_0M OSdc_dM].*[DC_0 DC_0M DC_dM],1)./sum([DC_0 DC_0M DC_dM],1)];
    OSdc_LsolH = [OSdc_0 OSdc_0H OSdc_dH; sum([OSdc_0 OSdc_0H OSdc_dH].*[DC_0 DC_0H DC_dH],1)./sum([DC_0 DC_0H DC_dH],1)];
    CR_LsolC = [CR_0 CR_0C CR_dC];
    CR_LsolM = [CR_0 CR_0M CR_dM];
    CR_LsolH = [CR_0 CR_0H CR_dH];
    ThetaW_LsolC = [ThetaW_0 ThetaW_0C ThetaW_dC];
    ThetaW_LsolM = [ThetaW_0 ThetaW_0M ThetaW_dM];
    ThetaW_LsolH = [ThetaW_0 ThetaW_0H ThetaW_dH];
    TimeDay_Lsol = [TimeDay_0 TimeDay_0(end)+TimeDay_0C TimeDay_0(end)+TimeDay_0C(end)+TimeDay_dC];
    MassPart_LsolC = [MassPart_0 MassPart_0C MassPart_dC];
    MassPart_LsolM = [MassPart_0 MassPart_0M MassPart_dM];
    MassPart_LsolH = [MassPart_0 MassPart_0H MassPart_dH];
    MassTot_LsolC = [MassTot_0 MassTot_0C MassTot_dC];
    MassTot_LsolM = [MassTot_0 MassTot_0M MassTot_dM];
    MassTot_LsolH = [MassTot_0 MassTot_0H MassTot_dH];
    TimeT = TimeDay_Lsol(end);

% % Plot data (pretreatment)
    colorC = [0 0 0];
    colorH= [0.8 0 0];
    colorM = [0 0.6 0];
    colorR = [1 0.6 0.1];
    colorK = [0.4 0.7 1];
    sizeF = 20;
    
figure(3); subplot(3,2,1); hold on; box on; set(gca,'fontsize',sizeF); subplot(3,2,2); hold on; box on; set(gca,'fontsize',sizeF); subplot(3,2,3); hold on; box on; set(gca,'fontsize',sizeF); subplot(3,2,4); hold on; box on; set(gca,'fontsize',sizeF); subplot(3,2,5); hold on; box on; set(gca,'fontsize',sizeF); subplot(3,2,6); hold on; box on; set(gca,'fontsize',sizeF); 
           set(gcf,'Position',get(0,'Screensize'));
        subplot(3,2,1); hh1 = plot(TimeDay_Lsol,ThetaW_LsolC,'Color',colorC,'LineStyle','-','LineWidth',2);
                        hh3 = plot(TimeDay_Lsol,ThetaW_LsolM,'Color',colorM,'LineStyle','--','LineWidth',2);
                        hh2 = plot(TimeDay_Lsol,ThetaW_LsolH,'Color',colorH,'LineStyle',':','LineWidth',2);
                        ylabel('\theta [cm^3w/cm^3]'); xlim([0 TimeT]); ylim([0 0.34]);
                        legend([hh2 hh3 hh1],{'MOD','SEV','CONTROL'},'Orientation','horizontal');
                        xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'A','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');
        subplot(3,2,2); plot(TimeDay_Lsol,CR_LsolC,'Color',colorC,'LineStyle','-','LineWidth',2);
                        plot(TimeDay_Lsol,CR_LsolM,'Color',colorM,'LineStyle','--','LineWidth',2);
                        plot(TimeDay_Lsol,CR_LsolH,'Color',colorH,'LineStyle',':','LineWidth',2);
                        ylabel('CR [mg/cm^3]'); xlim([0 TimeT]); ylim([0 0.38]);  
                        xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'B','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');
        subplot(3,2,3); plot(TimeDay_Lsol,DOC_LsolC.*ThetaW_LsolC,'Color',colorC,'LineStyle','-','LineWidth',2);
                        plot(TimeDay_Lsol,DOC_LsolM.*ThetaW_LsolM,'Color',colorM,'LineStyle','--','LineWidth',2);
                        plot(TimeDay_Lsol,DOC_LsolH.*ThetaW_LsolH,'Color',colorH,'LineStyle',':','LineWidth',2);
                        ylabel('DOC [mg/cm^3]'); xlim([0 TimeT]); ylim([0 0.8]);
                        xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'C','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');
        subplot(3,2,4); plot(TimeDay_Lsol,OSac_LsolC(end,:)*1.4e-3.*AC_LsolC(end,:),'Color',colorC,'LineStyle','-','LineWidth',2);
                        plot(TimeDay_Lsol,OSac_LsolM(end,:)*1.4e-3.*AC_LsolM(end,:),'Color',colorM,'LineStyle','--','LineWidth',2);
                        plot(TimeDay_Lsol,OSac_LsolH(end,:)*1.4e-3.*AC_LsolH(end,:),'Color',colorH,'LineStyle',':','LineWidth',2);
                        ylabel('OS_a_c [mg/cm^3]'); xlim([0 TimeT]); ylim([0 0.4]);
                        xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'D','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');
        subplot(3,2,5); plot(TimeDay_Lsol,AC_LsolC(end,:),'Color',colorC,'LineStyle','-','LineWidth',2);
                        plot(TimeDay_Lsol,AC_LsolM(end,:),'Color',colorM,'LineStyle','--','LineWidth',2);
                        plot(TimeDay_Lsol,AC_LsolH(end,:),'Color',colorH,'LineStyle',':','LineWidth',2);
                        ylabel('AC [mg/cm^3]'); xlim([0 TimeT]); ylim([0 1.8]); xlabel('Time [d]');
                        xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'E','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');
        subplot(3,2,6); plot(TimeDay_Lsol,DC_LsolC(end,:),'Color',colorC,'LineStyle','-','LineWidth',2);
                        plot(TimeDay_Lsol,DC_LsolM(end,:),'Color',colorM,'LineStyle','--','LineWidth',2);
                        plot(TimeDay_Lsol,DC_LsolH(end,:),'Color',colorH,'LineStyle',':','LineWidth',2);
                        ylabel('DC [mg/cm^3]'); xlim([0 TimeT]); ylim([0 1.9]); xlabel('Time [d]');
                        xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'F','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');

figure(4); subplot(3,2,1); hold on; box on; set(gca,'fontsize',sizeF); subplot(3,2,2); hold on; box on; set(gca,'fontsize',sizeF); 
       subplot(3,2,3); hold on; box on; set(gca,'fontsize',sizeF); subplot(3,2,4); hold on; box on; set(gca,'fontsize',sizeF);
       set(gcf,'Position',get(0,'Screensize'));
    subplot(3,2,1); hh3 = plot(TimeDay_Lsol,AC_LsolH(1,:),'Color',colorR,'LineStyle','-','LineWidth',2);
                   hh4 = plot(TimeDay_Lsol,AC_LsolH(2,:),'Color',colorK,'LineStyle','-','LineWidth',2);
                   hh1 = plot(TimeDay_Lsol,AC_LsolH(end,:),'Color',colorH,'LineStyle',':','LineWidth',1.5);
                   hh2 = plot(0,0,'Color',colorM,'LineStyle','--','LineWidth',1.5);
                   ylabel('AC [mg/cm^3]'); xlim([0 TimeT]); ylim([0 1.8]); title('MOD','FontWeight','Normal','FontAngle','Italic');
                   legend([hh1 hh2 hh3 hh4],{'MOD','SEV','FG-strategists','DR-strategists'},'Orientation','horizontal');
                    xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'A','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');
    subplot(3,2,3); plot(TimeDay_Lsol,AC_LsolM(1,:),'Color',colorR,'LineStyle','-','LineWidth',2);
                   plot(TimeDay_Lsol,AC_LsolM(2,:),'Color',colorK,'LineStyle','-','LineWidth',2);
                   plot(TimeDay_Lsol,AC_LsolM(end,:),'Color',colorM,'LineStyle','--','LineWidth',1.5);
                   ylabel('AC [mg/cm^3]'); xlim([0 TimeT]); ylim([0 1.8]); xlabel('Time [d]'); title('SEV','FontWeight','Normal','FontAngle','Italic');
                    xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'C','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');
    subplot(3,2,2); plot(TimeDay_Lsol,DC_LsolH(1,:),'Color',colorR,'LineStyle','-','LineWidth',2);
                   plot(TimeDay_Lsol,DC_LsolH(2,:),'Color',colorK,'LineStyle','-','LineWidth',2);
                   plot(TimeDay_Lsol,DC_LsolH(end,:),'Color',colorH,'LineStyle',':','LineWidth',1.5);
                   ylabel('DC [mg/cm^3]'); xlim([0 TimeT]); ylim([0 1.9]);
                    xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'B','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');
    subplot(3,2,4); plot(TimeDay_Lsol,DC_LsolM(1,:),'Color',colorR,'LineStyle','-','LineWidth',2);
                   plot(TimeDay_Lsol,DC_LsolM(2,:),'Color',colorK,'LineStyle','-','LineWidth',2);
                   plot(TimeDay_Lsol,DC_LsolM(end,:),'Color',colorM,'LineStyle','--','LineWidth',1.5);
                   ylabel('DC [mg/cm^3]'); xlim([0 TimeT]); ylim([0 1.9]); xlabel('Time [d]');
                    xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'D','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');
            
%% % % % % % % % % RW of soils after pretreatment
    
    fprintf('\n \n Computing the response to rewetting...')
    
% % CONTROL
    TimeD = 8.25;   % Length of the second part [d]
    modeRW = 'RW-C';
    Whc0 = ThetaW_LsolC(end)*100/ThetaFC;   % Pretreatment water holding capacity [%]
    FirstRW = 0.25;     % Time of the first RW event [d]
    WhcRW = 50;         % Water holding capacity after rewetting [%]
    Dt_RW = NaN;        % Interval between RW events [d]
    ChT_FC = NaN;       % Characteristic time of drainage (only relevant after rainfall events) [d]
    ChT_0 = NaN;        % Characteristic time of evaporation (only relevant after rainfall events) [d]
    T_hist = TimeDay_Lsol;      % Moisture history of soils (previous treatment): time [d]
    Whc_hist = ThetaW_LsolC*100/ThetaFC;  % Moisture history of soils (previous treatment): WHC [%]
    [POC_C,DOC_C,AC_C,DC_C,EZ_C,OSac_C,OSdc_C,CR_C,Psi_C,ThetaW_C,~,RespAC_C,RespCR_C,RespOS_C,RespDC_C,Growth_C,MassTot_C,MassPart_C,TimeDay_C] = EcoSMMARTS_community(TimeD,DT,modeRW,Whc0,FirstRW,WhcRW,Dt_RW,ChT_FC,ChT_0,T_hist,Whc_hist,...
                                                                                                         ThetaS,ThetaR,a,n,Gamma,POC_LsolC(end),DOC_LsolC(end),EZ_LsolC(end),AC_LsolC(1:end-1,end),DC_LsolC(1:end-1,end),CR_LsolC(end),OSac_LsolC(1:end-1,end),OSdc_LsolC(1:end-1,end),MassTot_LsolC(:,end),MassPart_LsolC(:,end),...
                                                                                                         ADD,Lambda_EE,Y_M,MuCA,MuCZ,NuPOC,Yos,Zre,multD,K_C,K_EE,KdAC,KdACs,...
                                                                                                         KdDC,KdEZ,KdCR,Lambda_r,Lambda_z,Tau_i,Tau_r,Tau_OS,tMemS,tMemB,RelCinlet,vec);
% % SEV
    modeRW = 'RW-C';
    Whc0 = ThetaW_LsolM(end)*100/ThetaFC;   % Pretreatment water holding capacity [%]
    FirstRW = 0.25;     % Time of the first RW event [d]
    WhcRW = 50;     % Water holding capacity after rewetting [%]
    Dt_RW = NaN;    % Interval between RW events [d]
    ChT_FC = NaN;   % Characteristic time of drainage (only relevant after rainfall events) [d]
    ChT_0 = NaN;    % Characteristic time of evaporation (only relevant after rainfall events) [d]
    T_hist = TimeDay_Lsol;      % Moisture history of soils (previous treatment): time [d]
    Whc_hist = ThetaW_LsolM*100/ThetaFC;  % Moisture history of soils (previous treatment): WHC [%]
    [POC_M,DOC_M,AC_M,DC_M,EZ_M,OSac_M,OSdc_M,CR_M,Psi_M,ThetaW_M,~,RespAC_M,RespCR_M,RespOS_M,RespDC_M,Growth_M,MassTot_M,MassPart_M,TimeDay_M] = EcoSMMARTS_community(TimeD,DT,modeRW,Whc0,FirstRW,WhcRW,Dt_RW,ChT_FC,ChT_0,T_hist,Whc_hist,...
                                                                                                         ThetaS,ThetaR,a,n,Gamma,POC_LsolM(end),DOC_LsolM(end),EZ_LsolM(end),AC_LsolM(1:end-1,end),DC_LsolM(1:end-1,end),CR_LsolM(end),OSac_LsolM(1:end-1,end),OSdc_LsolM(1:end-1,end),MassTot_LsolM(:,end),MassPart_LsolM(:,end),...
                                                                                                         ADD,Lambda_EE,Y_M,MuCA,MuCZ,NuPOC,Yos,Zre,multD,K_C,K_EE,KdAC,KdACs,...
                                                                                                         KdDC,KdEZ,KdCR,Lambda_r,Lambda_z,Tau_i,Tau_r,Tau_OS,tMemS,tMemB,RelCinlet,vec);
% % MOD
    modeRW = 'RW-C';
    Whc0 = ThetaW_LsolH(end)*100/ThetaFC;   % Pretreatment water holding capacity [%]
    FirstRW = 0.25;     % Time of the first RW event [d]
    WhcRW = 50;     % Water holding capacity after rewetting [%]
    Dt_RW = NaN;    % Interval between RW events [d]
    ChT_FC = NaN;   % Characteristic time of drainage (only relevant after rainfall events) [d]
    ChT_0 = NaN;    % Characteristic time of evaporation (only relevant after rainfall events) [d]
    T_hist = TimeDay_Lsol;      % Moisture history of soils (previous treatment): time [d]
    Whc_hist = ThetaW_LsolH*100/ThetaFC;    % Moisture history of soils (previous treatment): WHC [%]
    [POC_H,DOC_H,AC_H,DC_H,EZ_H,OSac_H,OSdc_H,CR_H,Psi_H,ThetaW_H,~,RespAC_H,RespCR_H,RespOS_H,RespDC_H,Growth_H,MassTot_H,MassPart_H,TimeDay_H] = EcoSMMARTS_community(TimeD,DT,modeRW,Whc0,FirstRW,WhcRW,Dt_RW,ChT_FC,ChT_0,T_hist,Whc_hist,...
                                                                                                         ThetaS,ThetaR,a,n,Gamma,POC_LsolH(end),DOC_LsolH(end),EZ_LsolH(end),AC_LsolH(1:end-1,end),DC_LsolH(1:end-1,end),CR_LsolH(end),OSac_LsolH(1:end-1,end),OSdc_LsolH(1:end-1,end),MassTot_LsolH(:,end),MassPart_LsolH(:,end),...
                                                                                                         ADD,Lambda_EE,Y_M,MuCA,MuCZ,NuPOC,Yos,Zre,multD,K_C,K_EE,KdAC,KdACs,...
                                                                                                         KdDC,KdEZ,KdCR,Lambda_r,Lambda_z,Tau_i,Tau_r,Tau_OS,tMemS,tMemB,RelCinlet,vec);

% % Sum of multi-compounds
    RespAC_C = [RespAC_C; sum(RespAC_C,1)];
    RespAC_M = [RespAC_M; sum(RespAC_M,1)];
    RespAC_H = [RespAC_H; sum(RespAC_H,1)];
    RespOS_C = [RespOS_C; sum(RespOS_C,1)];
    RespOS_M = [RespOS_M; sum(RespOS_M,1)];
    RespOS_H = [RespOS_H; sum(RespOS_H,1)];
    RespDC_C = [RespDC_C; sum(RespDC_C,1)];
    RespDC_M = [RespDC_M; sum(RespDC_M,1)];
    RespDC_H = [RespDC_H; sum(RespDC_H,1)];
    RespT_C = RespAC_C(end,:) + RespCR_C + RespOS_C(end,:) + RespDC_C(end,:);
    RespT_M = RespAC_M(end,:) + RespCR_M + RespOS_M(end,:) + RespDC_M(end,:);
    RespT_H = RespAC_H(end,:) + RespCR_H + RespOS_H(end,:) + RespDC_H(end,:);
    Growth_C = [Growth_C; sum(Growth_C,1)];
    Growth_M = [Growth_M; sum(Growth_M,1)];
    Growth_H = [Growth_H; sum(Growth_H,1)];

% % Plot results (Rewetting)
figure(5); subplot(2,3,1); hold on; box on; set(gca,'fontsize',sizeF); subplot(2,3,2); hold on; box on; set(gca,'fontsize',sizeF); subplot(2,3,4); hold on; box on; set(gca,'fontsize',sizeF); subplot(2,3,5); hold on; box on; set(gca,'fontsize',sizeF); 
           set(gcf,'Position',get(0,'Screensize'));  
    subplot(2,3,1); hh1 = plot(TimeDay_C-FirstRW,ThetaW_C,'Color',colorC,'LineStyle','-','LineWidth',2);
                    hh3 = plot(TimeDay_C-FirstRW,ThetaW_M,'Color',colorM,'LineStyle','--','LineWidth',2);
                    hh2 = plot(TimeDay_C-FirstRW,ThetaW_H,'Color',colorH,'LineStyle',':','LineWidth',2);
                    ylabel('\theta [cm^3w/cm^3]'); xlim([0 TimeD]-FirstRW); ylim([0 0.34]);
                    legend([hh2 hh3 hh1],{'MOD','SEV','CONTROL'},'Orientation','horizontal');
                    xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'A','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');
    subplot(2,3,5); plot(TimeDay_C-FirstRW,RespT_C,'Color',colorC,'LineStyle','-','LineWidth',2);
                    plot(TimeDay_C-FirstRW,RespT_M,'Color',colorM,'LineStyle','--','LineWidth',2);
                    plot(TimeDay_C-FirstRW,RespT_H,'Color',colorH,'LineStyle',':','LineWidth',2);
                    ylabel('Respiration [mg/cm^3/d]'); xlim([0 TimeD]-FirstRW); ylim([0 0.37]); xlabel('Time [d]'); 
                    xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'C','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');
    subplot(2,3,2); plot(TimeDay_C-FirstRW,Growth_C(end,:)./(Growth_C(end,:)+RespT_C),'Color',colorC,'LineStyle','-','LineWidth',2);
                    plot(TimeDay_C-FirstRW,Growth_M(end,:)./(Growth_M(end,:)+RespT_M),'Color',colorM,'LineStyle','--','LineWidth',2);
                    plot(TimeDay_C-FirstRW,Growth_H(end,:)./(Growth_H(end,:)+RespT_H),'Color',colorH,'LineStyle',':','LineWidth',2);
                    ylabel('CUE [-]'); xlim([0 TimeD]-FirstRW); ylim([0 1]); 
                    xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'B','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');
    subplot(2,3,4); plot(TimeDay_C-FirstRW,Growth_C(end,:),'Color',colorC,'LineStyle','-','LineWidth',2);
                    plot(TimeDay_C-FirstRW,Growth_M(end,:),'Color',colorM,'LineStyle','--','LineWidth',2);
                    plot(TimeDay_C-FirstRW,Growth_H(end,:),'Color',colorH,'LineStyle',':','LineWidth',2);
                    ylabel('Growth [mg/cm^3/d]'); xlim([0 TimeD]-FirstRW); ylim([0 0.25]); xlabel('Time [d]'); 
                    xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'D','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');

figure(6); subplot(3,4,1); hold on; box on; set(gca,'fontsize',sizeF); subplot(3,4,2); hold on; box on; set(gca,'fontsize',sizeF); subplot(3,4,3); hold on; box on; set(gca,'fontsize',sizeF); subplot(3,4,4); hold on; box on; set(gca,'fontsize',sizeF);
           set(gcf,'Position',get(0,'Screensize'));
    subplot(3,4,1); hh1 = plot(TimeDay_C-FirstRW,RespAC_C(end,:),'Color',colorC,'LineStyle','-','LineWidth',2);
                    hh3 = plot(TimeDay_C-FirstRW,RespAC_M(end,:),'Color',colorM,'LineStyle','--','LineWidth',2);
                    hh2 = plot(TimeDay_C-FirstRW,RespAC_H(end,:),'Color',colorH,'LineStyle',':','LineWidth',2);
                    ylabel('Resp. AC [mg/cm^3/d]'); xlim([0 TimeD]-FirstRW); ylim([0 0.21]); xlabel('Time [d]');
                    legend([hh2 hh3 hh1],{'MOD','SEV','CONTROL'},'Orientation','horizontal');
                    xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'A','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');
    subplot(3,4,2); plot(TimeDay_C-FirstRW,RespCR_C(end,:),'Color',colorC,'LineStyle','-','LineWidth',2);
                    plot(TimeDay_C-FirstRW,RespCR_M(end,:),'Color',colorM,'LineStyle','--','LineWidth',2);
                    plot(TimeDay_C-FirstRW,RespCR_H(end,:),'Color',colorH,'LineStyle',':','LineWidth',2);
                    ylabel('Resp. CR [mg/cm^3/d]'); xlim([0 TimeD]-FirstRW); ylim([0 0.21]); xlabel('Time [d]');
                    xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'B','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top'); 
    subplot(3,4,3); plot(TimeDay_C-FirstRW,RespOS_C(end,:),'Color',colorC,'LineStyle','-','LineWidth',2);
                    plot(TimeDay_C-FirstRW,RespOS_M(end,:),'Color',colorM,'LineStyle','--','LineWidth',2);
                    plot(TimeDay_C-FirstRW,RespOS_H(end,:),'Color',colorH,'LineStyle',':','LineWidth',2);
                    ylabel('Resp. OS [mg/cm^3/d]'); xlim([0 TimeD]-FirstRW); ylim([0 0.21]); xlabel('Time [d]'); 
                    xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'C','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');
    subplot(3,4,4); plot(TimeDay_C-FirstRW,RespDC_C(end,:),'Color',colorC,'LineStyle','-','LineWidth',2);
                    plot(TimeDay_C-FirstRW,RespDC_M(end,:),'Color',colorM,'LineStyle','--','LineWidth',2);
                    plot(TimeDay_C-FirstRW,RespDC_H(end,:),'Color',colorH,'LineStyle',':','LineWidth',2);
                    ylabel('Resp. DC [mg/cm^3/d]'); xlim([0 TimeD]-FirstRW); ylim([0 0.21]); xlabel('Time [d]'); 
                    xL=xlim; yL=ylim; text(0.01*xL(2),0.99*yL(2),'D','FontSize', 18,'HorizontalAlignment','left','VerticalAlignment','top');
    
fprintf('\n \n Done!!')
beep

%% % % % % % % % % Model EcoSMMARTS_Community

function [POC,DOC,AC,DC,EZ,OSac,OSdc,CR,Psi,ThetaW,ThetaFC,RespAC,RespCR,RespOS,RespDC,Growth,MassT,MassP,TimeDay] = EcoSMMARTS_community(TimeD,DT,modeRW,Whc0,First_RW,WhcRW,Dt_RW,ChT_FC,ChT_0,T_hist,Whc_hist,...
                                                                                                             ThetaS,ThetaR,a,n,Gamma,POC_in,C_in,EZ_in,AC_in,DC_in,CR_in,OSac_in,OSdc_in,MassTot_in,MassPart_in,...
                                                                                                             ADD,Lambda_EZ,Y_M,Mu_Ca,Mu_Cz,Mu_POC,Yos,Zre,multD,K_C,K_EZ,KdAC,KdACs,...
                                                                                                             KdDC,KdEZ,KdCR,Lambda_r,Lambda_z,Tau_i,Tau_r,Tau_OS,tMemS,tMemB,RelCinlet,vec)      
% % Initialization of variables and some parameters required for the computation
    % Note that all water potentials are defined directly as positive values (as suction) and in cm (1KPa = 10.197 cmH2O)

    TargetMin = TimeD*24*60;                % Temporal variables convertion to min
    t_hist = (T_hist-T_hist(end))*24*60;
    first_RW = First_RW*24*60;
    dt_RW = Dt_RW*24*60;
    chT_FC = ChT_FC*24*60;
    chT_0 = ChT_0*24*60;
    timeRW = [0 first_RW (first_RW+dt_RW):dt_RW:TargetMin 1e99];    % Identification of rewetting events
    Time = 0:DT:TargetMin;      % Time vector [min]
    lv = length(Time);
    
    vecM = max([ones(1,length(vec));vec]); % Pool active (1) or not (0)
    XiWmemB = zeros(vecM(3),1);
    vecK = min([ones(1,length(vec));vec]); % Number of elements in each pool
    ThetaEff = ThetaS-ThetaR;
    PsiFC = 348;             	% Water potential at the field capacity [cm H20] -> Richards and Weaver (1944)
    PsiD = 1e6;                 % Water potential at air-dry [cm H20] -> Modified here to better capture stress at very low potentials
    frThetaD = 3/100;           % WHC at air-dry (fraction of WHC that is the minimum supporting activity)
    PsiMos = 1.5e5;             % Maximum water potential compensated by osmoregulation (no worth at higher): Important for the Y
    wb = 1.4e-3;                % Unit convertion factor in the Van’t Hoff relation [J/cm^4]
    d1 = 1e-4;                  % Approximation to adjust the units in the Van’t Hoff relation [J/cm^4]
    d2 = 60000;                 % Molecular weight of a representative osmolyte [mg/mol] -> Manzoni et al. (2014)
    pi_b = 1022.7;              % Molecular turgor pressure [cm] -> Manzoni et al. (2014)
    
    countRW = 1;
    
    TOL = 1e-35;     % solver tolerance
    opts = optimoptions(@fsolve,'Algorithm','trust-region-dogleg','MaxIter',100000,'MaxFunctionEvaluations',100000,'TolFun',TOL,'TolX',TOL,'Display','off'); % solver options definition
    
    if vecK(1)==0; ADD = 0; end     % redefine some variables according to vecK (activates/deactivates model compartments)
    if vecK(3)==0; KdAC = 0; KdACs = 0; Tau_i = 0; Tau_r = 0; end
    Lambda_EPS = 0; K_EPS = 1; KdEPS = 0; Mu_EPS = 0; EPS_in = 0;  % EPS always inactivated in this model
    if vecK(5)==0; Lambda_EZ = 0; KdEZ = 0; end
    if vecK(6)==0; KdDC = 0; Tau_i = 0; Tau_r = 0; end
    if vecK(7)==0; Tau_OS = 0; KdACs = 0; end
    if vecK(8)==0; Tau_OS = 0; KdACs = 0; end
    if vecK(9)==0; Lambda_z = 0; KdCR = 0; end
    
    POC = zeros(1,lv); POC(:,1) = POC_in; POC_new = POC_in;     % initialize result variables
    DOC = zeros(1,lv);
    EZ = zeros(1,lv);
    AC = zeros(vecM(3),lv); AC(:,1) = AC_in; AC_new = AC(:,1);
    DC = zeros(vecM(6),lv); DC(:,1) = DC_in; DC_new = DC(:,1);
    EPS = zeros(1,lv); EPS(:,1) = EPS_in; EPS_new = EPS_in;
    OSac = zeros(vecM(7),lv); 
    OSdc = zeros(vecM(8),lv);
    CR = zeros(1,lv); CR(:,1) = CR_in; CR_new = CR_in;
    Psi = zeros(1,lv); ThetaW = zeros(1,lv);
    RespAC = zeros(vecM(3),lv); RespCR = zeros(1,lv); RespOS = zeros(vecM(7),lv); RespDC = zeros(vecM(6),lv); Growth = zeros(vecM(3),lv); MassT = zeros(5,lv); MassP = zeros(5,lv); MassPused = zeros(5,lv);

% % Function definitions
    funMon = @(C_,K_) (C_)/(C_+K_);     % Normal and reverse Michaelis-Menten -> Michaelis and Menten (1913)
    funWvgn = @(Psi_) ThetaEff.*(1+(a*Psi_).^n).^(1/n-1) + ThetaR;      % van Genuchten equation -> van Genuchten (1980)
    dTh_dPsi = @(x_) ThetaEff*(1/n-1)*(1+(a*x_)^(n))^(1/n-2) * n*a*(a*x_)^(n-1);    % derivative of van Genuchten equation
    dPsi_dLogPsi = @(x_) x_*log(10);
    funPsiX = @(x_) -(funWvgn(x_)-funWvgn(PsiFC)*frThetaD)/(log10(PsiD)-log10(x_)) - dTh_dPsi(x_)*dPsi_dLogPsi(x_);     % Webb's log-linear expression (Webb, 2000)
    [PsiX,~,~,~]  = fsolve(funPsiX,1000,opts);      % Estimation of the matching point in Webb's log-linear expression
    funW = @(Psi_) max([(dTh_dPsi(PsiX)*dPsi_dLogPsi(PsiX)*log10(Psi_)-log10(PsiD)*dTh_dPsi(PsiX)*dPsi_dLogPsi(PsiX)+funWvgn(PsiFC)*frThetaD).*(Psi_>PsiX) funWvgn(Psi_).*(Psi_<=PsiX)]); % New WRC
    funWcell = @(Cell_) wb*Cell_;       % function to transform mass of cells to volume of cells
    funOS = @(Psi_) d1*d2/298/8.314*(pi_b+Psi_);        % Van't Hoff equation
    OS_Max = funOS(PsiMos);     % Max concentration of OS possible
    funOSeq_rest = @(Psi_) min(funOS(Psi_),OS_Max);     % Concentration of OS needed (takes into account that there is a maximum)
    funTort = @(ThetaW_) (ThetaW_/ThetaS)^Gamma;        % Tortuosity model Hamamoto (2010)
    funAct = @(ThetaW_,ThetaFC_) 1./(ThetaW_-ThetaFC_*frThetaD).*exp((-(log((ThetaW_-ThetaFC_*frThetaD)/(ThetaFC_/2-ThetaFC_*frThetaD))).^2+2*log(ThetaW_-ThetaFC_*frThetaD))/2); % Coefficient of moisture
    funG = @(t_,Mem_) exp((-t_.^2)/2/Mem_^2)/Mem_/sqrt(2*pi);       % Gaussian kernel
    funActMem = @(ThetaW_,ThetaFC_) (ThetaW_<=ThetaFC_/2).*funAct(ThetaW_,ThetaFC_) + (ThetaW_>ThetaFC_/2)*1;       % Coefficient of moisture for the memory function
    funMemW = @(ThetaW_,t_,Mem_,ThetaFC_) sum(funG(t_,Mem_).*funActMem(ThetaW_,ThetaFC_))./sum(funG(t_,Mem_));      % Memory function

% % Initialization of the numerical computation
    ThetaFC = funW(PsiFC);              % Estimation of the field capacity
    w_hist = Whc_hist/100*ThetaFC;      % History of moisture (from WHC to VWC)
    Wc0 = Whc0/100*ThetaFC;
    WcRW = WhcRW/100*ThetaFC;
    Psi0 = fsolve(@(x) funW(x)-Wc0,1000,opts);
    if ~isequal(modeRW,'C')
        PsiRW = fsolve(@(x) funW(x)-WcRW,1000,opts);   
    end
    if timeRW(2) == 0
        if isequal(modeRW,'RW-D')
            funPsi = @(t_) (t_<=chT_FC)*(PsiFC/chT_FC*t_) + (t_>chT_FC)*min((PsiFC*10^((t_-chT_FC)/(chT_0-chT_FC)*log10(PsiRW/PsiFC))),PsiRW);
        elseif isequal(modeRW,'RW-C')
            funPsi = @(t_) PsiRW;
        elseif isequal(modeRW,'D')
            funPsi = @(t_) min((Psi0*10^(t_/chT_0*log10(PsiRW/PsiFC))),PsiRW);
        elseif isequal(modeRW,'C')
            funPsi = @(t_) Psi0;
        end
        countRW = countRW + 1;
    else
        funPsi = @(t_) Psi0;
    end
    Psi_new = funPsi(0);
    Psi(1) = Psi_new;
    ThetaW_new = funW(Psi_new);
    ThetaW(1) = ThetaW_new;
    ThetaAC_new = funWcell(AC_new);
    ThetaDC_new = funWcell(DC_new);
    C_new = C_in*Wc0/ThetaW_new; DOC(:,1) = C_new;
    EZ_new = EZ_in*Wc0/ThetaW_new; EZ(:,1) = EZ_new;
    XiEZ = funMon(EZ_new,K_EZ);
    XiTort = funTort(ThetaW_new);
    XiAct = funAct(ThetaW_new,ThetaFC);
    XiC = funMon(C_new,K_C);
     
    OSeq_new = funOSeq_rest(fsolve(@(Psi_) funW(Psi_) - w_hist(end),1,opts));
    if ~isnumeric(OSac_in) && ~isnumeric(OSdc_in)
        OSac(:,1) = OSeq_new;        OSac_new = OSeq_new;
        OSdc(:,1) = OSeq_new;        OSdc_new = OSeq_new;
    else
        OSac(:,1) = OSac_in;        OSac_new = OSac_in;
        OSdc(:,1) = OSdc_in;        OSdc_new = OSdc_in;
    end
    XiWmemS = funMemW([w_hist ThetaW_new],[t_hist 0],tMemS,ThetaFC);
    for jj=1:vecM(3)
        XiWmemB(jj,1) = funMemW([w_hist ThetaW_new],[t_hist 0],tMemB(jj),ThetaFC); 
    end
    KdACs_old = KdACs.*abs(OSac_new-funOS(Psi_new))/OS_Max;
    
    Y = vecK(7)*Y_M.*(1-min(1,abs(OSac_new-funOS(Psi_new))/OS_Max)).*XiWmemB + (vecK(7)==0)*Y_M.*XiWmemB;
    RespAC(:,1) = (1-Y).*Mu_Ca.*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C).*AC_new;
    RespCR(1) = Mu_Cz.*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C)*CR_new;
    RespOS(:,1) = (1-Yos).*Tau_OS.*ThetaAC_new.*(OSac_new-OSeq_new).*(OSac_new>OSeq_new);
    RespDC(:,1) = Zre*Tau_r.*XiWmemB*funMon(C_new,K_C).*DC_in;
    Growth(:,1) = Y.*Mu_Ca*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C).*AC_new;
    DecompNorm_old = Mu_POC*XiTort*XiEZ*POC_new;
    DecompDist_old = max(0,Mu_POC*multD*(1-XiWmemS)*XiTort*XiEZ*POC_new);
    RecycDec_old = Lambda_r*(sum(KdAC.*AC_new+(1-Lambda_z)*KdACs_old.*AC_new+KdDC.*DC_new)+KdCR*XiAct*CR_new+KdEPS*EPS_new+KdEZ*ThetaW_new*EZ_new) + sum(ThetaAC_new.*OSac_new.*(KdAC+KdACs_old) + ThetaDC_new.*OSdc_new.*KdDC);
    RecycOsm_old = sum(Yos.*Tau_OS.*ThetaAC_new.*(OSac_new-OSeq_new).*(OSac_new>OSeq_new));
    Used_old = sum(Mu_Ca*XiAct*XiC.*AC_new) + Mu_Cz*XiAct*XiC*CR_new + sum(Zre*Tau_r.*XiWmemB*XiC.*DC_new);
    MassTot_new = MassTot_in;
    MassPart_old = MassPart_in;
    MassT(:,1) = MassTot_new;
    MassP(:,1) = MassPart_old;
    MassPused(:,1) = MassPart_old.*Used_old;
    
    ii = 1;
    Time_aux = Time(1);
    TimeNext = Time(2);
    DTX = DT/10;

while TimeNext<=TargetMin
    
    POC_old = POC_new;      % refresh variables (current time step)
    C_old = C_new;
    EZ_old = EZ_new;
    AC_old = AC_new;
    DC_old = DC_new;
    EPS_old = EPS_new;
    OSac_old = OSac_new;
    OSdc_old = OSdc_new;
    OSeq_old = OSeq_new;
    CR_old = CR_new;
    Psi_old = Psi_new;
    ThetaW_old = ThetaW_new;
    ThetaAC_old = ThetaAC_new;
    ThetaDC_old = ThetaDC_new;
    XiTort = funTort(ThetaW_old);
    XiAct = funAct(ThetaW_old,ThetaFC);
    XiC = funMon(C_old,K_C);
    XiEZ = funMon(EZ_old,K_EZ);
    XiEPS = funMon(EPS_old,K_EPS);
    XiWmemS = funMemW([w_hist ThetaW(1:ii) ThetaW_old],[t_hist Time(1:ii) Time_aux] - Time_aux,tMemS,ThetaFC);
    for jj = 1:vecM(3)
        XiWmemB(jj,1) = funMemW([w_hist ThetaW(1:ii) ThetaW_old],[t_hist Time(1:ii) Time_aux] - Time_aux,tMemB(jj),ThetaFC);
    end
    KdACs_old = KdACs.*abs(OSac_old-funOS(Psi_old))/OS_Max;
    Y = vecK(7)*Y_M.*(1-min(1,abs(OSac_old-funOS(Psi_old))/OS_Max)).*XiWmemB + (vecK(7)==0)*Y_M.*XiWmemB;
    
    difTw = (timeRW-Time_aux);
    dt = min([DT TimeNext-Time_aux DTX difTw(difTw>0)]);
    countsol = 0;
    while countsol==0       % Finite differences

        POC_new = POC_old + vecK(1)*dt*(ADD - Mu_POC*(1+multD*(1-XiWmemS))*XiTort*XiEZ*POC_old + (1-Lambda_r)*(sum(KdAC.*AC_old+(1-Lambda_z)*KdACs_old.*AC_old+KdDC.*DC_old)+KdCR*XiAct*CR_old+KdEPS*EPS_old+KdEZ*ThetaW_old*EZ_old));
        C_new = C_old + vecK(2)*dt/ThetaW_old*(Mu_POC*(1+multD*(1-XiWmemS))*XiTort*XiEZ*POC_old - sum(Mu_Ca*XiAct*XiC.*AC_old) - Mu_Cz*XiAct*XiC*CR_old + sum(Yos.*Tau_OS.*ThetaAC_old.*(OSac_old-OSeq_old).*(OSac_old>OSeq_old) - Zre*Tau_r.*XiWmemB*XiC.*DC_old) + Lambda_r*(sum(KdAC.*AC_old+(1-Lambda_z)*KdACs_old.*AC_old+KdDC.*DC_old)+KdCR*XiAct*CR_old+KdEPS*EPS_old+KdEZ*ThetaW_old*EZ_old) + sum(ThetaAC_old.*OSac_old.*(KdAC+KdACs_old) + ThetaDC_old.*OSdc_old.*KdDC));
        AC_new = AC_old + vecK(3)*dt*(Y.*(1-Lambda_EPS*(1-XiWmemB)-Lambda_EZ*(1-XiEZ)).*Mu_Ca*XiAct*XiC.*AC_old - Tau_i*(1-XiAct*XiC).*AC_old + Tau_r.*XiWmemB*XiC.*DC_old - KdAC.*AC_old - KdACs_old.*AC_old - Tau_OS.*ThetaAC_old.*(OSeq_old-OSac_old).*(OSac_old<OSeq_old));
        EPS_new = EPS_old + vecK(4)*dt*(sum(Y*Lambda_EPS.*(1-XiWmemB).*Mu_Ca*XiAct*XiC.*AC_old - Mu_EPS*(1-XiAct*XiC)*XiEPS*XiEZ*AC_old) - KdEPS*EPS_old);
        EZ_new = EZ_old + vecK(5)*dt/ThetaW_old*(sum(Y*Lambda_EZ*(1-XiEZ).*Mu_Ca*XiAct*XiC.*AC_old) - KdEZ*ThetaW_old*EZ_old);
        DC_new = DC_old + vecK(6)*dt*(Tau_i*(1-XiAct*XiC).*AC_old-Tau_r.*XiWmemB*XiC.*DC_old - KdDC.*DC_old);
        OSac_new = OSac_old + vecK(7)*dt./ThetaAC_old.*(-Tau_OS.*ThetaAC_old.*(OSac_old-OSeq_old) - ThetaAC_old.*OSac_old.*Tau_i*(1-XiAct*XiC) + ThetaDC_old.*OSdc_old.*Tau_r.*XiWmemB*XiC - ThetaAC_old.*OSac_old.*(KdAC+KdACs_old));
        OSdc_new = OSdc_old + vecK(8)*dt./ThetaDC_old.*(ThetaAC_old.*OSac_old.*Tau_i*(1-XiAct*XiC) - ThetaDC_old.*OSdc_old.*Tau_r.*XiWmemB*XiC - ThetaDC_old.*OSdc_old.*KdDC);
        CR_new = CR_old + vecK(9)*dt*(sum(Lambda_z*KdACs_old.*AC_old) - KdCR*XiAct*CR_old);
         
        if all([POC_new C_new AC_new' EPS_new EZ_new DC_new' OSac_new' OSdc_new' CR_new]>=0) % just a control to ensure valid solutions
            countsol = 1;
            Time_aux = Time_aux + dt;
            if Time_aux~=timeRW(countRW+1)
                DTX = min(DT,1.01*DTX);
            else
                if countRW == 1
                    if isequal(modeRW,'RW-D')
                        funPsi = @(t_) (t_<=chT_FC)*(PsiFC/chT_FC*t_) + (t_>chT_FC)*min((PsiFC*10^((t_-chT_FC)/(chT_0-chT_FC)*log10(PsiRW/PsiFC))),PsiRW);
                    elseif isequal(modeRW,'RW-C')
                        funPsi = @(t_) PsiRW;
                    elseif isequal(modeRW,'D')
                        funPsi = @(t_) min((Psi0*10^(t_/chT_0*log10(PsiRW/PsiFC))),PsiRW);
                    elseif isequal(modeRW,'C')
                        funPsi = @(t_) Psi0;
                    end
                end
                countRW = countRW + 1;
                DTX = DT/10;
            end
        else
            dt = dt/10;
            DTX = min(DT/10,DTX/10);
            countsol = 0;
        end
    end
    
    XiEZ = funMon(EZ_new,K_EZ);
    XiC = funMon(C_old,K_C);
    DecompNorm_new = Mu_POC*XiTort*XiEZ*POC_new;
    DecompDist_new = max(0,Mu_POC*multD*(1-XiWmemS)*XiTort*XiEZ*POC_new);
    RecycDec_new = Lambda_r*(sum(KdAC.*AC_new+(1-Lambda_z)*KdACs_old.*AC_new+KdDC.*DC_new)+KdCR*XiAct*CR_new+KdEPS*EPS_new+KdEZ*ThetaW_new*EZ_new) + sum(ThetaAC_new.*OSac_new.*(KdAC+KdACs_old) + ThetaDC_new.*OSdc_new.*KdDC);
    RecycOsm_new = sum(Yos.*Tau_OS.*ThetaAC_new.*(OSac_new-OSeq_new).*(OSac_new>OSeq_new));
    Used_new = sum(Mu_Ca*XiAct*XiC.*AC_new) + Mu_Cz*XiAct*XiC*CR_new + sum(Zre*Tau_r.*XiWmemB*XiC.*DC_new);
    MassTot_new = [MassTot_new(1)+(DecompNorm_new+DecompNorm_old)/2*dt MassTot_new(2)+(DecompDist_new+DecompDist_old)/2*dt MassTot_new(3)+(RecycDec_new+RecycDec_old)/2*dt MassTot_new(4)+(RecycOsm_new+RecycOsm_old)/2*dt MassTot_new(5)]'; % Add decomposition and recycling inputs
    MassPart_new = [MassPart_old(1)+(DecompNorm_new+DecompNorm_old)/2*dt MassPart_old(2)+(DecompDist_new+DecompDist_old)/2*dt MassPart_old(3)+(RecycDec_new+RecycDec_old)/2*dt MassPart_old(4)+(RecycOsm_new+RecycOsm_old)/2*dt MassPart_old(5)]'; % Add decomposition and recycling inputs
    MassPart_new = MassPart_new - (Used_new+Used_old)/2*dt*MassPart_new/sum(MassPart_new); % Substract C_use (proportional to pools size)
    
    Psi_new = funPsi(Time_aux-timeRW(countRW));
    ThetaW_new = funW(Psi_new);
    ThetaAC_new = funWcell(AC_new);
    ThetaDC_new = funWcell(DC_new);
    OSeq_new = funOSeq_rest(Psi_new);
        
    if ThetaW_new<ThetaW_old && ThetaW_new<ThetaFC      % Changes of DOC and EZ (concentration increase) due to evaporation 
        C_new = vecK(2)*C_new*min(ThetaW_old,ThetaFC)/ThetaW_new + (vecK(2)==0)*C_new;
        EZ_new = vecK(5)*EZ_new*min(ThetaW_old,ThetaFC)/ThetaW_new + (vecK(5)==0)*EZ_new;
    elseif ThetaW_new>ThetaW_old                        % Change of DOC and EZ concentrations after a moistrue increase
        C_new = vecK(2)*(C_new*ThetaW_old+C_new*RelCinlet*(ThetaW_new-ThetaW_old))/ThetaW_new + (vecK(2)==0)*C_new;
     	EZ_new = vecK(5)*(EZ_new*ThetaW_old)/ThetaW_new + (vecK(5)==0)*EZ_new;
        DTX = DT/10;  
    end
    DecompNorm_old = DecompNorm_new;
    DecompDist_old = DecompDist_new;
    RecycDec_old = RecycDec_new;
    RecycOsm_old = RecycOsm_new;
    Used_old = Used_new;
    MassPart_new = MassPart_new*C_new*ThetaW_new/sum(MassPart_new); % here both drainage (losses) and no drainage (accumulation) are taken into account -> proportional to the current concentration
    MassPart_old = MassPart_new;
    
    if Time_aux >= TimeNext     % storage of output results
        ii = ii+1;
        Time_aux = TimeNext;
        TimeNext = Time(ii) + DT;
        POC(ii) = POC_new;
        DOC(ii) = C_new;
        EZ(ii) = EZ_new;
        AC(:,ii) = AC_new;
        DC(:,ii) = DC_new;
        EPS(ii) = EPS_new;
        OSac(:,ii) = OSac_new;
        OSdc(:,ii) = OSdc_new;
        CR(ii) = CR_new;
        ThetaW(ii) = ThetaW_new;
        Psi(ii) = Psi_new;
        for jj = 1:vecM(3)
            XiWmemB(jj,1) = funMemW([w_hist ThetaW(1:ii)],[t_hist Time(1:ii)] - Time_aux,tMemB(jj),ThetaFC);
        end
        Y = vecK(7)*Y_M.*(1-min(1,abs(OSac_new-funOS(Psi_new))/OS_Max)).*XiWmemB + (vecK(7)==0)*Y_M.*XiWmemB;
        RespAC(:,ii) = (1-Y).*Mu_Ca*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C).*AC_new;
        RespCR(ii) = Mu_Cz*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C)*CR_new;
        RespOS(:,ii) = (1-Yos).*Tau_OS.*ThetaAC_new.*(OSac_new-OSeq_new).*(OSac_new>OSeq_new);
        RespDC(:,ii) = Zre.*Tau_r.*XiWmemB*funMon(C_new,K_C).*DC_new;
        Growth(:,ii) = Y.*Mu_Ca*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C).*AC_new;
        MassT(:,ii) = MassTot_new;
        MassP(:,ii) = MassPart_new;
        MassPused(:,ii) = MassPart_new.*Used_new;
        
    end
end
    
% % Output convertion to the required units
    minTOday = 1/60/24;
    TimeDay = Time*minTOday;
    RespAC = RespAC/minTOday;
    RespCR = RespCR/minTOday;
    RespOS = RespOS/minTOday;
    RespDC = RespDC/minTOday;
    Growth = Growth/minTOday;
    
end
