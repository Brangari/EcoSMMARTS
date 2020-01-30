function [] = Model_EcoSMMARTS_v1()
% % % % Simulation of the EcoSMMARTS model
% % % % Reproduces the respiration data presented in Miller et al. (2005) - Other data can easily be used
% % % % Includes comparison to the models from Lawrence et al. (2009)
% % % % By Albert C Brangarí. Lund University, Centre for Environmental and Climate Research, 2019.

close all;
    
% % Define the temporal frame
    lv = 2000; % Number of elements in vectors (temporal discretization of the results). Used to define DT
    TimeD = 111; % Length of the simulation [d]  !!!!!->You can only run the simulation if moisture is definied until TimeD: TimeD<=LabW4w(end) (in line 17)

% % Data from Miller et al. (2005) or other (you can substitute it by other data)
    LabTimeR = [0 2 4 5 6 7 8 11 13 14 15 16 22 23 24 25 26 29 30 31 38 39 40 43 44 45 51 52 53 54 55 56 57 59 64 66 68 69 70 71 73 80 82 84 85 86 87 92 94 97 98 100 101 106 107 111];
    LabResp4w = [117.26 51.206 43.773 24.527 24.877 23.91 9.542 3.268 13.2 1.42 1.007 0.034 1.96 1 103.53 71.21 69.98 50.4 32.2 15.225 0 0.184 0.51 0.12 0.08 0 0.08 62.02 57.39 65.76 37.89 37.7 34.76 11.84 0.25 0.24 0.3 0 0 0.68 0.12 81.81 41.78 NaN 25.47 10.59 NaN NaN 2.15 0.24 0.07 0.12 0 80.46 27.32 14.96]/1000;
    LabTimeW4w = [0 1.9 4 6 7 23 23.1 25 26 29 31 37.9 51 51.1 53 54 59 63 64.1 79 79.1 82 84 86 87 92 105 105.1 111]; % must be defined at least until TimeD: LabW4w(end)>=TimeD
    LabW4w = [59.9 57.9 41.8 12.3 4.1 4 59.6 57.9 57.7 40.9 12.5 4.1 4 57.5 57.5 53.6 12.5 5 4 4 58.4 58.4 34.1 31 12.3 4.1 4 56.5 25.5];
    T_hist = linspace(0,90,90/TimeD*(lv-1));     % Legacy of the soils (previous treatment): stored during 90 days 
    W_hist = 4*ones(1,length(T_hist));      % Legacy of the soils (previous treatment): here stored under air-dry conditions 4%WHC !!All saturations must be always >= frThetaD (line 269)
% % Parameter definition for Lawrence models
    ka = [0.075 0.05 0.05 0.05]; % [1/d]
    ks = [0.035 0.04 0.05 0.007]; % [1/d]
    kp = 0.00002; % [1/d]
    kd = 0.1; % [1/d]
    kb = 0.35; % [1/d]
    ke = 0.05; % [1/d]
    ep = 0.05; % [1/d]
    Y = 0.35; % [-]
    rm = 0.02; % [1/d]
    sm = 0.4; % [-]
    Mic0 = [0.7 0.35 0.35 0.30]; % [mg/cm^3]
    Slw0 = 3.9; % [mg/cm^3]
    Pas0 = 18.4; % [mg/cm^3]
    Doc0 = [0 0.35 0.34 0.34]; % [mg/cm^3]
    Enz0 = [0 0 0.01 0.01]; % [mg/cm^3] 
    Cbi0 = [0 0 0 0.05]; % [mg/cm^3]

% % Parameter definition for EcoSMMARTS
    % % % % C compartments: POM (particulate organic matter), DOC (dissolved organic carbon), AC (active cells), EPS (extracellular polymeric substances), EZ  (enzymes), DC (dormant cells), OSac (osmolytes in AC), OSdc (osmolytes in DC), ZC (cell residues or "zombie cells")
    POM_in = 15;        % Initial concentration of POM [mg/cm^3] -> Brangarí et al. (2018)
    C_in = 0.05/0.0946; % Initial mass of DOC [mg/cm^3] -> Lawrence et al. (2009) (concentration / VWC)
    EZ_in = 6e-3;       % Initial concentration of EZ [mg/cm^3] -> Brangarí et al. (2018)
    AC_in = 7e-1;       % Initial concentration of AC [mg/cm^3] -> Lawrence et al. (2009)
    DC_in = AC_in;      % Initial concentration of DC [mg/cm^3] -> Blagodatskaya and Kuzyakov (2013)
    ZC_in = 0;          % Initial concentration of ZC [mg/cm^3] -> Assumed
    EPS_in = 0;         % Initial concentration of EPS [mg/cm^3] -> Excluded in this simulation (see Brangarí et al (2018))
    ThetaS = 0.32;      % Saturated water content [-] -> from Miller's data
    ThetaR = 0;         % Residual water content [-] -> from Miller's data
    a = 0.0172;         % Empirical coefficient WRC [1/cm] -> from calibration of Miller's data
    n = 1.3824;         % Empirical coefficient WRC [1/cm] -> from calibration of Miller's data
    Gamma = 1.8;        % Factor of tortuosity [-] -> Hamamoto et al. (2010)
    wb = 1.4e-3;        % Unit convertion factor in the Van’t Hoff relation [J/cm^4]
    a3 = 60000;         % Molecular weight of a representative osmolyte [mg/mol] -> Manzoni et al. (2014)
    
    ADD = 0;            % Inlet [mg/cm^3] -> Assumed 0
    Lambda_EPS = 0;     % Coefficient of carbon allocation diverted towards EPS [-] -> EPS is not included in this simulation (see Brangarí et al (2018) for details)
    Lambda_EZ = 0.01;   % Coefficient of carbon allocation diverted towards enzymes [-] -> Brangarí et al. (2018)
    Y_M = 0.7;          % Yield coefficient of carbon uptake [-]
    Mu_Ca = 3.5e-4;     % Maximum specific uptake rate [1/min] -> From calibration
    Mu_Cz = Mu_Ca/2;    % Maximum specific mineralization rate by cell residues [1/min] -> From calibration
    Mu_POM = 1.7e-5; 	% Maximum specific POM decomposition rate [1/min] -> From calibration
    Mu_EPS = 0;         % Maximum specific EPS decomposition rate [1/min] -> EPS excluded in this simulation (see Brangarí et al (2018))
    Yos = 0.6;         	% Yield coefficient of OS elimination [-] -> Calibrated (if =0: no assimilation -> all to respiration)
    rAC = 1e-2;         % Yield coefficient of reactivation	[-] -> Calibrated
    multDry = 2;        % Coefficient of aggregate disruption [-] -> Calibrated
    K_C = 5;            % Half-saturation constant of DOC [mg/cm^3] (never set to 0) -> Calibrated
    K_EZ = EZ_in/10;    % Half-saturation constant of EZ [mg/cm^3] (never set to 0) -> Assumed
    K_EPS = 1;          % Half-saturation constant of EPS [mg/cm^3] (never set to 0)  -> EPS excluded in this simulation. (see Brangarí et al (2018))
    KdAC = 1e-6;        % Constant rate of decay for AC [1/min] -> Calibrated
    KdACs = KdAC*10;    % Constant rate of decay for AC by stress [1/min] -> Calibrated
    KdDC = KdAC/10;     % Constant rate of decay for DC [1/min] -> Assumed
    KdEZ = KdAC;        % Constant rate of decay for EZ [1/min] -> Assumed
    KdEPS = 0;          % Constant rate of decay for EPS [1/min] -> Excluded in this simulation (see Brangarí et al (2018))
    KdZC = KdAC;        % Constant rate of decay for ZC [1/min] -> Assumed
    Lambda_r = 0.5;     % Coefficient of recycling [-] -> Assumed
    Lambda_z = 0.8;     % Coefficient of cell residues permanence [-] -> Assumed
    Tau_i = 2e-6;   	% Maximum specific inactivation rate [1/min] -> Calibrated
    Tau_a = Tau_i/1.43; % Maximum specific reactivation rate [1/min] -> Konopka (2000)
    Tau_OS = 6e-4;      % Maximum specific rate of osmoregulation [1/min] -> Calibrated
    tMem = [10 3]*24*60;% Coefficients of memory length [min] -> Assumed
    RelCinlet = 0;      % Concentration of DOC and EZ in the inlet water (0: if deonized/rainfall water, 1: water carring compounds (same concentration as in control volume))
    
    vecK = [0 1 1 0 1 1 1 1 1];     % Vector of parameteres that activates (1) or deactivates (0) a specific model compartment [POM DOC AC EPS EZ DC OSac OSdc ZC] (in this simulation POM and EPS are inactive/constant)
   
% % Run simulations (Lawrence)
    Time4w = NaN*zeros(5,lv);
    Resp4w = NaN*zeros(5,lv); 
    Grow4w = NaN*zeros(5,lv);
    wbar = waitbar(0,'Calculating...','Name','LAWRENCE MODELS');
    
    [Time4w(1,:),~,~,~,~,~,~,Resp4w(1,:),Grow4w(1,:)] = LawrenceM1(TimeD,ka(1),ks(1),kp,Y,Mic0(1),Slw0,Pas0,LabTimeW4w,LabW4w,lv);
        waitbar(0.25);
    [Time4w(2,:),~,~,~,~,~,~,Resp4w(2,:),Grow4w(2,:)] = LawrenceM2(TimeD,ka(2),ks(2),kp,kd,rm,sm,Y,Mic0(2),Slw0,Pas0,Doc0(2),LabTimeW4w,LabW4w,lv);
        waitbar(0.5);
    [Time4w(3,:),~,~,~,~,~,~,Resp4w(3,:),Grow4w(3,:)] = LawrenceM3(TimeD,ka(3),ks(3),kp,kd,ke,ep,rm,sm,Y,Mic0(3),Slw0,Pas0,Doc0(3),Enz0(3),LabTimeW4w,LabW4w,lv);
        waitbar(0.75);
    [Time4w(4,:),~,~,~,~,~,~,Resp4w(4,:),Grow4w(4,:)] = LawrenceM4(TimeD,ka(4),ks(4),kp,kd,kb,ke,ep,rm,sm,Y,Mic0(4),Slw0,Pas0,Doc0(4),Enz0(4),Cbi0(4),LabTimeW4w,LabW4w,lv);
        waitbar(1);
        close(wbar)
    
% % Run simulations (EcoSMMARTS)
    [POM,DOC,AC,DC,EPS,EZ,OSac,OSdc,ZC,Psi,ThetaW,RespAC,RespZC,RespOS,RespDC,Growth,TimeDay] = EcoSMMARTS(TimeD,LabTimeW4w,LabW4w,lv,T_hist,W_hist,...
                                                                                                         ThetaS,ThetaR,a,n,Gamma,wb,POM_in,C_in,EZ_in,AC_in,DC_in,ZC_in,EPS_in,a3,ADD,Lambda_EPS,...
                                                                                                         Lambda_EZ,Y_M,Mu_Ca,Mu_Cz,Mu_POM,Mu_EPS,Yos,rAC,multDry,K_C,K_EZ,K_EPS,KdAC,KdACs,...
                                                                                                         KdDC,KdEZ,KdEPS,KdZC,Lambda_r,Lambda_z,Tau_i,Tau_a,Tau_OS,tMem,RelCinlet,vecK);
    Time4w(5,:) = TimeDay;
    Resp4w(5,:) = RespAC+RespZC+RespOS+RespDC;
    Grow4w(5,:) = Growth;                     
    
% % plot solutions
 	figure(1);
        subplot(2,2,1); hold on; box on; set(gca,'fontsize',20); subplot(2,2,2); hold on; box on; set(gca,'fontsize',20); subplot(2,2,3); hold on; box on; set(gca,'fontsize',20); subplot(2,2,4); hold on; box on; set(gca,'fontsize',20); 
        set(gcf,'Position',get(0,'Screensize')); 
        subplot(2,2,1); plot(LabTimeW4w,LabW4w,'k','LineWidth',2);
                        ylabel('WHC [%]');
        subplot(2,2,3); hh(3:6)=plot(Time4w(1:4,:)',Resp4w(1:4,:)','--','LineWidth',1);
                        hh(2)=plot(Time4w(5,:)',Resp4w(5,:)','k','LineWidth',2);
                        hh(1)=plot(LabTimeR,LabResp4w,'o','MarkerEdgeColor',[0 0.5 0],'LineWidth',2);
                        ylabel('Resp. [mg/cm^3/d]');
                        legend(hh,{'data','EcoSMMARTS','FO1','FO2','EC1','EC2'});
        subplot(2,2,4); plot(Time4w(1:4,:)',(Grow4w(1:4,:)./(Resp4w(1:4,:)+Grow4w(1:4,:)))','--','LineWidth',1);
                        plot(Time4w(5,:)',(Grow4w(5,:)./(Resp4w(5,:)+Grow4w(5,:)))','k','LineWidth',2);
                        xlabel('Time [d]'); ylabel('CUE [-]');
        subplot(2,2,2); plot(Time4w(1:4,:)',Grow4w(1:4,:)','--','LineWidth',1);
                        plot(Time4w(5,:)',Grow4w(5,:)','k','LineWidth',2);
                        ylabel('Growth [mg/cm^3/d]');
                         
beep
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Models Lawrence
function [Time,Mic,Slw,Pas,Doc,Enz,Cbi,Resp,Grow] = LawrenceM1(TimeD,ka,ks,kp,Y,Mic0,Slw0,Pas0,TimeW,W,lv)
    
    funW = @(t_) (interp1(TimeW,W,t_)-3.99999)/(60-3.99999); % allows reducing respiration rates under "dry condiitons" (here min assumed to be ~4% WHC)
    
    Time = linspace(0,TimeD,lv);
    dt = TimeD/(lv-1);
    Mic = zeros(1,lv);
    Slw = zeros(1,lv);
    Pas = zeros(1,lv);
    Doc = NaN*zeros(1,lv);
    Enz = NaN*zeros(1,lv);
    Cbi = NaN*zeros(1,lv);
    Resp = zeros(1,lv);
    Grow = zeros(1,lv);
    Mic(1) = Mic0; Slw(1) = Slw0; Pas(1) = Pas0;
    Resp(1) = (1-Y)*ks*Slw(1)*funW(0) + (1-Y)*kp*Pas(1)*funW(0);
    Grow(1) = Y*ks*Slw(1)*funW(0) + Y*kp*Pas(1)*funW(0);

    for ii = 1:lv-1
        Mic(ii+1) = dt*(Y*ks*Slw(ii)*funW(Time(ii)) + Y*kp*Pas(ii)*funW(Time(ii)) - ka*Mic(ii)*funW(Time(ii))) + Mic(ii);
        Slw(ii+1) = dt*(ka*Mic(ii)*funW(Time(ii)) - ks*Slw(ii)*funW(Time(ii))) + Slw(ii);
        Pas(ii+1) = dt*(-kp*Pas(ii)*funW(Time(ii))) + Pas(ii);
        Resp(ii+1) = (1-Y)*ks*Slw(ii+1)*funW(Time(ii+1)) + (1-Y)*kp*Pas(ii+1)*funW(Time(ii+1));
        Grow(ii+1) = Y*ks*Slw(ii+1)*funW(Time(ii+1)) + Y*kp*Pas(ii+1)*funW(Time(ii+1));
    end
end

function [Time,Mic,Slw,Pas,Doc,Enz,Cbi,Resp,Grow] = LawrenceM2(TimeD,ka,ks,kp,kd,rm,sm,Y,Mic0,Slw0,Pas0,Doc0,TimeW,W,lv)
    
    funW = @(t_) (interp1(TimeW,W,t_)-3.99999)/(60-3.99999);
    
    Time = linspace(0,TimeD,lv);
    dt = TimeD/(lv-1);
    Mic = zeros(1,lv);
    Slw = zeros(1,lv);
    Pas = zeros(1,lv);
    Doc = zeros(1,lv);
    Enz = NaN*zeros(1,lv);
    Cbi = NaN*zeros(1,lv);
    Resp = zeros(1,lv);
    Grow = zeros(1,lv);
    Mic(1) = Mic0; Slw(1) = Slw0; Pas(1) = Pas0; Doc(1) = Doc0;
    Resp(1) = (1-Y)*ks*Slw(1)*funW(0) + (1-Y)*kp*Pas(1)*funW(0) + (1-Y)*kd*Doc(1)*funW(0) + rm*Mic(1)*funW(0);
    Grow(1) = Y*ks*Slw(1)*funW(0) + Y*kp*Pas(1)*funW(0) + Y*kd*Doc(1)*funW(0);

    for ii = 1:lv-1
        Mic(ii+1) = dt*(Y*ks*Slw(ii)*funW(Time(ii)) + Y*kp*Pas(ii)*funW(Time(ii)) + Y*kd*Doc(ii)*funW(Time(ii)) - ka*Mic(ii)*funW(Time(ii)) - rm*Mic(ii)*funW(Time(ii))) + Mic(ii);
        Slw(ii+1) = dt*((1-sm)*ka*Mic(ii)*funW(Time(ii)) - ks*Slw(ii)*funW(Time(ii))) + Slw(ii);
        Pas(ii+1) = dt*(-kp*Pas(ii)*funW(Time(ii))) + Pas(ii);
        Doc(ii+1) = dt*(sm*ka*Mic(ii)*funW(Time(ii)) - kd*Doc(ii)*funW(Time(ii))) + Doc(ii);
        Resp(ii+1) = (1-Y)*ks*Slw(ii+1)*funW(Time(ii+1)) + (1-Y)*kp*Pas(ii+1)*funW(Time(ii+1)) + (1-Y)*kd*Doc(ii+1)*funW(Time(ii+1)) + rm*Mic(ii+1)*funW(Time(ii+1));
        Grow(ii+1) = Y*ks*Slw(ii+1)*funW(Time(ii+1)) + Y*kp*Pas(ii+1)*funW(Time(ii+1)) + Y*kd*Doc(ii+1)*funW(Time(ii+1));
    end
end

function [Time,Mic,Slw,Pas,Doc,Enz,Cbi,Resp,Grow] = LawrenceM3(TimeD,ka,ks,kp,kd,ke,ep,rm,sm,Y,Mic0,Slw0,Pas0,Doc0,Enz0,TimeW,W,lv)
    
    funW = @(t_) (interp1(TimeW,W,t_)-3.99999)/(60-3.99999);
    funE = @(E_) E_/(0.001+E_);
    
    Time = linspace(0,TimeD,lv);
    dt = TimeD/(lv-1);
    Mic = zeros(1,lv);
    Slw = zeros(1,lv);
    Pas = zeros(1,lv);
    Doc = zeros(1,lv);
    Enz = zeros(1,lv);
    Cbi = NaN*zeros(1,lv);
    Resp = zeros(1,lv);
    Grow = zeros(1,lv);
    Mic(1) = Mic0; Slw(1) = Slw0; Pas(1) = Pas0; Doc(1) = Doc0; Enz(1) = Enz0;
    Resp(1) = (1-Y)*ks*Slw(1)*funE(Enz0)*funW(0) + (1-Y)*kp*Pas(1)*funW(0) + (1-Y)*kd*Doc(1)*funE(Enz0)*funW(0) + rm*Mic(1)*funW(0) + (1-Y)*ep*Mic(1);
    Grow(1) = Y*ks*Slw(1)*funE(Enz0)*funW(0) + Y*kp*Pas(1)*funW(0) + Y*kd*Doc(1)*funE(Enz0)*funW(0);

    for ii = 1:lv-1
        Mic(ii+1) = dt*(Y*ks*Slw(ii)*funE(Enz(ii))*funW(Time(ii)) + Y*kp*Pas(ii)*funW(Time(ii)) + Y*kd*Doc(ii)*funE(Enz(ii))*funW(Time(ii)) - ka*Mic(ii)*funW(Time(ii)) - rm*Mic(ii)*funW(Time(ii)) - ep*Mic(ii)) + Mic(ii);
        Slw(ii+1) = dt*((1-sm)*ka*Mic(ii)*funW(Time(ii)) - ks*Slw(ii)*funE(Enz(ii))*funW(Time(ii))+ 0.5*ke*Enz(ii)) + Slw(ii);
        Pas(ii+1) = dt*(-kp*Pas(ii)*funW(Time(ii))) + Pas(ii);
        Doc(ii+1) = dt*(sm*ka*Mic(ii)*funW(Time(ii)) - kd*Doc(ii)*funE(Enz(ii))*funW(Time(ii)) + 0.5*ke*Enz(ii)) + Doc(ii);
        Enz(ii+1) = dt*(Y*ep*Mic(ii)-ke*Enz(ii)) + Enz(ii);
        Resp(ii+1) = (1-Y)*ks*Slw(ii+1)*funE(Enz(ii))*funW(Time(ii+1)) + (1-Y)*kp*Pas(ii+1)*funW(Time(ii+1)) + (1-Y)*kd*Doc(ii+1)*funE(Enz(ii))*funW(Time(ii+1)) + rm*Mic(ii+1)*funW(Time(ii+1)) + (1-Y)*ep*Mic(ii+1);
        Grow(ii+1) = Y*ks*Slw(ii+1)*funE(Enz(ii))*funW(Time(ii+1)) + Y*kp*Pas(ii+1)*funW(Time(ii+1)) + Y*kd*Doc(ii+1)*funE(Enz(ii))*funW(Time(ii+1));
    end
end

function [Time,Mic,Slw,Pas,Doc,Enz,Cbi,Resp,Grow] = LawrenceM4(TimeD,ka,ks,kp,kd,kb,ke,ep,rm,sm,Y,Mic0,Slw0,Pas0,Doc0,Enz0,Cbi0,TimeW,W,lv)
    
    funW = @(t_) (interp1(TimeW,W,t_)-3.99999)/(60-3.99999);
    funE = @(E_) E_/(0.001+E_);
    funA = @(A_,B_) A_/(B_/1000/A_+A_);
    
    Time = linspace(0,TimeD,lv);
    dt = TimeD/(lv-1);
    Mic = zeros(1,lv);
    Slw = zeros(1,lv);
    Pas = zeros(1,lv);
    Doc = zeros(1,lv);
    Enz = zeros(1,lv);
    Cbi = zeros(1,lv);
    Resp = zeros(1,lv);
    Grow = zeros(1,lv);
    Mic(1) = Mic0; Slw(1) = Slw0; Pas(1) = Pas0; Doc(1) = Doc0; Enz(1) = Enz0; Cbi(1) = Cbi0;
    Resp(1) = (1-Y)*kb*Cbi0*funA(Mic0,Cbi0)*funW(0) + rm*Mic0*funW(0) + (1-Y)*ep*Mic0;
    Grow(1) = Y*kb*Cbi0*funA(Mic0,Cbi0)*funW(0);

    for ii = 1:lv-1
        Mic(ii+1) = dt*(Y*kb*Cbi(ii)*funA(Mic(ii),Cbi(ii))*funW(Time(ii)) - ka*Mic(ii)*funW(Time(ii)) - rm*Mic(ii)*funW(Time(ii)) - ep*Mic(ii)) + Mic(ii);
        Slw(ii+1) = dt*((1-sm)*ka*Mic(ii)*funW(Time(ii)) - ks*Slw(ii)*funE(Enz(ii))+ 0.5*ke*Enz(ii)) + Slw(ii);
        Pas(ii+1) = dt*(-kp*Pas(ii)) + Pas(ii);
        Doc(ii+1) = dt*(sm*ka*Mic(ii)*funW(Time(ii)) - kd*Doc(ii)*funE(Enz(ii)) + 0.5*ke*Enz(ii)) + Doc(ii);
        Enz(ii+1) = dt*(Y*ep*Mic(ii)-ke*Enz(ii)) + Enz(ii);
        Cbi(ii+1) = dt*(ks*Slw(ii)*funE(Enz(ii)) + kp*Pas(ii) + kd*Doc(ii)*funE(Enz(ii)) - kb*Cbi(ii)*funA(Mic(ii),Cbi(ii))*funW(Time(ii))) + Cbi(ii);
        Resp(ii+1) = (1-Y)*kb*Cbi(ii)*funA(Mic(ii),Cbi(ii))*funW(Time(ii)) + rm*Mic(ii+1)*funW(Time(ii+1)) + (1-Y)*ep*Mic(ii+1);
        Grow(ii+1) = Y*kb*Cbi(ii)*funA(Mic(ii),Cbi(ii))*funW(Time(ii));
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Model EcoSMMARTS

function [POM,DOC,AC,DC,EPS,EZ,OSac,OSdc,ZC,Psi,ThetaW,RespAC,RespZC,RespOS,RespDC,Growth,TimeDay] = EcoSMMARTS(TimeD,TimeW,W,lv,T_hist,W_hist,...
                                                                                                             ThetaS,ThetaR,a,n,Gamma,wb,POM_in,C_in,EZ_in,AC_in,DC_in,ZC_in,EPS_in,a3,ADD,Lambda_EPS,...
                                                                                                             Lambda_EZ,Y_M,Mu_Ca,Mu_Cz,Mu_POM,Mu_EPS,Yos,rAC,multDry,K_C,K_EZ,K_EPS,KdAC,KdACs,...
                                                                                                             KdDC,KdEZ,KdEPS,KdZC,Lambda_r,Lambda_z,Tau_i,Tau_a,Tau_OS,tMem,RelCinlet,vecK)        
% % Initialization of variables and some parameters required for the computation
% Note that all water potentials are defined directly as positive values (as suction) and in cm (1KPa = 10.197 cmH2O)
    wbar = waitbar(0,'Drying and rewetting a virtual soil requires some time...','Name','EcoSMMARTS MODEL');
    wbart = round(linspace(0,lv,100));
    
    Time = linspace(0,TimeD*24*60,lv);  % Time vector [min]
    DT = TimeD/(lv-1)*24*60; % Delta t [min]
    TargetMin = TimeD*24*60;
    ThetaEff = ThetaS-ThetaR;
    tMemS = tMem(1);
    tMemB = tMem(2);
    PsiFC = 348;             	% Water potential at the field capacity [cm H20] -> Richards and Weaver (1944)
    PsiD = 1.5e5;               % Water potential at air-dry [cm H20] -> Manzoni and Katul (2014)
    frThetaD = 3/100;           % WHC at air-dry (fraction of WHC that is the minimum supporting activity)
    PsiMos = 1.5e5;             % Maximum water potential compensated by osmoregulation (no worth at higher): Important for the Y
    t_hist = (T_hist-T_hist(end))*24*60; 
    posRW = diff([W 0])>0;      % identification of rewetting events
    timeRW = TimeW(posRW)*24*60;
    
    TOL = 1e-35;     % solver tolerance
    opts = optimoptions(@fsolve,'Algorithm','trust-region-dogleg','MaxIter',100000,'MaxFunctionEvaluations',100000,'TolFun',TOL,'TolX',TOL,'Display','off'); % solver options definition
    
    if vecK(1)==0; ADD = 0; end     % redefine some variables according to vecK (activates/deactivates model compartments)
    if vecK(3)==0; KdAC = 0; KdACs = 0; Tau_i = 0; Tau_a = 0; end
    if vecK(4)==0; Lambda_EPS = 0; KdEPS = 0; Mu_EPS = 0; end
    if vecK(5)==0; Lambda_EZ = 0; KdEZ = 0; end
    if vecK(6)==0; KdDC = 0; Tau_i = 0; Tau_a = 0; end
    if vecK(7)==0; Tau_OS = 0; KdACs = 0; end
    if vecK(8)==0; Tau_OS = 0; KdACs = 0; end
    if vecK(9)==0; Lambda_z = 0; KdZC = 0; end
    
    POM = zeros(1,lv); POM(:,1) = POM_in; POM_new = POM_in;     % initialize result variables
    DOC = zeros(1,lv); DOC(:,1) = C_in; C_new = C_in;
    EZ = zeros(1,lv); EZ(:,1) = EZ_in; EZ_new = EZ_in;
    AC = zeros(1,lv); AC(:,1) = AC_in; AC_new = AC_in;
    DC = zeros(1,lv); DC(:,1) = DC_in; DC_new = DC_in;
    EPS = zeros(1,lv); EPS(:,1) = EPS_in; EPS_new = EPS_in;
    OSac = zeros(1,lv); OSdc = zeros(1,lv);
    ZC = zeros(1,lv); ZC(:,1) = ZC_in; ZC_new = ZC_in;
    Psi = zeros(1,lv); ThetaW = zeros(1,lv);
    RespAC = zeros(1,lv); RespZC = zeros(1,lv); RespOS = zeros(1,lv); RespDC = zeros(1,lv); Growth = zeros(1,lv);

% % Function definitions
    funMon = @(C_,K_) (C_)/(C_+K_);     % Normal and reverse Michaelis-Menten -> Michaelis and Menten (1913)
    funWvgn = @(Psi_) ThetaEff.*(1+(a*Psi_).^n).^(1/n-1) + ThetaR;      % van Genuchten equation -> van Genuchten (1980)
    dTh_dPsi = @(x_) ThetaEff*(1/n-1)*(1+(a*x_)^(n))^(1/n-2) * n*a*(a*x_)^(n-1);    % derivative of van Genuchten equation
    dPsi_dLogPsi = @(x_) x_*log(10);
    funPsiX = @(x_) -(funWvgn(x_)-funWvgn(PsiFC)*frThetaD)/(log10(PsiD)-log10(x_)) - dTh_dPsi(x_)*dPsi_dLogPsi(x_);     % Webb's log-linear expression (Webb, 2000)
    [PsiX,~,~,~]  = fsolve(funPsiX,1000,opts);      % Estimation of the matching point in Webb's log-linear expression
    funW = @(Psi_) (dTh_dPsi(PsiX)*dPsi_dLogPsi(PsiX)*log10(Psi_)-log10(PsiD)*dTh_dPsi(PsiX)*dPsi_dLogPsi(PsiX)+funWvgn(PsiFC)*frThetaD).*(Psi_>PsiX) + funWvgn(Psi_).*(Psi_<=PsiX); % New WRC
    ThetaFC = funW(PsiFC);      % Estimation of the field capacity

    funWHC = @(t_) interp1(TimeW*24*60,W/100*ThetaFC,t_);       % function to transform WHC from experimental data to volumetric water content
    w_hist = W_hist/100*ThetaFC;      % water pretreatment (from WHC to VWC)
    funWcell = @(Cell_) wb*Cell_;       % function to transform mass of cells to volume of cells
    funOS = @(Psi_) a3*1e-4/298/8.314*(1022.7+Psi_);        % Van't Hoff equation
    OS_Max = funOS(PsiMos);     % Max concentration of OS possible
    funOSeq_rest = @(Psi_) min(funOS(Psi_),OS_Max);     % Concentration of OS needed (takes into account that there is a maximum)
    funTort = @(ThetaW_) (ThetaW_/ThetaS)^Gamma;        % Tortuosity model Hamamoto (2010)
    funAct = @(ThetaW_,ThetaFC_) 1./(ThetaW_-ThetaFC_*frThetaD).*exp((-(log((ThetaW_-ThetaFC_*frThetaD)/(ThetaFC_/2-ThetaFC_*frThetaD))).^2+2*log(ThetaW_-ThetaFC_*frThetaD))/2); % Coefficient of moisture
    funG = @(t_,Mem_) exp((-t_.^2)/2/Mem_^2)/Mem_/sqrt(2*pi);       % Gaussian kernel
    funActMem = @(ThetaW_,ThetaFC_) (ThetaW_<=ThetaFC_/2).*funAct(ThetaW_,ThetaFC_) + (ThetaW_>ThetaFC_/2)*1;       % Coefficient of moisture for the memory function
    funMemW = @(ThetaW_,t_,Mem_,ThetaFC_) sum(funG(t_,Mem_).*funActMem(ThetaW_,ThetaFC_))./sum(funG(t_,Mem_));      % Memory function

% % Initialization of the numerical computation
    ThetaW_new = funWHC(0);
    Psi_new = fsolve(@(Psi_) funW(Psi_) - ThetaW_new,1,opts);
    ThetaW(1) = ThetaW_new;
    ThetaAC_new = funWcell(AC_new);
    ThetaDC_new = funWcell(DC_new);
    Psi(1) = Psi_new;
    OSeq_new = funOSeq_rest(fsolve(@(Psi_) funW(Psi_) - w_hist(end),1,opts));
        OSac(:,1) = OSeq_new;        OSac_new = OSeq_new;
        OSdc(:,1) = OSeq_new;        OSdc_new = OSeq_new;
    XiWmemB = funMemW([w_hist ThetaW_new],[t_hist 0],tMemB,ThetaFC);
    
    Y = vecK(7)*Y_M*(1-min(1,abs(OSac_new-funOS(Psi_new))/OS_Max))*XiWmemB + (vecK(7)==0)*Y_M*XiWmemB;
    RespAC(1) = (1-Y)*Mu_Ca*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C)*AC_new;
    RespZC(1) = Mu_Cz*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C)*ZC_new;
    RespOS(1) = (1-Yos)*Tau_OS*ThetaAC_new*(OSac_new-OSeq_new)*(OSac_new>OSeq_new);
    RespDC(1) = rAC*Tau_a*XiWmemB*funMon(C_new,K_C)*DC_in;
    Growth(1) = Y*Mu_Ca*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C)*AC_new;

    ii = 1;
    Time_aux = Time(1);
    TimeNext = Time(2);
    DTX = DT/10;

while TimeNext<=TargetMin
    
    POM_old = POM_new;      % refresh variables (current time step)
    C_old = C_new;
    EZ_old = EZ_new;
    AC_old = AC_new;
    DC_old = DC_new;
    EPS_old = EPS_new;
    OSac_old = OSac_new;
    OSdc_old = OSdc_new;
    OSeq_old = OSeq_new;
    ZC_old = ZC_new;
    Psi_old = Psi_new;
    ThetaW_old = ThetaW_new;
    ThetaAC_old = ThetaAC_new;
    ThetaDC_old = ThetaDC_new;
    Psi_new = fsolve(@(Psi_) funW(Psi_) - ThetaW_new,Psi_old,opts);
    XiTort = funTort(ThetaW_old);
    XiAct = funAct(ThetaW_old,ThetaFC);
    XiC = funMon(C_old,K_C);
    XiEZ = funMon(EZ_old,K_EZ);
    XiEPS = funMon(EPS_old,K_EPS);
    XiWmemS = funMemW([w_hist ThetaW(1:ii) ThetaW_old],[t_hist Time(1:ii) Time_aux] - Time_aux,tMemS,ThetaFC);
    XiWmemB = funMemW([w_hist ThetaW(1:ii) ThetaW_old],[t_hist Time(1:ii) Time_aux] - Time_aux,tMemB,ThetaFC);
    KdACs_old = KdACs*abs(OSac_old-funOS(Psi_old))/OS_Max;
    Y = vecK(7)*Y_M*(1-min(1,abs(OSac_old-funOS(Psi_old))/OS_Max))*XiWmemB + (vecK(7)==0)*Y_M*XiWmemB;
    
    difTw = (TimeW*24*60-Time_aux);
    dt = min([DT TimeNext-Time_aux DTX difTw(difTw>0)]);
    countsol = 0;
    
    while countsol==0       % Finite differences

        POM_new = POM_old + vecK(1)*dt*(ADD - Mu_POM*(1+multDry*(1-XiWmemS))*XiTort*XiEZ*POM_old + (1-Lambda_r)*(KdAC*AC_old+(1-Lambda_z)*KdACs_old*AC_old+KdDC*DC_old+KdZC*XiAct*ZC_old+KdEPS*EPS_old+KdEZ*ThetaW_old*EZ_old));
        C_new = C_old + vecK(2)*dt/ThetaW_old*(Mu_POM*(1+multDry*(1-XiWmemS))*XiTort*XiEZ*POM_old-Mu_Ca*XiAct*XiC*AC_old-Mu_Cz*XiAct*XiC*ZC_old+Yos*Tau_OS*ThetaAC_old*(OSac_old-OSeq_old)*(OSac_old>OSeq_old)-rAC*Tau_a*XiWmemB*XiC*DC_old+Lambda_r*(KdAC*AC_old+(1-Lambda_z)*KdACs_old*AC_old+KdDC*DC_old+KdZC*XiAct*ZC_old+KdEPS*EPS_old+KdEZ*ThetaW_old*EZ_old)+ThetaAC_old*OSac_old*(KdAC+KdACs_old)+ThetaDC_old*OSdc_old*KdDC);
        AC_new = AC_old + vecK(3)*dt*(Y*(1-Lambda_EPS*(1-XiWmemB)-Lambda_EZ*(1-XiEZ))*Mu_Ca*XiAct*XiC*AC_old+Y*Mu_EPS*(1-XiAct*XiC)*XiEPS*XiEZ*AC_old-Tau_i*(1-XiAct*XiC)*AC_old+Tau_a*XiWmemB*XiC*DC_old-KdAC*AC_old-KdACs_old*AC_old-Tau_OS*ThetaAC_old*(OSeq_old-OSac_old)*(OSac_old<OSeq_old));
        EPS_new = EPS_old + vecK(4)*dt*(Y*Lambda_EPS*(1-XiWmemB)*Mu_Ca*XiAct*XiC*AC_old-KdEPS*EPS_old);  % Inactive here. See Brangarí et al. (2018) for details
        EZ_new = EZ_old + vecK(5)*dt/ThetaW_old*(Y*Lambda_EZ*(1-XiEZ)*Mu_Ca*XiAct*XiC*AC_old-KdEZ*ThetaW_old*EZ_old);
        DC_new = DC_old + vecK(6)*dt*(Tau_i*(1-XiAct*XiC)*AC_old-Tau_a*XiWmemB*XiC*DC_old-KdDC*DC_old);
        OSac_new = OSac_old + vecK(7)*dt/ThetaAC_old*(-Tau_OS*ThetaAC_old*(OSac_old-OSeq_old)-ThetaAC_old*OSac_old*Tau_i*(1-XiAct*XiC)+ThetaDC_old*OSdc_old*Tau_a*XiWmemB*XiC-ThetaAC_old*OSac_old*(KdAC+KdACs_old));
        OSdc_new = OSdc_old + vecK(8)*dt/ThetaDC_old*(ThetaAC_old*OSac_old*Tau_i*(1-XiAct*XiC)-ThetaDC_old*OSdc_old*Tau_a*XiWmemB*XiC-ThetaDC_old*OSdc_old*KdDC);
        ZC_new = ZC_old + vecK(9)*dt*(Lambda_z*KdACs_old*AC_old-KdZC*XiAct*ZC_old);
         
        if all([POM_new C_new AC_new EPS_new EZ_new DC_new OSac_new OSdc_new ZC_new]>=0) % just a control to ensure valid solutions
            countsol = 1;
            Time_aux = Time_aux + dt;
            if Time_aux~=timeRW
                DTX = min(DT,1.01*DTX);
            else
                DTX = DT/10;
            end
        else
            dt = dt/10;
            DTX = min(DT/10,DTX/10);
            countsol = 0;
        end
    end
    
    
    ThetaW_new = funWHC(Time_aux);
    ThetaAC_new = funWcell(AC_new);
    ThetaDC_new = funWcell(DC_new);
    OSeq_new = funOSeq_rest(Psi_new);
    
    if ThetaW_new<ThetaW_old        % Changes of DOC and EZ (concentration increase) due to evaporation
        C_new = vecK(2)*C_new*ThetaW_old/ThetaW_new + (vecK(2)==0)*C_new;
        EZ_new = vecK(5)*EZ_new*ThetaW_old/ThetaW_new + (vecK(5)==0)*EZ_new;
    elseif ThetaW_new>ThetaW_old        % Change of DOC and EZ concentrations after a moistrue increase
        C_new = vecK(2)*(C_new*ThetaW_old+C_new*RelCinlet*(ThetaW_new-ThetaW_old))/ThetaW_new + (vecK(2)==0)*C_new;
     	EZ_new = vecK(5)*(EZ_new*ThetaW_old)/ThetaW_new + (vecK(5)==0)*EZ_new;
        DTX = DT/10;
    end

    if Time_aux >= TimeNext     % storage of output results
        ii = ii+1;
        Time_aux = TimeNext;
        TimeNext = Time(ii) + DT;
        POM(ii) = POM_new;
        DOC(ii) = C_new;
        EZ(ii) = EZ_new;
        AC(ii) = AC_new;
        DC(ii) = DC_new;
        EPS(ii) = EPS_new;
        OSac(ii) = OSac_new;
        OSdc(ii) = OSdc_new;
        ZC(ii) = ZC_new;
        ThetaW(ii) = ThetaW_new;
        Psi(ii) = Psi_new;
        XiWmemB = funMemW([w_hist ThetaW(1:ii)],[t_hist Time(1:ii)] - Time_aux,tMemB,ThetaFC);
        Y = vecK(7)*Y_M*(1-min(1,abs(OSac_new-funOS(Psi_new))/OS_Max))*XiWmemB + (vecK(7)==0)*Y_M*XiWmemB;
        RespAC(ii) = (1-Y)*Mu_Ca*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C)*AC_new;
        RespZC(ii) = Mu_Cz*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C)*ZC_new;
        RespOS(ii) = (1-Yos)*Tau_OS*ThetaAC_new*(OSac_new-OSeq_new)*(OSac_new>OSeq_new);
        RespDC(ii) = rAC*Tau_a*XiWmemB*funMon(C_new,K_C)*DC_new;
        Growth(ii) = Y*Mu_Ca*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C)*AC_new;
        if any(ii==wbart)
            waitbar(ii/lv);
        end
    end
end

% % Output convertion to the required units
    minTOday = 1/60/24;
    
    TimeDay = Time*minTOday;
    RespAC = RespAC/minTOday;
    RespZC = RespZC/minTOday;
    RespOS = RespOS/minTOday;
    RespDC = RespDC/minTOday;
    Growth = Growth/minTOday;
    close(wbar)
end
