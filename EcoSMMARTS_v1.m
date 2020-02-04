%%%%%%%%%%%%%%%%%%%%%%%%%%%% Model EcoSMMARTS_v1

% % % % By Albert C Brangarí. Lund University, Centre for Environmental and Climate Research, 2019.

function [POC,DOC,AC,DC,EZ,OSac,OSdc,CR,Psi,ThetaW,RespAC,RespCR,RespOS,RespDC,Growth,TimeDay] = EcoSMMARTS_v1(TimeD,TimeW,W,lv,T_hist,W_hist,...
                                                                                                             ThetaS,ThetaR,a,n,Gamma,POC_in,C_in,EZ_in,AC_in,DC_in,CR_in,ADD,...
                                                                                                             Lambda_EZ,Y_M,Mu_Ca,Mu_Cz,Mu_POC,Yos,Zre,multD,K_C,K_EZ,KdAC,KdACs,...
                                                                                                             KdDC,KdEZ,KdCR,Lambda_r,Lambda_z,Tau_i,Tau_r,Tau_OS,tMem,RelCinlet,vecK)        
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
    wb = 1.4e-3;                % Unit convertion factor in the Van’t Hoff relation [J/cm^4]
    a3 = 60000;                 % Molecular weight of a representative osmolyte [mg/mol] -> Manzoni et al. (2014)
    
    t_hist = (T_hist-T_hist(end))*24*60; 
    posRW = diff([W 0])>0;      % identification of rewetting events
    timeRW = TimeW(posRW)*24*60;
    
    TOL = 1e-35;     % solver tolerance
    opts = optimoptions(@fsolve,'Algorithm','trust-region-dogleg','MaxIter',100000,'MaxFunctionEvaluations',100000,'TolFun',TOL,'TolX',TOL,'Display','off'); % solver options definition
    
    if vecK(1)==0; ADD = 0; end     % redefine some variables according to vecK (activates/deactivates model compartments)
    if vecK(3)==0; KdAC = 0; KdACs = 0; Tau_i = 0; Tau_r = 0; end
    if vecK(4)==0; Lambda_EZ = 0; KdEZ = 0; end
    if vecK(5)==0; KdDC = 0; Tau_i = 0; Tau_r = 0; end
    if vecK(6)==0; Tau_OS = 0; KdACs = 0; end
    if vecK(7)==0; Tau_OS = 0; KdACs = 0; end
    if vecK(8)==0; Lambda_z = 0; KdCR = 0; end
    
    POC = zeros(1,lv); POC(:,1) = POC_in; POC_new = POC_in;     % initialize result variables
    DOC = zeros(1,lv); DOC(:,1) = C_in; C_new = C_in;
    EZ = zeros(1,lv); EZ(:,1) = EZ_in; EZ_new = EZ_in;
    AC = zeros(1,lv); AC(:,1) = AC_in; AC_new = AC_in;
    DC = zeros(1,lv); DC(:,1) = DC_in; DC_new = DC_in;
    OSac = zeros(1,lv); OSdc = zeros(1,lv);
    CR = zeros(1,lv); CR(:,1) = CR_in; CR_new = CR_in;
    Psi = zeros(1,lv); ThetaW = zeros(1,lv);
    RespAC = zeros(1,lv); RespCR = zeros(1,lv); RespOS = zeros(1,lv); RespDC = zeros(1,lv); Growth = zeros(1,lv);

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
    
    Y = vecK(6)*Y_M*(1-min(1,abs(OSac_new-funOS(Psi_new))/OS_Max))*XiWmemB + (vecK(6)==0)*Y_M*XiWmemB;
    RespAC(1) = (1-Y)*Mu_Ca*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C)*AC_new;
    RespCR(1) = Mu_Cz*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C)*CR_new;
    RespOS(1) = (1-Yos)*Tau_OS*ThetaAC_new*(OSac_new-OSeq_new)*(OSac_new>OSeq_new);
    RespDC(1) = Zre*Tau_r*XiWmemB*funMon(C_new,K_C)*DC_in;
    Growth(1) = Y*Mu_Ca*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C)*AC_new;

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
    OSac_old = OSac_new;
    OSdc_old = OSdc_new;
    OSeq_old = OSeq_new;
    CR_old = CR_new;
    Psi_old = Psi_new;
    ThetaW_old = ThetaW_new;
    ThetaAC_old = ThetaAC_new;
    ThetaDC_old = ThetaDC_new;
    Psi_new = fsolve(@(Psi_) funW(Psi_) - ThetaW_new,Psi_old,opts);
    XiTort = funTort(ThetaW_old);
    XiAct = funAct(ThetaW_old,ThetaFC);
    XiC = funMon(C_old,K_C);
    XiEZ = funMon(EZ_old,K_EZ);
    XiWmemS = funMemW([w_hist ThetaW(1:ii) ThetaW_old],[t_hist Time(1:ii) Time_aux] - Time_aux,tMemS,ThetaFC);
    XiWmemB = funMemW([w_hist ThetaW(1:ii) ThetaW_old],[t_hist Time(1:ii) Time_aux] - Time_aux,tMemB,ThetaFC);
    KdACs_old = KdACs*abs(OSac_old-funOS(Psi_old))/OS_Max;
    Y = vecK(6)*Y_M*(1-min(1,abs(OSac_old-funOS(Psi_old))/OS_Max))*XiWmemB + (vecK(6)==0)*Y_M*XiWmemB;
    
    difTw = (TimeW*24*60-Time_aux);
    dt = min([DT TimeNext-Time_aux DTX difTw(difTw>0)]);
    countsol = 0;
    
    while countsol==0       % Finite differences

        POC_new = POC_old + vecK(1)*dt*(ADD - Mu_POC*(1+multD*(1-XiWmemS))*XiTort*XiEZ*POC_old + (1-Lambda_r)*(KdAC*AC_old+(1-Lambda_z)*KdACs_old*AC_old+KdDC*DC_old+KdCR*XiAct*CR_old+KdEZ*ThetaW_old*EZ_old));
        C_new = C_old + vecK(2)*dt/ThetaW_old*(Mu_POC*(1+multD*(1-XiWmemS))*XiTort*XiEZ*POC_old-Mu_Ca*XiAct*XiC*AC_old-Mu_Cz*XiAct*XiC*CR_old+Yos*Tau_OS*ThetaAC_old*(OSac_old-OSeq_old)*(OSac_old>OSeq_old)-Zre*Tau_r*XiWmemB*XiC*DC_old+Lambda_r*(KdAC*AC_old+(1-Lambda_z)*KdACs_old*AC_old+KdDC*DC_old+KdCR*XiAct*CR_old+KdEZ*ThetaW_old*EZ_old)+ThetaAC_old*OSac_old*(KdAC+KdACs_old)+ThetaDC_old*OSdc_old*KdDC);
        AC_new = AC_old + vecK(3)*dt*(Y*(1-Lambda_EZ*(1-XiEZ))*Mu_Ca*XiAct*XiC*AC_old-Tau_i*(1-XiAct*XiC)*AC_old+Tau_r*XiWmemB*XiC*DC_old-KdAC*AC_old-KdACs_old*AC_old-Tau_OS*ThetaAC_old*(OSeq_old-OSac_old)*(OSac_old<OSeq_old));
        EZ_new = EZ_old + vecK(4)*dt/ThetaW_old*(Y*Lambda_EZ*(1-XiEZ)*Mu_Ca*XiAct*XiC*AC_old-KdEZ*ThetaW_old*EZ_old);
        DC_new = DC_old + vecK(5)*dt*(Tau_i*(1-XiAct*XiC)*AC_old-Tau_r*XiWmemB*XiC*DC_old-KdDC*DC_old);
        OSac_new = OSac_old + vecK(6)*dt/ThetaAC_old*(-Tau_OS*ThetaAC_old*(OSac_old-OSeq_old)-ThetaAC_old*OSac_old*Tau_i*(1-XiAct*XiC)+ThetaDC_old*OSdc_old*Tau_r*XiWmemB*XiC-ThetaAC_old*OSac_old*(KdAC+KdACs_old));
        OSdc_new = OSdc_old + vecK(7)*dt/ThetaDC_old*(ThetaAC_old*OSac_old*Tau_i*(1-XiAct*XiC)-ThetaDC_old*OSdc_old*Tau_r*XiWmemB*XiC-ThetaDC_old*OSdc_old*KdDC);
        CR_new = CR_old + vecK(8)*dt*(Lambda_z*KdACs_old*AC_old-KdCR*XiAct*CR_old);
         
        if all([POC_new C_new AC_new EZ_new DC_new OSac_new OSdc_new CR_new]>=0) % just a control to ensure valid solutions
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
        EZ_new = vecK(4)*EZ_new*ThetaW_old/ThetaW_new + (vecK(4)==0)*EZ_new;
    elseif ThetaW_new>ThetaW_old        % Change of DOC and EZ concentrations after a moistrue increase
        C_new = vecK(2)*(C_new*ThetaW_old+C_new*RelCinlet*(ThetaW_new-ThetaW_old))/ThetaW_new + (vecK(2)==0)*C_new;
     	EZ_new = vecK(4)*(EZ_new*ThetaW_old)/ThetaW_new + (vecK(4)==0)*EZ_new;
        DTX = DT/10;
    end

    if Time_aux >= TimeNext     % storage of output results
        ii = ii+1;
        Time_aux = TimeNext;
        TimeNext = Time(ii) + DT;
        POC(ii) = POC_new;
        DOC(ii) = C_new;
        EZ(ii) = EZ_new;
        AC(ii) = AC_new;
        DC(ii) = DC_new;
        OSac(ii) = OSac_new;
        OSdc(ii) = OSdc_new;
        CR(ii) = CR_new;
        ThetaW(ii) = ThetaW_new;
        Psi(ii) = Psi_new;
        XiWmemB = funMemW([w_hist ThetaW(1:ii)],[t_hist Time(1:ii)] - Time_aux,tMemB,ThetaFC);
        Y = vecK(6)*Y_M*(1-min(1,abs(OSac_new-funOS(Psi_new))/OS_Max))*XiWmemB + (vecK(6)==0)*Y_M*XiWmemB;
        RespAC(ii) = (1-Y)*Mu_Ca*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C)*AC_new;
        RespCR(ii) = Mu_Cz*funAct(ThetaW_new,ThetaFC)*funMon(C_new,K_C)*CR_new;
        RespOS(ii) = (1-Yos)*Tau_OS*ThetaAC_new*(OSac_new-OSeq_new)*(OSac_new>OSeq_new);
        RespDC(ii) = Zre*Tau_r*XiWmemB*funMon(C_new,K_C)*DC_new;
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
    RespCR = RespCR/minTOday;
    RespOS = RespOS/minTOday;
    RespDC = RespDC/minTOday;
    Growth = Growth/minTOday;
    close(wbar)
end
