%% Adiabatic temperature profiles with Tp

%% Constants
N_Avogadro = 6.02214076E23; % Avogadro consant
T0 = 300; % Ambient temperature
GasConstant = 8.3144598; 
Mg2SiO4 = 140.69 * 1E-3;
MgSiO3 =  (200.78 / 2) * 1E-3;

%% PREM data: PREM
% PREM(1,:): depth
% PREM(2,:): pressure
% PREM(3,:): density
PREM = ...
     [0 (3+3)/2 (15+15)/2 (24+24)/2 40 60 (80+80)/2 ... 
        115 150 185 (220+220)/2 265 310 355 400 ... 
        400 450 500 550 600 600 635 660 670 ...
        721 771 771 871 971 1071 1171 1271 1371 ...
        1471 1571 1671 1771 1871 1971 2071 2171 2271 ...
        2371 2471 2571 2671 2741 2741 2771 2871 2891; ...
        0 (0.0299+0.0303)/2 (0.3364+0.3370)/2 (0.6040+0.6043)/2 1.1239 1.7891 (2.4539+2.4546)/2 ...
        3.6183 4.7824 5.9466 (7.1108+7.1115)/2 8.6497 10.2027 11.7702 13.3520 ...
        13.3527 15.2251 17.1311 19.0703 21.0425 21.0426 22.4364 23.4348 23.8342 ...
        26.0783 28.2927 28.2928 32.7623 37.2852 41.8606 46.4882 51.1676 55.8991 ...
        60.6830 65.5202 70.4119 75.3598 80.3660 85.4332 90.5646 95.7641 101.0363 ...
        106.3864 111.8207 117.3465 122.9719 126.9741 126.9742 128.7067 134.5619 135.7509; ...
        982 (982+982)/2 (983+983)/2 (984+984)/2 984 985 (986+986)/2 ...
        987 988 989 (990+990)/2 992 995 995 997 ...
        997 998 999 1000 1000 1000 1001 1001 1001 ...
        1001 1000 1000 998 997 996 994 994 993 ...
        993 993 993 995 996 999 1002 1005 1010 ...
        1016 1023 1031 1041 1048 1048 1052 1065 1068; ...
        1.02 (1.02+2.60)/2 (2.60+2.9)/2 (2.9+3.38076)/2 3.37906 3.37688 (3.37471+3.37471)/2 ...
        3.37091 3.36710 3.36330 (3.35950+3.43578)/2 3.46264 3.48951 3.51639 3.54325 ...
        3.72378 3.78678 3.84980 3.91282 3.97584 3.97584 3.98339 4.27146 4.38071 ... 
        4.41241 4.44316 4.44317 4.50372 4.56307 4.62129 4.67844 4.73460 4.78983 ...
        4.84422 4.89783 4.95073 5.00299 5.05469 5.10590 5.15669 5.20713 5.25729...  
        5.30724 5.35706 5.40681 5.45657 5.49145 5.49145 5.50642 5.55641 5.56645];
PREM(3,:) = PREM(3,:) * 1E-2; % Conversion of the density from g cm^(-3) to kg cm^(-3)
PREM(4,:) = PREM(4,:) * 1E3; % ?????
PREM(5,:) = PREM(3,:) .* PREM(4,:) * 1E-9 * 1E3; % ?????

z_int = 10; % The depth interval for each step of calculation
z_top = [0 410 520 660]; % The depths of the tops of olivine, wadselyite, ringwoodite, and bridgmanite layers 
z_bottom = [410 520 660 2800]; % The depths of the bottoms of olivine, wadselyite, ringwoodite, and bridgmanite layers
DeltaT_trans = [60.0 43.0 -34.0]; % The temperature changes associated with the phase transitions

T_start_0 = 1820; % Temperature at D410
T_step = 1;

%% Preparation

% Fit the P-V-T data to obtain the pressure dependence of the isothermal bulk modulus, Grueneisen parameter, its volume dependence at HP-HT
% Phase number; 
%   olivine = 1, 
%   wadsleyite = 2,
%   ringwoodite = 3,
%   bridgmanite = 4

V00_lat_phase = [289.58 538.14 524.56 162.35]; % The unit-cell volumes in A^3
W_mol_phase = [Mg2SiO4 Mg2SiO4 Mg2SiO4 MgSiO3]; % the molar weights in kg mol-1 
KT0_phase = [127.4 169.2 188 254]; % The isothermal bulk moduli at 0-P in GPa
Theta0_phase = [760 814 830 1030]; % The Debye Temperatures in K at V/V0 = 1 [Isaak et al., 1989]
n_phase = [7 7 7 5]; % The number of atoms in the formula unit
Z_phase = [4 8 8 4]; % The number of unit formulas in the unit cell
V00_mol_phase = N_Avogadro .* V00_lat_phase .* 1E-30 ./ Z_phase ; % The molar volumes in m3 mol-1
rho00_phase = W_mol_phase ./ V00_mol_phase; % the density at ambient P-T in kg m-3

%% Read KT0P, gamma0, and q
    KT0P_phase = [4.21116896988574 4.26046830971026 3.98777150245615 4.08753052699457];
    gamma0_phase = [0.997258584556838 1.12120572353003 1.2839885493668 1.57667292653236];
    q_phase = [2.43854091371459 1.44073720960721 2.42067752768808 1.62939463749366];

cd('E:\Dropbox\BGI\Œ¤‹†\EOS\Geotherm 2020\Software\Geotherm-Tp');
FileName = 'Geotherm2022 Tp.xlsx';
cell0 = {'Depth km'}
xlswrite(FileName, cell0,1,'A1');
%xlswrite(FileName, z_mantle,1,'B1');
    
%% Geotherm calculation

T_minus = -600;
T_plus = 600;
for i_T = T_minus:T_step:T_plus
    T_start = T_start_0 + i_T*T_step;
    
for i_phase = 1 : length(n_phase)
    n = n_phase(i_phase); % the number of atoms in the formula unit of the concerned material
    Z = Z_phase(i_phase); % the number of unit formulas in a unit cell
    V00_mol = V00_mol_phase(i_phase); % the molar volume at ambient P-T in m3 mol-1 (Saxena et al., 1993)
    V00_lat = V00_lat_phase(i_phase); %  the unit cell volume in A^3
    M_mol = W_mol_phase(i_phase); % the molar weight in kg mol-1
    rho00 = rho00_phase(i_phase); % the density at ambient P-T in kg m-3
    KT0 = KT0_phase(i_phase); % The isothermal bulk modulus at standard temperature (fixed)
    Theta0 = Theta0_phase(i_phase); % the Debye Temperature in K at V/V0 = 1 [Isaak et al., 1989]
    
    % Set the temperature at the top of the layer of the region of interest for Monte Carlo Simulation.
    if i_phase == 1 % If the phase is olivine
        deltaT = DeltaT_trans(i_phase) .* 2 ./3;
        T_start_1 = T_start - deltaT; % set the starting temperature with subtracting 2/3 deltaT  
    elseif i_phase == 2 % if the phase is wadsleyite 
        deltaT = DeltaT_trans(i_phase-1) ./3;
        T_start_1 = T_start + deltaT; % set the starting temperature with adding 1/3 deltaT
    else
        T_start_1 = T_region(:, length(z_region)) + DeltaT_trans(i_phase-1); % T_region(:, length(z_region)) is the temperature of the bottom layer in the upper region
    end
    
    z_region = [z_top(i_phase): z_int: z_bottom(i_phase)]; % The depths in the region
    P_region = Earth_DP(PREM, z_region); % The pressures in the region
    dPdz_region = Earth_dPdz(PREM, z_region); % The pressure gradient with respect to the depth in GPa/km
    g_region = Earth_DG(PREM, z_region); % The gravity accelerations in m s-2. The division by 100 is due to the unit conversion from mgal to m s-2

    T_region = zeros(1, length(z_region));
    dTdz_region = zeros(1, length(z_region));
    alpha_region = zeros(1, length(z_region));
    CpW_region = zeros(1, length(z_region));
    CpV_region = zeros(1, length(z_region));
    rho_region = zeros(1, length(z_region));
    VV0_region = zeros(1, length(z_region));
    gamma = zeros(1, length(z_region)); % 
    Cv = zeros(1, length(z_region)); % 
    OneAlphaGammaT = zeros(1, length(z_region)); % 
    KT = zeros(1, length(z_region)); % for
    KS = zeros(1, length(z_region)); % 
    dTdP_region = zeros(1, length(z_region));
    
    KT0P = KT0P_phase(:,i_phase); % Set the replicas of K' in the Monte Carlo simulation
    gamma0 = gamma0_phase(:,i_phase); % Set the replicas of gamma0 in the Monte Carlo simulation
    q = q_phase(:,i_phase); % Set the replicas of q in the Monte Carlo simulation

    if i_phase == 1 
        start_layer = length(z_region);
        end_layer = 1;
        inc_layer = -1;
        z_inc = -z_int;
    else
        start_layer = 1;
        end_layer = length(z_region);
        inc_layer = 1;
        z_inc = +z_int;
    end
    
        T_region( start_layer) = T_start_1; % set the temperature of the first layer
        for  i_layer = start_layer:inc_layer:end_layer % if the calculation has not yet reached the bottom of the regions
            alpha_region( i_layer) = Calc_alpha(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, P_region(i_layer), T_region(i_layer));
%            if ((i_phase == 2) || (i_phase == 3))
%                alpha_region( i_layer) = alpha_region( i_layer) * 0.95;
%            end;
            VV0_region( i_layer) = Inv_MGD_EOS_RJ(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, P_region(i_layer), T_region(i_layer));
            CpW_region(i_layer) = Calc_Cp_wt(i_phase, KT0, KT0P, M_mol, V00_mol, Theta0, gamma0, q, n, P_region(i_layer), T_region(i_layer));
            rho_region( i_layer) = rho00 ./ VV0_region( i_layer);
  
            CpV_region(i_layer) = CpW_region(i_layer) .* rho_region( i_layer); 
            dTdP_region( i_layer) = alpha_region( i_layer) .* T_region( i_layer) ./ CpV_region( i_layer); % added for correction
            dTdz_region( i_layer) = dTdP_region( i_layer) * 1E9 .* dPdz_region(i_layer) * 1E3; 

            gamma( i_layer) =  gamma0_phase(i_phase) * (VV0_region( i_layer) .^ q_phase(i_phase)); % for debug
            OneAlphaGammaT(i_layer) = (1 + alpha_region( i_layer) .* gamma( i_layer) .*  T_region(i_layer)); % for debug
            Cv(i_layer) = CpW_region(i_layer) .* M_mol ./ OneAlphaGammaT(i_layer) ./ (3 * n .* GasConstant) ; % for debug
            KT( i_layer) = Calc_KT(KT0, KT0P_phase(i_phase), V00_mol, Theta0_phase(i_phase), gamma0_phase(i_phase), q_phase(i_phase), n, P_region(i_layer), T_region(i_layer)); % for debug
            KS( i_layer) = KT( i_layer) .* OneAlphaGammaT(i_layer); % for debug
            
%            if (dTdz_region( i_layer) > 1) | (dTdz_region( i_layer) < 0) % If the dT/dz is an abnormal value
%                Adiabat_abnormal= dTdz_region( i_layer) % show the dT/dz
%                alpha_abnormal = alpha_region( i_layer) % show the alpha
%                Cp_abnormal = CpW_region(i_layer) % show the Cp
%                1;
%            end
            if i_layer ~= end_layer
                T_region( i_layer + inc_layer) = T_region( i_layer) + dTdz_region( i_layer) .* z_inc; % calculate the next temperature
            end
        end

    
    if i_phase == 1 % if the stable phase is olivine,
        z_mantle = z_region; % 
        P_mantle = P_region;
        g_mantle = g_region;
        dPdz_mantle = dPdz_region;
        alpha_mantle = alpha_region;
        CpW_mantle = CpW_region;
        CpV_mantle = CpV_region;
        VV0_mantle = VV0_region;
        rho_mantle = rho_region;
        dTdz_mantle = dTdz_region;
        dTdP_mantle = dTdP_region;
        T_mantle = T_region;
        gamma_mantle = gamma; %for debug
        OneAlphaGammaT_mantle = OneAlphaGammaT; %for debug
        Cv_mantle = Cv; %for debug
        KT_mantle = KT; %for debug
        KS_mantle = KS; %for debug
        Tp = T_mantle(6) - (T_mantle(7)-T_mantle(6))*5;

    else % if the stable phase is not olivine but wadsleyite, ringwoodite, or bridgmanite
        z_mantle = [z_mantle z_region]; % add the depth layers of the region to those of the whole mantle
        P_mantle = [P_mantle P_region];
        g_mantle = [g_mantle g_region];
        dTdP_mantle = [dTdP_mantle dTdP_region];
        dPdz_mantle = [dPdz_mantle dPdz_region];
        alpha_mantle = [alpha_mantle alpha_region];
        CpW_mantle = [CpW_mantle CpW_region];
        CpV_mantle = [CpV_mantle CpV_region];
        VV0_mantle = [VV0_mantle VV0_region];
        rho_mantle = [rho_mantle rho_region];
        dTdz_mantle = [dTdz_mantle dTdz_region];
        T_mantle = [T_mantle T_region];       
        gamma_mantle = [gamma_mantle gamma];
        OneAlphaGammaT_mantle = [OneAlphaGammaT_mantle OneAlphaGammaT];
        Cv_mantle = [Cv_mantle Cv];
        KT_mantle = [KT_mantle KT];
        KS_mantle = [KS_mantle KS];     
    end
end;
  xlRange1 = sprintf('A%d',i_T-T_minus+2);
  xlswrite(FileName, Tp, 1, xlRange1);
  xlRange2 = sprintf('B%d',i_T-T_minus+2);
  xlswrite(FileName, T_mantle, 1, xlRange2);
  
end;



%%
% To obtain depth from pressure based on Table _ of PREM [Dziewonski and
% Anderson, PEPI 1984] 
function press = Earth_DP(PREM, depth)

    press0 = depth;
    for i0 = 1 : length(depth)
        i1 = 1;
        while PREM(1, i1) < depth(i0)
            i1 = i1 + 1;
        end
        if i1 > 1
            slope = (PREM(2, i1) - PREM(2, i1-1)) ./ (PREM(1, i1)-PREM(1, i1-1));
            press0(i0) = slope .* (depth(i0)-PREM(1, i1-1)) + PREM(2, i1-1);
        else
            slope = (PREM(2, i1+1) - PREM(2, i1)) ./ (PREM(1, i1+1)-PREM(1, i1));
            press0(i0) = slope .* (depth(i0)-PREM(1, i1)) + PREM(2, i1);
        end
    end
    press = press0;
end    


% To obtain the dP/dz (GPa/km) from depth based on Table _ of PREM [Dziewonski and Anderson, PEPI 1984] 
function dPdz = Earth_dPdz(PREM, depth)
    dPdz = depth;
    for i0 = 1 : length(depth)
        i1 = 1;
        while PREM(1, i1) <= depth(i0)
            i1 = i1 + 1;
        end
        if i1 > 1
        slope = (PREM(5, i1) - PREM(5, i1-1)) ./ (PREM(1, i1)-PREM(1, i1-1));
        dPdz(i0) = slope .* (depth(i0)-PREM(1, i1-1)) + PREM(5, i1-1);
        else
        slope = (PREM(5, i1+1) - PREM(5, i1)) ./ (PREM(1, i1+1)-PREM(1, i1));
        dPdz(i0) = slope .* (depth(i0)-PREM(1, i1)) + PREM(5, i1);
        end
    end
    dPdz = dPdz * 1E-3;
end

% To obtain depth from pressure based on Table _ of PREM [Dziewonski and
% Anderson, PEPI 1984] 
function gravity = Earth_DG(PREM, depth)

    gravity = depth;
    for i0 = 1 : length(depth)
        i1 = 1;
        while PREM(1, i1) < depth(i0)
            i1 = i1 + 1;
        end
        if i1 > 1
        slope = (PREM(3, i1) - PREM(3, i1-1)) ./ (PREM(1, i1)-PREM(1, i1-1));
        gravity(i0) = slope .* (depth(i0)-PREM(1, i1-1)) + PREM(3, i1-1);
        else
        slope = (PREM(3, i1+1) - PREM(3, i1)) ./ (PREM(1, i1+1)-PREM(1, i1));
        gravity(i0) = slope .* (depth(i0)-PREM(1, i1)) + PREM(3, i1);
        end
    end
    gravity = gravity;
end

% Calculation of thermal expansion coefficient using 3BM-MGD EOS at given P and T
function alpha = Calc_alpha(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, P, T)
    dT = ones(1, length(T)) * 10; % temperature increment
%    dT = ones(1, length(T)) * 10; % temperature increment
    
    T_p = T + dT;
    T_m = T - dT;
    
    VV_P0 = Inv_MGD_EOS_RJ(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, P, T);    
    VV0_p = Inv_MGD_EOS_RJ(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, P, T_p);
    VV0_m = Inv_MGD_EOS_RJ(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, P, T_m);
    
%    alpha = (VV0_p - VV0_m) ./ (T_p - T_m);
    alpha = (VV0_p - VV0_m) ./ (T_p - T_m) ./ VV_P0;
end


function VV0 = Inv_MGD_EOS_RJ(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, P, T)
% To obtain volume at pressure and ambient temperature
    % 3rd-order Birch-Murnaghan EOS
    % K: the bulk modulus at zero pressure
    % KP: the pressure derivative of the bulk modulus
    % V: V/V0
    % p: pressure
    
    VV0 = zeros(length(P), length(T));
    for i1 = 1:length(P)
        for i2 = 1:length(T)
            VV01 = 2;
            VV02 = 0.5;
            P1 = MGD_EOS_RJ(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, VV01, T(i2));
            P2 = MGD_EOS_RJ(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, VV02, T(i2));
            while abs(VV01-VV02) > 1E-7
                VV03 = (VV01 + VV02)/2;
                P3 = MGD_EOS_RJ(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, VV03, T(i2));
                if P(i1) > P3
                    VV01 = VV03;
                    P1 = P3;
                else
                    VV02 = VV03;
                    P2 = P3;
                end
            end
            VV0(i1,i2) = (VV01 + VV02)/2;
        end
    end
end


% Mie-Grueneisen-Debye equation of state based on Rigden and Jackson. [PEPI 1996] 
function P_total = MGD_EOS_RJ(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, VV0, T)
    % V0 : the unit cell volume in A^3
    % KT0 : the isothermal bulk modulus at zero pressure in GPa
    % KT0P : the pressure derivative of KT0
    % Theta0 : the Debye Temperature in K at V/V0 = 1
    % gamma0 : the Gr?neisen parameter at zero pressure
    % a : the volume dependence of Gr?neisen parameter 1
    % b : the volume dependence of Gr?neisen parameter 2
    % n : the number of atoms in the formula unit of the concerned material
    % Z : the number of unit formulas in a unit cell
    T0 = 300;
    
    V_mol =  V00_mol .* VV0; % the molar volume at compression
    gamma = gamma0 .*  (VV0.^q); % the Gr?neisen parameter at compression
    ThetaD = Theta0 .* exp(-(gamma-gamma0)/q); % the Debye temperature at compression
    
    Eth_T0 = ThermalEnergy(n, ThetaD, T0); % the thermal energy at standard temperature
    Eth_T = ThermalEnergy(n, ThetaD, T); % the thermal energy at target temperature
    P_T0 = BM3(KT0, KT0P, VV0); % the pressure by compression
    Pth = gamma./V_mol.* (Eth_T-Eth_T0) * 10^(-9); % thermal pressure in GPa
    P_total = Pth + P_T0; % the total pressure

end

% Thermal energy based on the Debye model
% The Debye function is approximated by a 6th-order polynominal
function E_th = ThermalEnergy(n, ThetaD, T)
 % n: the number of atoms in the formula unit of the concerned material
 % ThetaD:  the Debye temperature in K at compression
 % T: temperature in K
    df = [0.33333333 -0.12495878 0.016520035 0.00018490174 ...
        -0.00031034588 0.000035022653 -0.0000013397592]; % The polynominal to approximate the Debye function
    R = 8.314472; % The gas constant in J K-1mol-1 
    x = ThetaD ./ T; % the variable in the Debye function
    E_th = 9 * n .* R .* T .* (df(7)*x.^6 + df(6)*x.^5 + df(5)*x.^4 + df(4)*x.^3 +df(3)*x.^2 + df(2)*x + df(1));
end

function Cp = Calc_Brg_Cp(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, P, T)
    dT = 0.1;
    T_p = T + dT;
    T_m = T - dT;
    alpha = Calc_alpha(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, P, T);
    gamma = Calc_HP_gamma(KT0, KT0P, gamma0, q, P);
    ThetaD = Theta0 .* exp(-(gamma-gamma0)/q); % the Debye temperature at compression
    E_p = ThermalEnergy(n, ThetaD, T_p);
    E_m = ThermalEnergy(n, ThetaD, T_m);
    Cv = (E_p - E_m)/(T_p - T_m);
    OnePlusAlphaGammaT = 1 + alpha .* gamma .* T;
    Cp = OnePlusAlphaGammaT .* Cv;
end

% Calculation of the Grueneisen parameter
function gamma = Calc_HP_gamma(KT0, KT0P, gamma0, q, P)
    VV0_T0 = InvBM3(KT0, KT0P, P);
    gamma = gamma0 .*  (VV0_T0.^q); % the Grueneisen parameter at compression  
end

% The 3rd-order Birch-Murnaghan equation of state
function p = BM3(K, KP, v)
    % 3rd-order Birch-Murnaghan EOS
    % K: the bulk modulus at zero pressure
    % KP: the pressure derivative of the bulk modulus
    % v: V/V0
    % p: pressure
    
    f1 = (3/2)*K;
    f2 = v.^(-7/3)-v.^(-5/3);
    f3 = 1+ (3/4)*(KP-4)*(v.^(-2/3)-1);
    p = f1 .* f2 .*f3;
end

% Calculation of high-pressure isobaric heat capacity per volume at HP
function Cp_HP_wt = Calc_Cp_wt(i_Phase, KT0, KT0P, M_mol, V00_mol, Theta0, gamma0, q, n, P, T)
    if i_Phase >= 1 % If the phase is Brg
        Cp_HP = Calc_Brg_Cp(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, P, T);
    else
        Cp0 = Calc_Cp0(i_Phase, T); % The isobaric heat capacity per mol with T at zero P

        alpha = Calc_alpha(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, P, T);
        alpha0 = Calc_alpha(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, 0, T);
        gamma = Calc_HP_gamma(KT0, KT0P, gamma0, q, P); % The Grueneisen parameter at HP and RT

        OnePlusAlphaGammaT_0 = 1 + alpha0 .* gamma0 .* T; % 1 + alpha*gamma*T factor at RP and HT
        OnePlusAlphaGammaT_HP = 1 + alpha .* gamma .* T; % 1 + alpha*gamma*T factor at HP and HT
        OnePlusAlphaGammaT_ratio = OnePlusAlphaGammaT_HP ./ OnePlusAlphaGammaT_0; % 1 + alpha*gamma*T ratio at HP to RP
        Cp_HP = Cp0 .* OnePlusAlphaGammaT_ratio; % The isoparic heat capacity per mol at HP and HT
    end
    Cp_HP_wt = Cp_HP ./ M_mol; 
end
% Inverse of the 3rd-order Birch Murnaghan equation of state
function V = InvBM3(K, KP, P)
% To obtain volume at pressure and ambient temperature
    % 3rd-order Birch-Murnaghan EOS
    % K: the bulk modulus at zero pressure
    % KP: the pressure derivative of the bulk modulus
    % V: V/V0
    % p: pressure
    
    n_P = length(P);    
    V = P;
    for i0 = 1:n_P
        V1 = 1.1;
        V2 = 0.8;
        P1 = BM3(K, KP, V1);
        P2 = BM3(K, KP, V2);
        while abs(V1-V2) > 1E-6
            V3 = (V1 + V2)/2;
            P3 = BM3(K, KP, V3);
            if P(i0) > P3
                V1 = V3;
                P1 = P3;
            else
                V2 = V3;
                P2 = P3;
            end
        end
        V(i0) = (V1 + V2)/2;
    end
end

function KT = Calc_KT(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, P, T)
    dP =  0.1; % temperature increment
    
    P_p = P + dP;
    P_m = P - dP;
    
    VV0_0 = Inv_MGD_EOS_RJ(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, P, T);    
    VV0_p = Inv_MGD_EOS_RJ(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, P_p, T);
    VV0_m = Inv_MGD_EOS_RJ(KT0, KT0P, V00_mol, Theta0, gamma0, q, n, P_m, T);
    
    KT = - VV0_0 .* (P_p - P_m) ./ (VV0_p - VV0_m);
end


