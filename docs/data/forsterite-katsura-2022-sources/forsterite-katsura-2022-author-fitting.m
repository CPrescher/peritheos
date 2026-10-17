%% Preface
% Examination of the MGD fitting program using Tange et al.'s [2012] data
% with Wolf et al's [2015] parametes
% Different from ver 1, the K = 188 GPa of Rw is taken from Sinogeikin et al. (2003)
% acording to the reviewer's request

% Phase number; 
%   olivine = 1, 
%   wadsleyite = 2, 
%   ringwoodite = 3, 
%   bridgmanite = 4
V00_Phase = [289.58 538.14 524.56 162.373]; % The unit cell volume in A^3
KT0_Phase = [127.4 169.2 188 256.7]; % The isothermal bulk modulus at zero pressure in GPa
KT0P_Phase = [4 4.25817550991571 4 4.044]; % Temporary pressure derivative values of KT0
Theta0_Phase = [760 814 830 1030]; % The Debye Temperature in K at V/V0 = 1 [Isaak et al., 1989]
gamma0_Phase = [1.26 1.13340160150991 1.86 1.8]; % Temperory Grueneisen parameter values ats zero pressure
n_Phase = [7 7 7 5]; % The number of atoms in the formula unit of the concerned material
Z_Phase = [4 8 8 4]; % The number of unit formulas in a unit cell
q_Phase = [2.9 1.62063758409731 1.7 1.9]; % Temporary volume dependence values of Grueneisen parameter (Rigden and Jackson. [PEPI 1996])
T0 = 300; % the standard temperature

PhaseName = {'Olivine', 'Wadsleyite', 'Ringwoodite', 'Bridgmanite'};
for i_phase = 1:length(n_Phase)
    DataFile_Phase{i_phase} = strcat(PhaseName{i_phase}, '_PVT2021.xlsx');
    ResultFile_Phase{i_phase} = strcat(PhaseName{i_phase}, '_MDG2021.xlsx');
end

%% Select mineral
i_phase = 3; % Select phase, olivine = 1, wadsleyite = 2, ringwoodite = 3, bridgmanite = 4
V00 = V00_Phase(i_phase);
KT0 = KT0_Phase(i_phase);
KT0P = KT0P_Phase(i_phase);
Theta0 = Theta0_Phase(i_phase);
gamma0 = gamma0_Phase(i_phase);
n = n_Phase(i_phase);
Z = Z_Phase(i_phase);
q = q_Phase(i_phase);
DataFile = DataFile_Phase{i_phase};
% DataFile = 'Olivine_PVT.xlsx';
nboot = 1000; % the number of bootstrap samples

% prm(1): the pressure derivative of KT0
% prm(2): the Gr?neisen parameter at zero pressure
% prm(3): the temporary volume dependence of Grueneisen parameter
prm_intial = [KT0P gamma0 q];

%% Read data

% Read data
cd 'E:\Dropbox\BGI\Œ¤‹†\EOS\Geotherm 2020\Software\PVT2020-MGD 2021';
Data = xlsread(DataFile);
P = Data(:,1); % pressure
sP = Data(:,2); % error of pressure
VV0 = Data(:,3); % volume, V/V0
sVV0 = Data(:,4); % error of volume
T = Data(:,5); % temperature
sT = Data(:,6); % the errors of temperature
n_data = length(P); % the number of data points

%% Produce nboot data set replicas by the bootstrap method

% Create a uniformly-distributed random data index set of HT data
% The range of indexes is from 1 to the number of data points "n_data"

rand_index = randi(n_data, nboot, n_data); % The "nboot" sets
rand_normal_P = randn(nboot, n_data); % Create a normally-distributed random numbers for "P"
rand_normal_VV0 = randn(nboot, n_data); % Create a normally-distributed random numbers for "V"
rand_normal_T = randn(nboot, n_data); % Create a normally-distributed random numbers for "T"
%rand_normal_P = zeros(nboot, n_data); % For debugging
%rand_normal_VV0 = zeros(nboot, n_data); % For debugging
prm_boot = zeros(nboot, 3);% produce a matrix of nboot sets of 3 parameters
P_boot = zeros(nboot, n_data);
VV0_boot = zeros(nboot, n_data);
T_boot = zeros(nboot, n_data);

for i0 = 1:nboot
    % create a bootstrap P-V-T data set
%    if (i0 > 11)
%        i0
%    end
    
    for i1 = 1 : n_data 
        data_index = rand_index(i0, i1);
        % Create bootstrap data of pressure with an normally distributed error
%        P_index = P(data_index); % For debug
%        NormalDistrP = rand_normal_P(i0, i1); % For debug
%        sP_index = sP(data_index); % For debug
%        sP_temp = NormalDistrP * sP_index; % For debug
%        P_temp = P_index + sP_temp; % For debug
        P_boot(i0, i1) = P(data_index) + rand_normal_P(i0, i1) .* sP(data_index); 
        % Create bootstrap data of volume with an normally distributed error
%        VV0_index = VV0(data_index); % For debug
%        NormalDistrVV0 = rand_normal_VV0(i0, i1); % For debug
%        sVV0_index = sVV0(data_index); % For debug
%        sVV0_temp = NormalDistrVV0 * sVV0_index; % For debug
%        VV0_temp = VV0_index + sVV0_temp; % For debug
        VV0_boot(i0, i1) = VV0(data_index) + rand_normal_VV0(i0, i1) .* sVV0(data_index); 
        % Create bootstrap data of temperature without an error
        T_boot(i0, i1) = T(data_index) + rand_normal_T(i0, i1) .* sT(data_index);
    end
end

%% Evaluation of KT0', gamma, q
for i0 = 1:nboot
%    if (i0 > 11) 
%    i0
%    end
    
    prm1 = NlLSq(prm_intial, KT0, V00, Z, Theta0, n, T0, P_boot(i0,:), VV0_boot(i0,:), T_boot(i0,:));
    prm_boot(i0,:) = prm1;
    P_fit = MGD_EOS_RJ(KT0, prm1(1), V00, Z, Theta0, prm1(2), prm1(3), n, T0, VV0_boot(i0,:), T_boot(i0,:));
%    h = figure; % Plot for debugging
%    plot(VV0_boot(i0,:), P_boot(i0,:), 'b o', VV0_boot(i0,:), P_fit, 'r o'); % Plot for debugging
%    xlabel('\itV\rm/\itV\rm_0'); % Plot for debugging
%    ylabel('\itP\rm (GPa)'); % Plot for debugging
%    title('Fitting of olivine data to the 3BM & MGD EOS'); % Plot for debugging
%    legend('Data', 'Fit'); % Plot for debugging
%    close(h); % Plot for debugging
end

% Examination of Bootstrap
    h = figure; 
    plot(prm_boot(:,2), prm_boot(:,3), 'b o'); % gamma0-q
    ylabel('\itq\rm'); % Plot for debugging
    xlabel('\itƒÁ\rm_0'); % Plot for debugging
    str_title = strcat(PhaseName{i_phase}, ' data fitted to the 3BM & MGD EOS');
    title(str_title); % Plot for debugging
    close(h); % Plot for debugging

prm_ave = mean(prm_boot);
prm_sd = std(prm_boot);

P_fit = MGD_EOS_RJ(KT0, prm_ave(1), V00, Z, Theta0, prm_ave(2), prm_ave(3), n, T0, VV0, T);
h = figure;
plot(VV0, P, 'b o', VV0, P_fit, 'r o');
xlabel('\itV\rm/\itV\rm_0');
ylabel('\itP\rm (GPa)');
title('Fitting of olivine data to the 3BM & MGD EOS');
legend('Data', 'Fit');

%% Write a result Excel table

ResultFileName = ResultFile_Phase{i_phase};
delete(ResultFileName);

SheetName01 = 'Parameters';
cell0 = {'KT0 (fixed)' ; 'KT_prime' ; 'V00 (fixed)' ; 'n (fixed)'; 'Z (fixed)' ; 'Theta0 (fixed)' ; 'gamma0' ; 'q' };
cell1 = {KT0; prm_ave(1); V00; n; Z; Theta0; prm_ave(2); prm_ave(3)};
cell2 = {''; prm_sd(1,1); ''; ''; ''; ''; prm_sd(1,2); prm_sd(1,3)};
[ ~, ~] = xlswrite(ResultFileName, cell0, SheetName01,'A1');
[ ~, ~] = xlswrite(ResultFileName, cell1, SheetName01,'B1');
[ ~, ~] = xlswrite(ResultFileName, cell2, SheetName01,'C1');

SheetName02 = 'Data fit';
cell3 = {'P (GPa)', 'ƒÐP (GPa)', 'V/V0', 'ƒÐV/V0', 'T (K)', 'P_fit (GPa)', 'dif_P (GPa)', 'ThetaD (K)', 'gamma'};
[ ~, ~] = xlswrite(ResultFileName, cell3, SheetName02,'A1');
[ ~, ~] = xlswrite(ResultFileName, P, SheetName02,'A2');
[ ~, ~] = xlswrite(ResultFileName, sP, SheetName02,'B2');
[ ~, ~] = xlswrite(ResultFileName, VV0, SheetName02,'C2');
[ ~, ~] = xlswrite(ResultFileName, sVV0, SheetName02,'D2');
[ ~, ~] = xlswrite(ResultFileName, T, SheetName02,'E2');
[ ~, ~] = xlswrite(ResultFileName, P_fit, SheetName02,'F2');
diff_P = P - P_fit;
[ ~, ~] = xlswrite(ResultFileName, diff_P, SheetName02,'G2');
gamma = prm_ave(1,2) .*  (VV0.^prm_ave(1,3)); % the Gr?neisen parameter at compression
ThetaD = Theta0 .* exp(-(gamma - prm_ave(1,2))/prm_ave(1,3)); % the Debye temperature at compression
[ ~, ~] = xlswrite(ResultFileName, ThetaD, SheetName02,'H2');
[ ~, ~] = xlswrite(ResultFileName, gamma, SheetName02,'I2');

SheetName03 = 'V_V0';
cell4 = {'P (GPa) \ T (K)'};
P_array = [0:1:120]';
T_array = [T0:100:3000];
VV0_table = Inv_MGD_EOS_RJ(KT0, prm_ave(1,1), V00, Z, Theta0, prm_ave(1,2), prm_ave(1,3), n, P_array, T0, T_array);
[ ~, ~] = xlswrite(ResultFileName, cell4, SheetName03,'A1');
[ ~, ~] = xlswrite(ResultFileName, T_array, SheetName03,'B1');
[ ~, ~] = xlswrite(ResultFileName, P_array, SheetName03,'A2');
[ ~, ~] = xlswrite(ResultFileName, VV0_table, SheetName03,'B2');

SheetName04 = 'V_VT0';
VVT0_table = VV0_table;
for i1 = 1:length(P_array);
    VVT0_table(i1,:) = VV0_table(i1,:) ./ VV0_table(i1,1);
end
[ ~, ~] = xlswrite(ResultFileName, cell4, SheetName04,'A1');
[ ~, ~] = xlswrite(ResultFileName, T_array, SheetName04,'B1');
[ ~, ~] = xlswrite(ResultFileName, P_array, SheetName04,'A2');
[ ~, ~] = xlswrite(ResultFileName, VVT0_table, SheetName04,'B2');

SheetName05 = 'alpha';
dT = ones(1, length(T_array));
dT = dT * 10;
Tp_array = T_array + dT;
Tm_array = T_array - dT;
VV0p_table = Inv_MGD_EOS_RJ(KT0, prm_ave(1,1), V00, Z, Theta0, prm_ave(1,2), prm_ave(1,3), n, P_array, T0, Tp_array);
VV0m_table = Inv_MGD_EOS_RJ(KT0, prm_ave(1,1), V00, Z, Theta0, prm_ave(1,2), prm_ave(1,3), n, P_array, T0, Tm_array);
% alpha_table = (VV0p_table - VV0m_table) ./ dT / 2;
VV0_average = (VV0p_table + VV0m_table) / 2;
alpha_table = (VV0p_table - VV0m_table) ./ VV0_average ./ dT / 2;
[ ~, ~] = xlswrite(ResultFileName, cell4, SheetName05,'A1');
[ ~, ~] = xlswrite(ResultFileName, T_array, SheetName05,'B1');
[ ~, ~] = xlswrite(ResultFileName, P_array, SheetName05,'A2');
[ ~, ~] = xlswrite(ResultFileName, alpha_table, SheetName05,'B2');

SheetName06 = 'Prm_boot';
cell0 = {'KT_prime' , 'gamma0' , 'q' };
% cell1 = {KT0, prm_boot(:,1), V00, n, Z; Theta0, prm_boot(:,2), prm(:,3)};
[ ~, ~] = xlswrite(ResultFileName, cell0, SheetName06,'A1');
[ ~, ~] = xlswrite(ResultFileName, prm_boot, SheetName06,'A2');

1;
close(h);

%% Subroutine

% The parameter derivatives of high-pressure thermal expansion
function dP_dPrm_res = dP_dPrm(prm, KT0, V_lat0, Z, Theta0, n, T0, VV0, T)
    % prm(1): thermal expansion coefficient at standard temperature and zero
    %   pressure multiplied by 1E+5
    % prm(2): temperature derivative of thermal expansion coefficient at zero
    %   pressure multiplied by 1E+8
    % prm(3): Anderson-Grueneisen parameter 
    % VT0: the zero-temperature volume at a given pressure
    % T0: the standard temperature
    % T: the temperature
    n_prm = length(prm);
    prm_p = zeros(n_prm);
    prm_m = zeros(n_prm);
    n_data = length(T);
    P_p = zeros(n_prm, n_data);
    P_m = zeros(n_prm, n_data);
    dP_dPrm_res = zeros(n_prm, n_data);
    
    for i0 = 1 : n_prm
        prm_p = prm; % to copy prm to prm_p
        prm_m = prm; % to copy prm to prm_m
        prm_p(i0) = prm(i0) * 1.001; % the i0-th component of prm increased by 0.1%
        prm_m(i0) = prm(i0) * 0.999; % the i0-th component of prm decreased by 0.1%
        P_p(i0,:) = MGD_EOS_RJ(KT0, prm_p(1), V_lat0, Z, Theta0, prm_p(2), prm_p(3), n, T0, VV0, T); % P by the increased prm
        P_m(i0,:) = MGD_EOS_RJ(KT0, prm_m(1), V_lat0, Z, Theta0, prm_m(2), prm_m(3), n, T0, VV0, T); % P by the decreased prm
        dP_dPrm_res(i0,:) = (P_p(i0,:) - P_m(i0,:)) / (prm_p(i0) - prm_m(i0)); % differente P with respect to the n-th pareameter
    end
end


function beta = Beta_Array(prm, KT0, V_lat0, Z, Theta0, n, T0, P, VV0, T)
    n_prm = length(prm);  % the number of parameters    
    beta0 =  dP_dPrm(prm, KT0, V_lat0, Z, Theta0, n, T0, VV0, T);
    P_sim =  MGD_EOS_RJ(KT0, prm(1), V_lat0, Z, Theta0, prm(2), prm(3), n, T0, VV0, T);
    diff_P = P - P_sim;
    beta_sum = zeros(1, n_prm);
    for i1 = 1 : n_prm
        beta1 = diff_P .* beta0(i1, :);
        beta_sum(i1) = sum(beta1); 
    end
    beta = beta_sum;
end

% Hesse matrix according to Seto's home page
function hesse = HesseMatrix(prm, lambda, KT0, V_lat0, Z, Theta0, n, T0, P, VV0, T)
    % prm(1): thermal expansion coefficient at standard temperature and zero
    %   pressure multiplied by 1E+5
    % prm(2): temperature derivative of thermal expansion coefficient at zero
    %   pressure multiplied by 1E+8
    % prm(3): Anderson-Grueneisen parameter 
    % VT0: the zero-temperature volume at a given pressure
    % T0: the standard temperature
    % T: the temperature

    n_prm = length(prm); % the number of parameters   
    dP_dPrm1 = dP_dPrm(prm, KT0, V_lat0, Z, Theta0, n, T0, VV0, T); % the differentiation of P with respect to the parameters 
    alpha = zeros(n_prm, n_prm); % The components of the error matrix
    for i1 = 1 : n_prm
        for i2 = 1 : n_prm
            if i2>=i1
                alpha(i1,i2) = sum(dP_dPrm1(i1,:) .* dP_dPrm1(i2,:)); % the component of the error matrix is the sum of the product of the parameter derivatives
                if i1 == i2
                    alpha(i1,i2) = (1 + lambda) .* alpha(i1,i2); % the conversion parameter ƒÉ is included in the diagonal components
                else
                    alpha(i2,i1) = alpha(i1, i2);
                end
            end
        end
    end
hesse = alpha;
end

function chisq = ChisqPressure(prm, KT0, V_lat0, Z, Theta0, n, T0, P_dat, VV0, T)
  P_sim = MGD_EOS_RJ(KT0, prm(1), V_lat0, Z, Theta0, prm(2), prm(3), n, T0, VV0, T);
  delta_P_sq = (P_sim - P_dat) .^ 2;
  chisq = sum(delta_P_sq);
end

% Non-lnear least square fitting by the Marquardt method 
function result = NlLSq(prm_initial, KT0, V_lat0, Z, Theta0, n, T0, P, VV0, T)
    % prm0: initial values of parameters

    lambda = 0.01;
    mag_prm_inc = 0.01;
    prm0 = prm_initial;
    chisq0 = ChisqPressure(prm0, KT0, V_lat0, Z, Theta0, n, T0, P, VV0, T);
    while mag_prm_inc > 1E-4
        Hesse = HesseMatrix(prm0, lambda, KT0, V_lat0, Z, Theta0, n, T0, P, VV0, T);
        InvHesse = inv(Hesse); 
        beta = Beta_Array(prm0, KT0, V_lat0, Z, Theta0, n, T0, P, VV0, T);
        prm_inc = InvHesse * beta'; 
        prm1 = prm0 + prm_inc';
        chisq1 = ChisqPressure(prm1, KT0, V_lat0, Z, Theta0, n, T0, P, VV0, T);
        diff_chisq = chisq1 - chisq0;
        if diff_chisq < 0
            prm0 = prm1;
            chisq0 = chisq1;
            lambda = lambda ./ 2.0;
        else
            lambda = lambda .* 2.0;
        end
        prm_inc_norm = prm_inc' ./ prm0;
        mag_prm_inc = norm(prm_inc_norm);
    end
    result = prm0;
end

% Mie-Grueneisen-Debye equation of state based on Rigden and Jackson. [PEPI 1996] 
function P_total = MGD_EOS_RJ(KT0, KT0P, V_lat0, Z, Theta0, gamma0, q, n, T0, VV0, T)
    % V0 : the unit cell volume in A^3
    % KT0 : the isothermal bulk modulus at zero pressure in GPa
    % KT0P : the pressure derivative of KT0
    % Theta0 : the Debye Temperature in K at V/V0 = 1
    % gamma0 : the Gr?neisen parameter at zero pressure
    % a : the volume dependence of Gr?neisen parameter 1
    % b : the volume dependence of Gr?neisen parameter 2
    % n : the number of atoms in the formula unit of the concerned material
    % Z : the number of unit formulas in a unit cell

    NA = 6.0221367E+23; % Avogadro number

    V_lat = VV0 .* V_lat0; % the lattice volume at compression
    V_mol = (NA ./ Z) .* V_lat * 10^(-30); % the molar volume at compression
    gamma = gamma0 .*  (VV0.^q); % the Gr?neisen parameter at compression
    ThetaD = Theta0 .* exp(-(gamma-gamma0)/q); % the Debye temperature at compression
    
    Eth_T0 = ThermalEnergy(n, ThetaD, T0); % the thermal energy at standard temperature
    Eth_T = ThermalEnergy(n, ThetaD, T); % the thermal energy at target temperature
    P_T0 = BM_EOS(KT0, KT0P, VV0); % the pressure by compression
    Pth = gamma./V_mol.* (Eth_T-Eth_T0) * 10^(-9); % thermal pressure in GPa
    P_total = Pth + P_T0; % the total pressure

end

% Inverse of the 3rd-order Birch Murnaghan equation of state
function VV0 = Inv_MGD_EOS_RJ(KT0, KT0P, V_lat0, Z, Theta0, gamma0, q, n, P, T0, T)
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
            P1 = MGD_EOS_RJ(KT0, KT0P, V_lat0, Z, Theta0, gamma0, q, n, T0, VV01, T(i2));
            P2 = MGD_EOS_RJ(KT0, KT0P, V_lat0, Z, Theta0, gamma0, q, n, T0, VV02, T(i2));
            while abs(VV01-VV02) > 1E-7
                VV03 = (VV01 + VV02)/2;
                P3 = MGD_EOS_RJ(KT0, KT0P, V_lat0, Z, Theta0, gamma0, q, n, T0, VV03, T(i2));
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

% Mie-Gr?neisen-Debye equation of state based on Tange et al. [JGR 2009] 
function P_total = MGD_EOS_Tange(KT0, KT0P, V_lat0, Z, Theta0, gamma0, a, b, n, T0, VV0, T)
    % V0 : the unit cell volume in A^3
    % KT0 : the isothermal bulk modulus at zero pressure in GPa
    % KT0P : the pressure derivative of KT0
    % Theta0 : the Debye Temperature in K at V/V0 = 1
    % gamma0 : the Gr?neisen parameter at zero pressure
    % a : the volume dependence of Gr?neisen parameter 1
    % b : the volume dependence of Gr?neisen parameter 2
    % n : the number of atoms in the formula unit of the concerned material
    % Z : the number of unit formulas in a unit cell

    NA = 6.0221367E+23; % Avogadro number

    V_lat = VV0 .* V_lat0; % the lattice volume at compression
    V_mol = (NA ./ Z) .* V_lat * 10^(-30); % the molar volume at compression
    gamma = gamma0 .* (1 + a .* (VV0.^b -1)); % the Gr?neisen parameter at compression
    ThetaD = Theta0 .* VV0.^(-(1-a).*gamma0) .* exp(-(gamma-gamma0)/b); % the Debye temperature at compression
    
    Eth_T0 = ThermalEnergy(n, ThetaD, T0); % the thermal energy at standard temperature
    Eth_T = ThermalEnergy(n, ThetaD, T); % the thermal energy at target temperature
    P_T0 = BM_EOS(KT0, KT0P, VV0); % the pressure by compression
    Pth = gamma./V_mol.* (Eth_T-Eth_T0) * 10^(-9); % thermal pressure in GPa
    P_total = Pth + P_T0; % the total pressure

end

% Thermal energy based on the Debye model
% The Debye function is approximated by a 6th-order polynominal
function E_th = ThermalEnergy(n, ThetaD, T);
 % n: the number of atoms in the formula unit of the concerned material
 % ThetaD:  the Debye temperature in K at compression
 % T: temperature in K
    df = [0.33333333 -0.12495878 0.016520035 0.00018490174 ...
        -0.00031034588 0.000035022653 -0.0000013397592]; % The polynominal to approximate the Debye function
    R = 8.314472; % The gas constant in J K-1mol-1 
    x = ThetaD ./ T; % the variable in the Debye function
    E_th = 9 * n .* R .* T .* (df(7)*x.^6 + df(6)*x.^5 + df(5)*x.^4 + df(4)*x.^3 +df(3)*x.^2 + df(2)*x + df(1));
end

% The 3rd-order Birch-Murnaghan equation of state
function P_T0 = BM_EOS(KT0, KT0P, VV0)
% VV0: the relative volume
% P_T0: the pressure at compression
    P_T0 = 3/2 * KT0 .* ((VV0).^(-7/3) - (VV0).^(-5/3)) .* (1 - 3/4*(4 - KT0P) .* ((VV0).^(-2/3) - 1));
end

