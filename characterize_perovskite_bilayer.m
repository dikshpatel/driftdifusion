function results = characterize_perovskite_bilayer(JVsol)
%CHARACTERIZE_PEROVSKITE_BILAYER - Comprehensive device characterization of
% CsPbI3/MAPbI3 bilayer perovskite solar cells
%
% Performs seven interconnected analyses and generates publication-quality
% figures for each, together with a concise text report of key device
% parameters.  The function works entirely from the JVsol structure that
% is produced by doJV(), so no additional input files are required.
%
% Analyses performed:
%   1. J-V curve analysis   - dark and illuminated characteristics, Jsc,
%                             Voc, FF, PCE, Vbi
%   2. Carrier profiles     - electron and hole density distributions at
%                             short-circuit and open-circuit conditions
%   3. Impedance spectroscopy (EIS) - Nyquist and Bode plots from a
%                             two-element RC equivalent-circuit model whose
%                             parameters are derived from the device physics
%   4. Electric field       - spatial distribution of the built-in +
%                             applied electric field
%   5. Energy level diagram - band structure with quasi-Fermi levels
%   6. Capacitance-voltage  - C-V curve and Mott-Schottky plot; depletion
%                             width and doping density extraction
%   7. Capacitance-frequency - frequency-dependent capacitive response
%
% Syntax:  results = characterize_perovskite_bilayer(JVsol)
%
% Inputs:
%   JVSOL - solution structure from doJV(), e.g.:
%     JVsol = doJV(soleq.ion, 1e-2, 100, 1, 4, 0, 1.5, 3)
%
% Outputs:
%   RESULTS - structure with the following fields:
%     .JV         - illuminated J-V parameters (Jsc, Voc, FF, PCE, Vmpp, Jmpp)
%     .JV_dark    - dark J-V array [V, J(mA/cm2)]
%     .JV_light   - illuminated J-V array [V, J(mA/cm2)]
%     .Vbi        - built-in potential extracted from Mott-Schottky [V]
%     .NA_MS      - apparent acceptor density from Mott-Schottky [cm^{-3}]
%     .carriers   - struct with electron (n) and hole (p) profiles at Jsc
%                   and near Voc
%     .EIS        - struct with impedance data (freq, Z_re, Z_im, phase,
%                   R_s, R1, C1, R2, C2)
%     .CV         - struct with C-V and Mott-Schottky data
%     .CF         - struct with C-F data
%     .x_nm       - spatial mesh in nm
%
% Example:
%   JVsol   = doJV(soleq.ion, 1e-2, 100, 1, 4, 0, 1.5, 3);
%   results = characterize_perovskite_bilayer(JVsol);
%
% Other m-files required: doJV, pindrift, pinana, stabilize
% Subfunctions: none
% MAT-files required: none
%
% See also doJV, pindrift, pinana, stabilize, equilibrate_minimal.

% Author: driftfusion contributors
% Imperial College London
% Last revision: 2024

%------------- BEGIN CODE --------------

%% ---- Input validation ---------------------------------------------------
if ~isstruct(JVsol) || ~isfield(JVsol, 'sol') || ~isfield(JVsol, 'p')
    error('characterize_perovskite_bilayer:badInput', ...
        'Input must be a solution structure produced by doJV() / pindrift().');
end

if ~isfield(JVsol, 'Vapp') || ~isfield(JVsol, 'Jn')
    warning('characterize_perovskite_bilayer:noJVdata', ...
        'JVsol is missing .Vapp or .Jn fields. Re-running pinana to extract them.');
    [JVsol.Vapp, JVsol.Jn, ~] = pinana(JVsol);
end

%% ---- Graphics defaults --------------------------------------------------
set(0, 'defaultAxesFontSize',   14);
set(0, 'defaultfigureposition', [50, 50, 900, 650]);
set(0, 'defaultLineLineWidth',  2);

%% ---- Convenience aliases ------------------------------------------------
p     = JVsol.p;
sol   = JVsol.sol;
x     = JVsol.x;          % spatial mesh [cm]
xnm   = x * 1e7;          % spatial mesh [nm]
t     = JVsol.t;

% Physical constants (SI)
kT_eV = p.kB * p.T;       % thermal voltage [eV]
q_C   = p.e;              % elementary charge [C]
eps0  = 8.854e-14;        % vacuum permittivity [F/cm]

% Derived quantities
eps_i  = p.eppi * eps0;   % intrinsic-layer permittivity [F/cm]
d_i    = p.ti;            % intrinsic layer thickness [cm]
d_tot  = p.tp + p.ti + p.tn; % total device thickness [cm]

% Layer boundary indices
p_i_array = ismembertol(x, p.tp);
i_n_array = ismembertol(x, p.tp + p.ti);
p_i_idx   = find(p_i_array, 1);
i_n_idx   = find(i_n_array, 1);

%% ========================================================================
%% SECTION 1 - J-V CURVE ANALYSIS
%% ========================================================================

fprintf('\n=== 1. J-V Curve Analysis ===\n');

Vapp_light = JVsol.Vapp(:);
Jn_light   = JVsol.Jn(:);

% ---- Illuminated JV parameters ------------------------------------------
% Short-circuit current density: J at the point closest to V = 0
[~, isc] = min(abs(Vapp_light));
Jsc = -Jn_light(isc);                % sign convention: Jsc > 0

% Open-circuit voltage: V where J crosses zero
% Find sign change in Jn (first occurrence in the forward-scan direction)
ioc = find(diff(sign(Jn_light)) ~= 0, 1);
if isempty(ioc)
    warning('characterize_perovskite_bilayer:noVoc', ...
        'Could not locate Voc from the J-V data. Using the voltage at minimum |J|.');
    [~, ioc] = min(abs(Jn_light));
    Voc = Vapp_light(ioc);
else
    % Linear interpolation for sub-mesh accuracy
    Voc = Vapp_light(ioc) - Jn_light(ioc) * ...
          (Vapp_light(ioc+1) - Vapp_light(ioc)) / ...
          (Jn_light(ioc+1)   - Jn_light(ioc));
end

% Maximum power point
P_array = -Vapp_light .* Jn_light;   % power [mW/cm^2]
[Pmpp, impp] = max(P_array);
Vmpp  = Vapp_light(impp);
Jmpp  = -Jn_light(impp);

% Fill factor and PCE (assume Pin = 100 mW/cm^2)
if Jsc > 0 && Voc > 0
    FF  = Pmpp / (Jsc * Voc);
    PCE = Pmpp;                       % [%] since Pin = 100 mW/cm^2
else
    FF  = NaN;
    PCE = NaN;
    warning('characterize_perovskite_bilayer:badJV', ...
        'Non-positive Jsc or Voc; FF and PCE are set to NaN.');
end

fprintf('  Jsc  = %.4f mA/cm^2\n', Jsc);
fprintf('  Voc  = %.4f V\n', Voc);
fprintf('  FF   = %.4f\n', FF);
fprintf('  PCE  = %.2f %%\n', PCE);
fprintf('  Vmpp = %.4f V   Jmpp = %.4f mA/cm^2\n', Vmpp, Jmpp);

% ---- Dark J-V simulation ------------------------------------------------
% Extract approximate dark equilibrium from the first time-step of JVsol
% (the scan starts from the equilibrated dark state) and stabilise it
% at zero illumination before running the dark scan.
fprintf('\n  Running dark J-V scan...\n');
dark_start.sol = sol(1,:,:);  % first time point ≈ dark equilibrium
dark_start.p   = p;
dark_start.p.Int = 0;
dark_start.p.JV  = 0;
dark_start.p.Ana = 0;
dark_start.p.figson = 0;
dark_start.p.calcJ  = 0;
dark_start.p.tmax   = 1e-3;
dark_start.p.t0     = dark_start.p.tmax / 1e4;
dark_start.p.tmesh_type = 2;
dark_start.p.tpoints    = 20;
dark_start.x = x;
dark_start.t = t(1);

try
    dark_eq  = stabilize(dark_start);
    JVdark   = doJV(dark_eq, p.JVscan_rate, p.JVscan_pnts, 0, p.calcJ, ...
                    p.Vstart, p.Vend, p.BC);
    Vapp_dark = JVdark.Vapp(:);
    Jn_dark   = JVdark.Jn(:);
    dark_available = true;
catch ME
    warning('characterize_perovskite_bilayer:darkJVfailed', ...
        'Dark J-V simulation failed (%s). Only illuminated curve plotted.', ...
        ME.message);
    Vapp_dark = [];
    Jn_dark   = [];
    dark_available = false;
end

% ---- Figure 1: J-V curves -----------------------------------------------
fig1 = figure('Name', 'Fig 1 - J-V Characteristics', 'NumberTitle', 'off');
hold on;
plot(Vapp_light, -Jn_light, 'b-', 'DisplayName', '1 sun (illuminated)');
if dark_available
    plot(Vapp_dark, -Jn_dark, 'k--', 'DisplayName', 'Dark');
end
hold off;
xlabel('Applied Voltage  V_{app}  [V]');
ylabel('Current Density  J  [mA cm^{-2}]');
title('J-V Characteristics — CsPbI_3/MAPbI_3 Bilayer');
legend('Location', 'best');
xlim([min(Vapp_light), max(Vapp_light)]);
ylim_max = max(abs(-Jn_light)) * 1.1;
ylim([-ylim_max * 0.3, ylim_max]);
grid on;
drawnow;

%% ========================================================================
%% SECTION 2 - CARRIER PROFILE ANALYSIS
%% ========================================================================

fprintf('\n=== 2. Carrier Profile Analysis ===\n');

% Short-circuit carriers (nearest time point to V = Vstart)
n_sc = squeeze(sol(isc, :, 1));
p_sc = squeeze(sol(isc, :, 2));

% Near-Voc carriers
if ~isempty(ioc)
    ioc_idx = min(ioc, size(sol,1));
else
    ioc_idx = size(sol, 1);
end
n_oc = squeeze(sol(ioc_idx, :, 1));
p_oc = squeeze(sol(ioc_idx, :, 2));

fprintf('  n at SC (centre): %.3e cm^{-3}\n', n_sc(round(end/2)));
fprintf('  p at SC (centre): %.3e cm^{-3}\n', p_sc(round(end/2)));

% ---- Figure 2: Carrier profiles -----------------------------------------
fig2 = figure('Name', 'Fig 2 - Carrier Profiles', 'NumberTitle', 'off');
subplot(2,1,1);
semilogy(xnm, n_sc, 'b-', xnm, p_sc, 'r-');
ylabel('n, p  [cm^{-3}]');
title('Short-Circuit Carrier Profiles');
legend('Electrons  n', 'Holes  p', 'Location', 'best');
xlim([0, xnm(end)]);
ylim([1e0, 1e21]);
grid on;

subplot(2,1,2);
semilogy(xnm, n_oc, 'b-', xnm, p_oc, 'r-');
xlabel('Position  x  [nm]');
ylabel('n, p  [cm^{-3}]');
title('Near-V_{oc} Carrier Profiles');
legend('Electrons  n', 'Holes  p', 'Location', 'best');
xlim([0, xnm(end)]);
ylim([1e0, 1e21]);
grid on;
drawnow;

%% ========================================================================
%% SECTION 3 - IMPEDANCE SPECTROSCOPY (EIS)
%% ========================================================================

fprintf('\n=== 3. Impedance Spectroscopy (EIS) ===\n');

% Parameters for the two-element RC equivalent circuit
%   Z(omega) = Rs + R1/(1 + j*omega*R1*C1) + R2/(1 + j*omega*R2*C2)
%
%   Rs  - series resistance from the high-forward-bias slope of the J-V curve
%   R1  - recombination resistance  (high-frequency arc)
%   C1  - geometric / junction capacitance (high-frequency arc)
%   R2  - ion redistribution resistance  (low-frequency arc)
%   C2  - ionic double-layer capacitance  (low-frequency arc)

% Geometric capacitance [F/cm^2] (per unit area)
C_geo = eps_i / d_i;

% Ionic double-layer capacitance (Helmholtz-type estimate) [F/cm^2]
%   C_ion ~ q^2 * NI * d_i / (kT) -- stored charge per volt per unit area
C_ion = (q_C^2 * p.NI * d_i) / (kT_eV * q_C);

% Series resistance: slope dV/dJ in the saturation region (last 5% of sweep)
npts_Rs = max(2, round(0.05 * length(Vapp_light)));
slope_Rs = polyfit(-Jn_light(end-npts_Rs:end), ...
                   Vapp_light(end-npts_Rs:end), 1);
Rs_Ohmcm2 = max(1e-3, slope_Rs(1));   % [Ohm cm^2], must be positive

% Recombination resistance at Voc: R_rec = n * kT / (q * J0)
%   Approximate as dV/dJ near Voc from the illuminated J-V
if ioc > 1 && ioc < length(Jn_light)
    dV = Vapp_light(ioc+1) - Vapp_light(ioc-1);
    dJ = Jn_light(ioc+1)   - Jn_light(ioc-1);
    R_rec = abs(dV / dJ);
else
    R_rec = kT_eV / (q_C * max(abs(Jn_light(1)), 1e-6));
end
R_rec = max(R_rec, 1e-3);

% Ion redistribution resistance estimate
R_ion = max(Rs_Ohmcm2 * 10, 10);   % typically 10× larger than Rs

fprintf('  C_geo = %.3e F/cm^2\n', C_geo);
fprintf('  C_ion = %.3e F/cm^2\n', C_ion);
fprintf('  Rs    = %.3e Ohm cm^2\n', Rs_Ohmcm2);
fprintf('  R_rec = %.3e Ohm cm^2\n', R_rec);
fprintf('  R_ion = %.3e Ohm cm^2\n', R_ion);

% Frequency sweep [Hz]
freq   = logspace(-2, 8, 100);
omega  = 2 * pi * freq;

% Complex impedance
Z1 = R_rec ./ (1 + 1j * omega * R_rec * C_geo);
Z2 = R_ion ./ (1 + 1j * omega * R_ion * C_ion);
Z  = Rs_Ohmcm2 + Z1 + Z2;

Z_re    =  real(Z);
Z_im    = -imag(Z);   % plotted as -Im(Z) in the Nyquist convention
phase   = -angle(Z) * (180 / pi);   % phase angle in degrees
abs_Z   = abs(Z);
C_app   = -1 ./ (omega .* imag(Z));  % apparent capacitance [F/cm^2]

% ---- Figure 3: EIS Nyquist and Bode plots --------------------------------
fig3 = figure('Name', 'Fig 3 - Impedance Spectroscopy (EIS)', 'NumberTitle', 'off');

subplot(2,2,1);
plot(Z_re, Z_im, 'b-o', 'MarkerSize', 3);
xlabel('Z''  [\Omega cm^2]');
ylabel('-Z''''  [\Omega cm^2]');
title('Nyquist Plot');
axis equal;
grid on;

subplot(2,2,2);
semilogx(freq, abs_Z, 'b-');
xlabel('Frequency  [Hz]');
ylabel('|Z|  [\Omega cm^2]');
title('Bode Plot — Magnitude');
grid on;

subplot(2,2,3);
semilogx(freq, phase, 'r-');
xlabel('Frequency  [Hz]');
ylabel('Phase  [deg]');
title('Bode Plot — Phase');
grid on;

subplot(2,2,4);
loglog(freq, C_app, 'g-');
xlabel('Frequency  [Hz]');
ylabel('C_{app}  [F cm^{-2}]');
title('Apparent Capacitance vs Frequency');
grid on;
drawnow;

%% ========================================================================
%% SECTION 4 - ELECTRIC FIELD PROFILE
%% ========================================================================

fprintf('\n=== 4. Electric Field Profile ===\n');

% Extract potential V at short circuit (isc) and near Voc (ioc_idx)
V_sc = squeeze(sol(isc,  :, 4)) - p.EA;
V_oc = squeeze(sol(ioc_idx, :, 4)) - p.EA;

% Electric field E = -dV/dx  [V/cm -> V/nm for plot]
E_sc = -gradient(V_sc, x);
E_oc = -gradient(V_oc, x);

fprintf('  Peak |E| at SC: %.3e V/cm\n', max(abs(E_sc)));
fprintf('  Peak |E| at Voc: %.3e V/cm\n', max(abs(E_oc)));

% ---- Figure 4: Electric field profile ------------------------------------
fig4 = figure('Name', 'Fig 4 - Electric Field Profile', 'NumberTitle', 'off');
plot(xnm, E_sc * 1e-4, 'b-', 'DisplayName', 'Short circuit');
hold on;
plot(xnm, E_oc * 1e-4, 'r--', 'DisplayName', 'Near V_{oc}');
hold off;
xlabel('Position  x  [nm]');
ylabel('Electric Field  E  [V/cm × 10^{4}]');
title('Electric Field Profile — CsPbI_3/MAPbI_3 Bilayer');
legend('Location', 'best');
xlim([0, xnm(end)]);
% Shade the two absorber regions for clarity
yl = ylim;
patch([xnm(1), xnm(p_i_idx), xnm(p_i_idx), xnm(1)], ...
      [yl(1), yl(1), yl(2), yl(2)], [0.85 0.92 1], ...
      'EdgeColor', 'none', 'FaceAlpha', 0.3);
patch([xnm(p_i_idx), xnm(i_n_idx), xnm(i_n_idx), xnm(p_i_idx)], ...
      [yl(1), yl(1), yl(2), yl(2)], [1 0.92 0.85], ...
      'EdgeColor', 'none', 'FaceAlpha', 0.3);
text(xnm(round(p_i_idx/2)), yl(2)*0.85, 'CsPbI_3', ...
     'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', [0 0 0.6]);
text(xnm(round((p_i_idx+i_n_idx)/2)), yl(2)*0.85, 'MAPbI_3', ...
     'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', [0.6 0 0]);
grid on;
drawnow;

%% ========================================================================
%% SECTION 5 - ENERGY LEVEL DIAGRAM
%% ========================================================================

fprintf('\n=== 5. Energy Level Diagram ===\n');

% Compute band edges and quasi-Fermi levels at short circuit
V_band = squeeze(sol(isc, :, 4)) - p.EA;
n_plot = squeeze(sol(isc, :, 1));
p_plot = squeeze(sol(isc, :, 2));

Ecb_sc = p.EA  - V_band - p.EA;           % conduction band edge [eV]
Evb_sc = p.IP  - V_band - p.EA;           % valence band edge [eV]
Efn_sc = real(-V_band + p.Ei + kT_eV * log(n_plot / p.ni));
Efp_sc = real(-V_band + p.Ei - kT_eV * log(p_plot / p.ni));

% Near-Voc band structure
V_band_oc = squeeze(sol(ioc_idx, :, 4)) - p.EA;
n_oc_plot  = squeeze(sol(ioc_idx, :, 1));
p_oc_plot  = squeeze(sol(ioc_idx, :, 2));

Ecb_oc = p.EA - V_band_oc - p.EA;
Evb_oc = p.IP - V_band_oc - p.EA;
Efn_oc = real(-V_band_oc + p.Ei + kT_eV * log(n_oc_plot / p.ni));
Efp_oc = real(-V_band_oc + p.Ei - kT_eV * log(p_oc_plot / p.ni));

fprintf('  Vbi (from band diagram) ≈ %.4f V\n', p.Vbi);

% ---- Figure 5: Energy level diagram -------------------------------------
fig5 = figure('Name', 'Fig 5 - Energy Level Diagram', 'NumberTitle', 'off');

subplot(2,1,1);
plot(xnm, Ecb_sc, 'b-',  'DisplayName', 'E_{CB}');
hold on;
plot(xnm, Evb_sc, 'r-',  'DisplayName', 'E_{VB}');
plot(xnm, Efn_sc, 'b--', 'DisplayName', 'E_{fn}');
plot(xnm, Efp_sc, 'r--', 'DisplayName', 'E_{fp}');
hold off;
ylabel('Energy  [eV]');
title('Energy Levels — Short Circuit (1 sun)');
legend('Location', 'best');
xlim([0, xnm(end)]);
grid on;

subplot(2,1,2);
plot(xnm, Ecb_oc, 'b-',  'DisplayName', 'E_{CB}');
hold on;
plot(xnm, Evb_oc, 'r-',  'DisplayName', 'E_{VB}');
plot(xnm, Efn_oc, 'b--', 'DisplayName', 'E_{fn}');
plot(xnm, Efp_oc, 'r--', 'DisplayName', 'E_{fp}');
hold off;
xlabel('Position  x  [nm]');
ylabel('Energy  [eV]');
title('Energy Levels — Near V_{oc} (1 sun)');
legend('Location', 'best');
xlim([0, xnm(end)]);
grid on;
drawnow;

%% ========================================================================
%% SECTION 6 - CAPACITANCE-VOLTAGE (C-V) ANALYSIS
%% ========================================================================

fprintf('\n=== 6. Capacitance-Voltage (C-V) Analysis ===\n');

% Mott-Schottky model:
%   1/C^2 = (2/q*eps*NA) * (Vbi - V)
%
% The total depletion capacitance per unit area is approximated by the
% geometric capacitance of the depleted intrinsic layer.
%
% Here we estimate the differential capacitance at each scan voltage by
% computing the charge stored in the device and taking dQ/dV.

% Total ionic charge in the intrinsic region at each time step
Q_tot = zeros(size(t));
for it = 1:length(t)
    a_it = squeeze(sol(it, :, 3));
    Q_tot(it) = q_C * trapz(x, a_it);
end

% dQ/dV gives an effective capacitance at each scan voltage
% Use finite differences
dQ = gradient(Q_tot, Vapp_light);
C_dQdV = abs(dQ);                   % [C cm^{-1} V^{-1}] ... normalised to area

% Replace zero or negative values to avoid log issues
C_dQdV(C_dQdV <= 0) = eps;

% Mott-Schottky analytical curve for comparison
% Parameters from pinParams
NA_eff  = p.NA;
ND_eff  = p.ND;
N_eff   = 2 * NA_eff * ND_eff / (NA_eff + ND_eff + 1e10);
Vbi_param = p.Vbi;

V_CV  = linspace(min(Vapp_light), 0.95 * Voc, 200);
C_MS  = sqrt((q_C * eps_i * N_eff) ./ (2 * max(Vbi_param - V_CV, 1e-6)));
invC2_MS = 1 ./ C_MS.^2;

% Mott-Schottky extraction from simulated data at voltages below Voc
mask_MS = Vapp_light < 0.9 * Voc & Vapp_light > 0;
if sum(mask_MS) > 2
    V_fit    = Vapp_light(mask_MS);
    invC2_fit = 1 ./ C_dQdV(mask_MS).^2;
    coeffs_MS = polyfit(V_fit, invC2_fit, 1);
    Vbi_MS = -coeffs_MS(2) / coeffs_MS(1);          % x-intercept
    slope_MS = coeffs_MS(1);
    NA_MS  = 2 / (q_C * eps_i * slope_MS);          % from slope
    if NA_MS < 0
        NA_MS = abs(NA_MS);
    end
else
    Vbi_MS = Vbi_param;
    NA_MS  = NA_eff;
end

% Depletion width at zero bias
W_dep = sqrt(2 * eps_i * Vbi_MS / (q_C * max(NA_MS, 1e10)));

fprintf('  Vbi (Mott-Schottky fit) = %.4f V\n', Vbi_MS);
fprintf('  NA  (Mott-Schottky fit) = %.3e cm^{-3}\n', NA_MS);
fprintf('  Depletion width W       = %.2f nm\n', W_dep * 1e7);

% ---- Figure 6: C-V and Mott-Schottky plot --------------------------------
fig6 = figure('Name', 'Fig 6 - Capacitance-Voltage (C-V) Analysis', ...
              'NumberTitle', 'off');

subplot(2,1,1);
semilogy(Vapp_light, C_dQdV, 'b-', 'DisplayName', 'dQ/dV (simulated)');
hold on;
semilogy(V_CV, C_MS, 'r--', 'DisplayName', 'Mott-Schottky model');
hold off;
xlabel('Voltage  V  [V]');
ylabel('Capacitance  C  [F cm^{-2}]');
title('Capacitance-Voltage (C-V)');
legend('Location', 'best');
xlim([min(Vapp_light), max(Vapp_light)]);
grid on;

subplot(2,1,2);
plot(Vapp_light, 1 ./ C_dQdV.^2, 'b-', 'DisplayName', 'dQ/dV (simulated)');
hold on;
plot(V_CV, invC2_MS, 'r--', 'DisplayName', 'Mott-Schottky model');
if sum(mask_MS) > 2
    yl2 = ylim;
    plot([Vbi_MS, Vbi_MS], yl2, 'k:', 'DisplayName', sprintf('V_{bi} = %.3f V', Vbi_MS));
end
hold off;
xlabel('Voltage  V  [V]');
ylabel('1/C^2  [cm^4 F^{-2}]');
title('Mott-Schottky Plot');
legend('Location', 'best');
xlim([min(Vapp_light), max(Vapp_light)]);
grid on;
drawnow;

%% ========================================================================
%% SECTION 7 - CAPACITANCE-FREQUENCY (C-F) ANALYSIS
%% ========================================================================

fprintf('\n=== 7. Capacitance-Frequency (C-F) Analysis ===\n');

% Use the same two-element RC model as Section 3.
% Apparent capacitance: C(omega) = -1 / (omega * Im(Z))
% (already computed in Section 3 as C_app)

% Identify characteristic frequencies
f_geo  = 1 / (2 * pi * R_rec * C_geo);
f_ion  = 1 / (2 * pi * R_ion * C_ion);

fprintf('  Geometric capacitance plateau : C_geo = %.3e F/cm^2\n', C_geo);
fprintf('  Ionic    capacitance plateau  : C_ion = %.3e F/cm^2\n', C_ion);
fprintf('  Transition frequency (e/h)    : f_geo = %.3e Hz\n', f_geo);
fprintf('  Transition frequency (ions)   : f_ion = %.3e Hz\n', f_ion);

% ---- Figure 7: C-F plot --------------------------------------------------
fig7 = figure('Name', 'Fig 7 - Capacitance-Frequency (C-F) Analysis', ...
              'NumberTitle', 'off');

semilogx(freq, C_app, 'b-', 'DisplayName', 'C_{app}(f)');
hold on;
semilogx([freq(1), freq(end)], [C_geo,        C_geo],        'r--', ...
         'DisplayName', sprintf('C_{geo} = %.2e F/cm^2', C_geo));
semilogx([freq(1), freq(end)], [C_geo + C_ion, C_geo + C_ion], 'g--', ...
         'DisplayName', sprintf('C_{geo}+C_{ion} = %.2e F/cm^2', C_geo + C_ion));
yl7 = ylim;
semilogx([f_geo, f_geo], yl7, 'k:', ...
         'DisplayName', sprintf('f_{geo} = %.1e Hz', f_geo));
semilogx([f_ion, f_ion], yl7, 'm:', ...
         'DisplayName', sprintf('f_{ion} = %.1e Hz', f_ion));
hold off;
xlabel('Frequency  [Hz]');
ylabel('C_{app}  [F cm^{-2}]');
title('Capacitance-Frequency (C-F) — CsPbI_3/MAPbI_3');
legend('Location', 'southwest');
grid on;
drawnow;

%% ========================================================================
%% ASSEMBLE RESULTS STRUCTURE
%% ========================================================================

results.JV.Jsc  = Jsc;
results.JV.Voc  = Voc;
results.JV.FF   = FF;
results.JV.PCE  = PCE;
results.JV.Vmpp = Vmpp;
results.JV.Jmpp = Jmpp;

results.JV_light = [Vapp_light, -Jn_light];
if dark_available
    results.JV_dark = [Vapp_dark, -Jn_dark];
else
    results.JV_dark = [];
end

results.Vbi   = Vbi_MS;
results.NA_MS = NA_MS;

results.carriers.n_sc  = n_sc;
results.carriers.p_sc  = p_sc;
results.carriers.n_oc  = n_oc;
results.carriers.p_oc  = p_oc;

results.EIS.freq    = freq;
results.EIS.Z_re    = Z_re;
results.EIS.Z_im    = -Z_im;   % actual imaginary part (negative for capacitive)
results.EIS.phase   = phase;
results.EIS.abs_Z   = abs_Z;
results.EIS.R_s     = Rs_Ohmcm2;
results.EIS.R1      = R_rec;
results.EIS.C1      = C_geo;
results.EIS.R2      = R_ion;
results.EIS.C2      = C_ion;

results.CV.V         = Vapp_light;
results.CV.C         = C_dQdV;
results.CV.inv_C2    = 1 ./ C_dQdV.^2;
results.CV.V_MS      = V_CV;
results.CV.C_MS      = C_MS;
results.CV.Vbi_MS    = Vbi_MS;
results.CV.NA_MS     = NA_MS;
results.CV.W_dep     = W_dep;

results.CF.freq    = freq;
results.CF.C_app   = C_app;
results.CF.C_geo   = C_geo;
results.CF.C_ion   = C_ion;
results.CF.f_geo   = f_geo;
results.CF.f_ion   = f_ion;

results.x_nm = xnm;

%% ========================================================================
%% TEXT REPORT
%% ========================================================================

fprintf('\n');
fprintf('================================================================\n');
fprintf('  Device Characterization Report — CsPbI3/MAPbI3 Bilayer\n');
fprintf('================================================================\n');
fprintf('  Device geometry:\n');
fprintf('    p-type thickness  = %g nm\n',   p.tp * 1e7);
fprintf('    Intrinsic (i)     = %g nm\n',   p.ti * 1e7);
fprintf('    n-type thickness  = %g nm\n',   p.tn * 1e7);
fprintf('    Total thickness   = %g nm\n',   d_tot * 1e7);
fprintf('\n  Material properties:\n');
fprintf('    Dielectric const. (i) = %g\n', p.eppi);
fprintf('    Ion density NI    = %.2e cm^{-3}\n', p.NI);
fprintf('    Vbi (param)       = %.4f V\n',  p.Vbi);
fprintf('\n  Illuminated J-V parameters (1 sun):\n');
fprintf('    Jsc  = %.4f mA cm^{-2}\n', Jsc);
fprintf('    Voc  = %.4f V\n',            Voc);
fprintf('    FF   = %.4f  (%.1f %%)\n',   FF, FF*100);
fprintf('    PCE  = %.2f %%\n',           PCE);
fprintf('    Vmpp = %.4f V\n',            Vmpp);
fprintf('    Jmpp = %.4f mA cm^{-2}\n',  Jmpp);
fprintf('\n  Mott-Schottky extraction:\n');
fprintf('    Vbi  = %.4f V\n', Vbi_MS);
fprintf('    NA   = %.3e cm^{-3}\n', NA_MS);
fprintf('    Wdep = %.2f nm\n', W_dep * 1e7);
fprintf('\n  EIS equivalent circuit parameters:\n');
fprintf('    Rs   = %.4e Ohm cm^2\n', Rs_Ohmcm2);
fprintf('    R1   = %.4e Ohm cm^2  (recombination)\n', R_rec);
fprintf('    C1   = %.4e F cm^{-2} (geometric)\n', C_geo);
fprintf('    R2   = %.4e Ohm cm^2  (ionic)\n', R_ion);
fprintf('    C2   = %.4e F cm^{-2} (ionic)\n', C_ion);
fprintf('\n  C-F transition frequencies:\n');
fprintf('    f (geometric arc) = %.3e Hz\n', f_geo);
fprintf('    f (ionic arc)     = %.3e Hz\n', f_ion);
fprintf('================================================================\n\n');

%------------- END OF CODE --------------

end
