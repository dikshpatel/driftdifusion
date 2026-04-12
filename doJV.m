function JVsol = doJV(sol, scan_rate, scan_points, Int, calcJ, Vstart, Vend, BC)
%DOJV - Perform a current-voltage (J-V) scan simulation using pindrift
% Wraps PINDRIFT with parameters configured for a linear voltage sweep
% from VSTART to VEND at the given SCAN_RATE, returning the full solution
% structure including the applied voltage and current density arrays.
%
% Syntax:  JVsol = doJV(sol, scan_rate, scan_points, Int, calcJ, Vstart, Vend, BC)
%
% Inputs:
%   SOL        - initial solution structure as created by PINDRIFT,
%                EQUILIBRATE or EQUILIBRATE_MINIMAL; the equilibrated dark
%                solution with mobile ions is recommended as the starting
%                point (e.g. sol_i_eq_SR from equilibrate_minimal).
%   SCAN_RATE  - voltage scan rate [V s^{-1}], e.g. 1e-2 for 10 mV/s.
%   SCAN_POINTS - number of (uniformly spaced) time/voltage points in the
%                sweep; more points give a smoother J-V curve but increase
%                computation time.
%   INT        - background illumination intensity [suns]; use 0 for dark
%                and 1 for AM1.5 equivalent (1 sun).
%   CALCJ      - integer flag selecting the current-density calculation
%                method passed to pindrift:
%                  0 = no current calculation
%                  1 = full drift-diffusion at every grid point
%                  2 = continuity-equation method
%                  3 = drift-diffusion at the right boundary only
%                  4 = extraction-velocity boundary method (BC=3 only,
%                      recommended)
%   VSTART     - starting voltage of the scan [V].
%   VEND       - ending voltage of the scan [V].
%   BC         - boundary condition type; must match the value used when
%                SOL was obtained (typically 3 for ohmic contacts).
%
% Outputs:
%   JVSOL - solution structure from PINDRIFT with additional fields:
%     .sol    - solution array (time x position x variable) where
%               variable 1=electrons, 2=holes, 3=ions, 4=potential
%     .p      - parameters structure used for this simulation
%     .t      - time array [s]
%     .x      - spatial mesh [cm]
%     .Vapp   - applied voltage at each time point [V]  (set by PINANA)
%     .Jn     - total current density [mA cm^{-2}] at each time point
%               (set by PINANA when calcJ ~= 0)
%
% Example:
%   % Illuminated forward scan, 10 mV/s, 100 points, 1 sun, BC=3
%   JVsol = doJV(soleq.ion, 1e-2, 100, 1, 4, 0, 1.5, 3)
%
%   % Dark JV scan
%   JVsol_dark = doJV(soleq.ion, 1e-2, 100, 0, 4, 0, 1.5, 3)
%
% Other m-files required: pindrift, pinana
% Subfunctions: none
% MAT-files required: none
%
% See also pindrift, pinana, equilibrate, equilibrate_minimal, changeLight.

% Author: driftfusion contributors
% Imperial College London
% Last revision: 2024

%------------- BEGIN CODE --------------

p = sol.p;

% Suppress interactive figures during the scan to avoid disrupting the
% current workspace; re-enable after the call if needed.
p.figson = 0;

% Enable analysis so that PINDRIFT stores Vapp and Jn in the output struct.
p.Ana = 1;

% Set the current-calculation and boundary-condition modes.
p.calcJ = calcJ;
p.BC    = BC;

% Set illumination intensity (0 = dark, 1 = 1 sun AM1.5).
p.Int = Int;

% Disable pulse perturbation and open-circuit mode for a standard J-V scan.
p.pulseon = 0;
p.OC      = 0;

% Configure J-V scan mode (p.JV = 1 activates the linear voltage ramp in
% PINANA and sets the time mesh via p.tmax).
p.JV          = 1;
p.Vstart      = Vstart;
p.Vend        = Vend;
p.JVscan_rate = scan_rate;
p.JVscan_pnts = scan_points;

% Derive the total simulation time from the scan range and rate, then
% build a uniform (linear) time mesh.
p.tmax       = abs(Vend - Vstart) / scan_rate;
p.t0         = 0;
p.tmesh_type = 1;   % 1 = linear time mesh, required for uniform voltage step
p.tpoints    = scan_points;

% Run the drift-diffusion simulation.
JVsol = pindrift(sol, p);

%------------- END OF CODE --------------

end
