clear all;

%% ===== LOAD PARAMETERS =====
ZL = 25 + 1j*15;         % Load impedance (Ohm)
RL = real(ZL);           % Real part of load
XL = imag(ZL);           % Imaginary part of load

%% ===== LINE PARAMETERS =====
Z0 = 30;                 % Characteristic impedance (Ohm)
Y0 = 1 / Z0;             % Characteristic admittance (S)

%% ===== FREQUENCY PARAMETERS =====
f0 = 8e9;                % Design frequency (Hz)
c = 3e8;                 % Speed of light (m/s)
lambda0 = c / f0;        % Wavelength at f0 (m)
beta0 = 2 * pi / lambda0; % Propagation constant (rad/m)

%% ===== CALCULATE t AND B =====
t = ( XL + (-1).^[0 1] .* sqrt( RL * ((Z0 - RL)^2 + XL^2) / Z0 ) ) / ( RL - Z0 );
B = ( RL^2 .* t - (Z0 - XL .* t) .* (XL + Z0 .* t) ) ./ ...
    ( Z0 * (RL^2 + (XL + Z0 .* t).^2) );

%% ===== NORMALIZED STUB LENGTHS =====

% --- Compute normalized distance d from load to stub ---
% From: t = tan(2πd) => d = (1/2π) * atan(t)
norm_d = atan(t) / (2*pi);
norm_d(norm_d < 0) = norm_d(norm_d < 0) + 0.5;

% --- Compute normalized length ls of short-circuit stub ---
% From: B = Y0 * tan(2πl) => l = (1/2π) * atan(Y0/B)
norm_ls = atan(Y0 ./ B) / (2*pi);
norm_ls(norm_ls < 0) = norm_ls(norm_ls < 0) + 0.5;
norm_ls = mod(norm_ls, 0.5);  % restrict to [0, 0.5)

%% ===== PHYSICAL LENGTHS (meters) =====
d = norm_d * lambda0;   
ls = norm_ls * lambda0; 

%% ===== PRINT SUMMARY =====
nsol = length(norm_d);  % number of matching solutions
fprintf('\n========== Short-stub Shunt Tuner ==========\n');
fprintf('Total Matching Solutions: %d\n', nsol);

for i = 1:nsol
    fprintf('\n----------------------------------------------\n');
    fprintf('Solution %d\n', i);
    fprintf('----------------------------------------------\n');
    
    % Theoretical results
    fprintf('>> THEORETICAL CALCULATION\n');
    fprintf('t  = %12.5f\n', t(i));
    fprintf('B  = %12.5f S\n\n', B(i));
    
    % Physical results
    fprintf('>> PHYSICAL IMPLEMENTATION\n');
    fprintf('d  = %12.5f m   %8.5f λ\n', d(i), norm_d(i));
    fprintf('ls = %12.5f m   %8.5f λ\n', ls(i), norm_ls(i));
end

%% ===== FREQUENCY SWEEP ANALYSIS =====
f = linspace(0, 11e9, 1000);
lambda = c ./ f;
beta = 2*pi ./ lambda;
Gamma_short = zeros(nsol, length(f));

for i = 1:nsol
    Zin_d = Z0 * (ZL + 1j*Z0 * tan(beta * d(i))) ./ ...
                 (Z0 + 1j * ZL .* tan(beta * d(i)));
    Yin_stub = -1j * Y0 * cot(beta .* ls(i));
    Yin_total = 1 ./ Zin_d + Yin_stub;
    Zin_total = 1 ./ Yin_total;
    Gamma_short(i, :) = abs((Zin_total - Z0) ./ (Zin_total + Z0));
end

%% ===== PLOT REFLECTION COEFFICIENT =====
figure;
plot(f/1e9, Gamma_short(1,:), 'b', 'LineWidth', 2); hold on;
if nsol > 1
    plot(f/1e9, Gamma_short(2,:), 'r', 'LineWidth', 2);
end
yline(0.2, 'k--', '|\Gamma| = 0.2');
xlabel('Frequency (GHz)');
ylabel('|\Gamma|');
legend('Short Stub - Solution 1', 'Short Stub - Solution 2', 'Location', 'Best');
title('Short Stub Matching ');
grid on;

%% ===== FRACTIONAL BANDWIDTH CALCULATION =====
for i = 1:nsol
    idx = find(Gamma_short(i,:) < 0.2);
    if ~isempty(idx)
        fmin = f(idx(1));
        fmax = f(idx(end));
        fbw = (fmax - fmin) / f0;
        fprintf('\n[Short Stub Solution %d - FBW Summary]\n', i);
        fprintf('FBW = %.2f %%\n', fbw * 100);
    end
end
