clear all;

% กำหนดค่า R และ C ที่โหลด
R_load = 25;              % (Ohm) ตัวต้านทานที่โหลด
C_load = 1.32 * 10^(-12); % (F) ตัวเก็บประจุที่โหลด

% องค์ประกอบของวงจรที่ Matching Network (Solution 1) - C1_1 และ C1_2
C1_1 = 5.2 * 10^(-12);   % (F) ตัวเก็บประจุที่ Matching Network (Solution 1) - C1
C1_2 = 0.296 * 10^(-12); % (F) ตัวเก็บประจุที่ Matching Network (Solution 1) - C2

% องค์ประกอบของวงจรที่ Matching Network (Solution 2)
L2 = 1.33 * 10^(-9);   % (H) ตัวเหนี่ยวนำที่ Matching Network (Solution 2)
C2 = 0.759 * 10^(-12); % (F) ตัวเก็บประจุที่ Matching Network (Solution 2)

% ค่า Z0 (Impedance) ของวงจร
Z0 = 30; % (Ohm) characteristic impedance

% กำหนดค่า ZL (โหลด) ตามที่มึงให้มา
ZL = 25 + 15j;  % (Ohm) ค่าโหลดใหม่ที่กำหนด

% ช่วงความถี่ที่ต้องการแสดงผล (8 GHz)
f0 = 8e9;  % 8 GHz
f = linspace(0, 20e9, 1000); % ความถี่ (Hz) 7-9 GHz

%% Solution 1: Matching Network using only C1_1 and C1_2
ZC1_1 = 1 ./ (j*2*pi*f*C1_1);
ZC1_2 = 1 ./ (j*2*pi*f*C1_2);  % Impedance ของ C ใน Matching Network
ZL_C1_1 = ZL + ZC1_1    % series impedance ZL , ZC1_1
Zin1 = ZL_C1_1 .* ZC1_2 ./ (ZL_C1_1 + ZC1_2);  % Impedance เข้า (Solution 1)
Gamma1 = (Zin1 - Z0) ./ (Zin1 + Z0);  % Reflection Coefficient

%% Solution 2: Matching Network using L2 and C2
ZL2 = j*2*pi*f*L2;
ZC2 = 1 ./ (j*2*pi*f*C2);  % Impedance ของ C ใน Matching Network
ZL_C2 = ZL + ZC2    % series impedance ZL , ZC2
Zin2 = ZL_C2 .* ZL2 ./ (ZL_C2 + ZL2);  % Impedance เข้า (Solution 2)
Gamma2 = (Zin2 - Z0) ./ (Zin2 + Z0);  % Reflection Coefficient

% Plotting the reflection coefficients for both solutions
figure;
plot(f/1e9, abs(Gamma1), 'b', 'LineWidth', 2);  % Solution 1 (C1_1 or C1_2)
hold on;
plot(f/1e9, abs(Gamma2), 'r', 'LineWidth', 2);  % Solution 2
legend('Solution 1', 'Solution 2');
xlabel('Frequency (GHz)');
ylabel('|Gamma|');
title('Reflection Coefficient (|Gamma|) vs Frequency for L-Matching Networks');
grid on;

%% คำนวณ Fractional Bandwidth (FBW) สำหรับ Solution 1
idx = find(abs(Gamma1) <= 0.2);
fmin = f(idx(1));
fmax = f(idx(end));
FBW1 = (fmax - fmin) / f0 * 100;
fprintf('FBW1 = %.2f %%\n', FBW1);

%% คำนวณ Fractional Bandwidth (FBW) สำหรับ Solution 2
idx2 = find(abs(Gamma2) <= 0.2);
fmin2 = f(idx2(1));
fmax2 = f(idx2(end));
FBW2 = (fmax2 - fmin2) / f0 * 100;
fprintf('FBW2 = %.2f %%\n', FBW2);

%% check if not found Gamma = 0.2
% ===== Solution 1 =====
idx1 = find(abs(Gamma1) <= 0.2);
if ~isempty(idx1)
    fmin1 = f(idx1(1));
    fmax1 = f(idx1(end));
    FBW1 = (fmax1 - fmin1) / f0 * 100;
    fprintf('FBW1 = %.2f %%\n', FBW1);
else
    FBW1 = 0;
    fprintf('FBW1: No valid frequency range where |Gamma| <= 0.2\n');
end

% ===== Solution 2 =====
idx2 = find(abs(Gamma2) <= 0.2);
if ~isempty(idx2)
    fmin2 = f(idx2(1));
    fmax2 = f(idx2(end));
    FBW2 = (fmax2 - fmin2) / f0 * 100;
    fprintf('FBW2 = %.2f %%\n', FBW2);
else
    FBW2 = 0;
    fprintf('FBW2: No valid frequency range where |Gamma| <= 0.2\n');
end

%% เปรียบเทียบว่าทั้งสอง solution ไหนดีกว่า
if FBW1 > FBW2
    fprintf('Solution 1 has a higher Fractional Bandwidth.\n');
elseif FBW2 > FBW1
    fprintf('Solution 2 has a higher Fractional Bandwidth.\n');
else
    fprintf('Both solutions have equal Fractional Bandwidth.\n');
end