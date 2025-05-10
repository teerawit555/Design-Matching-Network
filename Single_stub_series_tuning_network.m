% adjustable parameters
clear all;
ZL = 25 + j*15; % input impedance (Ohm)
Z0 = 30; % characteristic impedance (Ohm), must be real number
f0 = 8e9;
lambda0 = 3e8/f0;

% program starts here %
YL = 1/ZL; % load ad
GL = real(YL); % conductance of the load
BL = imag(YL); % susceptance of the load
Y0 = 1/Z0; % characteristic admittance of the line

% obtain t
if ( GL ~= Y0 )
% two solution. Put them in a vector 't'
t = ( BL + (-1).^[1 0] * sqrt( GL*( (Y0-GL)^2 + BL^2 )/Y0 ) ) / ( GL - Y0 );
else
% one solution
t = -BL/(2*Y0);
end

% obtain X
X = ( GL^2*t - (Y0-BL*t).*(BL + Y0*t) ) ./ ( Y0*(GL^2 + (BL + Y0*t).^2 ) );

% obtain normalized length, norm_ls = ls/lambda, of the short-circuit stub

norm_ls = -atan( X/Z0 ) / (2*pi); 
norm_ls( norm_ls < 0 ) = norm_ls( norm_ls < 0 ) + 1/2;
ls = norm_ls * lambda0;

% obtain normalized length, norm_lo = lo/lambda, of the open-circuit stub
norm_lo = atan( Z0./X ) / (2*pi);
norm_lo( norm_lo < 0 ) = norm_lo( norm_lo < 0 ) + 1/2;
lo = norm_lo * lambda0;

% obtain the normalized distance, norm_d = d/lambda, of the stub
% Note that t can be a vector, so we can have multiple solutions of norm_d
norm_d = atan( t ) / (2*pi);
norm_d( t<0 ) = norm_d( t< 0 ) + 1/2;
d = norm_d * lambda0;

% Print out the solutionsf
nsol = length( norm_d ); % number of solution
fprintf(1, '[Single-stub series tuner] %d solution(s):', nsol );
for k=1:nsol
fprintf(1, '\nSolution #%d\n', k );
fprintf(1, ' Distance of the stub: d/lambda = %g\n', norm_d(k) );
fprintf(1, ' d= %g\n', d(k) );
fprintf(1, ' Short circuit: ls/lambda = %g\n', norm_ls(k) );
fprintf(1, ' ls = %g\n', ls(k) );
fprintf(1, ' Open circuit: lo/lambda = %g\n', norm_lo(k) );
fprintf(1, ' lo = %g\n', lo(k) );

end

% สร้างช่วงความถี่
f = linspace(5e9, 12e9, 1000); % ความถี่จาก 1 GHz ถึง 15 GHz
lambda = 3e8 ./ f; % ความยาวคลื่นสำหรับแต่ละความถี่
beta = 2 * pi ./ lambda;

% คำนวณ Reflection Coefficient สำหรับแต่ละความถี่
Gamma_short = zeros(nsol, length(f));

for i = 1:nsol
    
    % คำนวณ Reflection Coefficient
    Zin = Z0 .* (ZL + 1j * Z0 .* tan(beta .* d(i))) ./ (Z0 + 1j * ZL .* tan(beta .* d(i)));
    Z_stub = (1)*j * Z0 .* tan(beta * ls(i));

    Z_total = Zin + Z_stub;
    Gamma_short(i, :) = abs((Z_total - Z0) ./ (Z_total + Z0)); % Reflection Coefficient
end

% Plot กราฟ Reflection Coefficient vs Frequency
figure;
plot(f / 1e9, Gamma_short(1,:), 'LineWidth', 2); hold on;
plot(f / 1e9, Gamma_short(2,:), 'LineWidth', 2);
xlabel('Frequency (GHz)');
ylabel('Reflection Coefficient (|\Gamma|)');
title('Single Stub series tunning Short Ckt.');
legend('Short Solution 1', 'Short Solution 2');

grid on;


for i = 1:nsol
    
    % คำนวณ Reflection Coefficient
    Zin = Z0 .* (ZL + 1j * Z0 .* tan(beta .* d(i))) ./ (Z0 + 1j * ZL .* tan(beta .* d(i)));
    Z_stub = (-1)*j * Z0 .* cot(beta * lo(i));

    Z_total = Zin + Z_stub;
    Gamma_open(i, :) = abs((Z_total - Z0) ./ (Z_total + Z0)); % Reflection Coefficient
end

% Plot กราฟ Reflection Coefficient vs Frequency
figure;
plot(f / 1e9, Gamma_open(1,:), 'LineWidth', 2); hold on;
plot(f / 1e9, Gamma_open(2,:), 'LineWidth', 2);
xlabel('Frequency (GHz)');
ylabel('Reflection Coefficient (|\Gamma|)');
title('Single Stub series tunning Open Ckt.');
legend('Open Solution 1', 'Open Solution 2');

grid on;

%% คำนวณ Fractional Bandwidth (FBW) สำหรับ Short Ckt.
for i = 1:2
    idx = find(Gamma_short(i,:) <= 0.2);

    fmin = f(idx(1));
    fmax = f(idx(end));
    FBW = (fmax - fmin)/ f0;

    fprintf('\n Short Ckt. Solution %d FBW \n', i);
    fprintf('FBW = %.2f  %%\n', FBW*100);

end

%% คำนวณ Fractional Bandwidth (FBW) สำหรับ Open Ckt.
for i = 1:2
    idx = find(Gamma_open(i,:) <= 0.2);

    fmin = f(idx(1));
    fmax = f(idx(end));
    FBW = (fmax - fmin)/ f0;

    fprintf('\n Open Ckt. Solution %d FBW \n', i);
    fprintf('FBW = %.2f  %%\n', FBW*100);

end
