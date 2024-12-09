%% Task 1: compare single global polynomial fit with piecewise polynomial fits
% Data Generation
rng(100); % set seed
n = 30;
f = @(x) sin(2 * pi * x); % smooth function
x = linspace(-1,1,30); % 30 observations with equal distance

delta = 0.05; % noise level
y = f(x) + delta * (2 * rand(size(x)) - 1); % uniform noises


%% Global Polynomial Fit (For Comparison)
%A = fliplr(vander(y)); 
%A = A(:, 2:5); % columns for degree 2-5
A = fliplr(vander(x)); % create the vandermode matrix
A1 = A(:, 1:3);  % degree 1-2
A2 = A2(:, 1:6); % degree 1-5
A3 = A3(:, 1:26); % degree 1-25

% Degree 2
coeff_1 = A1 \ y(:); % Coefficients for degree 2
p1_fit = A1 * coeff_1; % Evaluate the polynomial for degree 2
rmse_p1 = sqrt(mean((f(x) - p1_fit).^2, 'all')); % RMSE for degree 2

%disp(A1)

% Degree 5
coeff_2 = A2 \ y(:); % Coefficients for degree 5
p2_fit = A2 * coeff_2; % Evaluate the polynomial for degree 5
rmse_p2 = sqrt(mean((f(x) - p2_fit).^2, 'all')); % RMSE for degree 5

% Degree 25
coeff_3 = A3 \ y(:); % Coefficients for degree 25
p3_fit = A3 * coeff_3; % Evaluate the polynomial for degree 25
rmse_p3 = sqrt(mean((f(x) - p3_fit).^2, 'all')); % RMSE for degree 25

% Display RMSE Values
disp('RMSE values for each fits:');
fprintf('degree m = 2: RMSE = %.4f\n', rmse_p1);
fprintf('degree m = 5: RMSE = %.4f\n', rmse_p2);
fprintf('degree m = 25: RMSE = %.4f\n', rmse_p3);


%% Piecewise polynomial fits
% Divide the interval [-1; 1] into k = 3; 6; 10 equal subintervals. For each subinterval: 
% 1. Use only the data points that fall within the subinterval.
% 2. Fit a polynomial of degree m = 2 to the data in that subinterval using the least-squares method.
% Combine the polynomials to create a piecewise function for the entire interval. 
% Compute the RMSE of the piecewise fit for each k = 3; 6; 10.

% k = 3
k = 3;
m = 2; % this is the degree
fits_k1 = zeros(size(y)); % array to store fitted value for k = 3
intv = linspace(-1,1,k+1); % so we have 4 pivot points that defines three equal spaced intervals
coeff_k1 = zeros(k, m + 1);
for j = 1:k
    
    %idx = x((x >= intv(1)) & (x <= intv(1)));
    x_min = intv(j);
    x_max = intv(j+1);
    idx = (x >= x_min) & (x <= x_max); % Logical indexing for points in the subinterval
    A_sub = A(idx, 1:(m+1));
    y_sub = y(idx);
    coeff_sub = A_sub \ y_sub(:);
    %coeff_k1(idx) = coeff_sub;
    fits_k1(idx) = A_sub * coeff_sub;
    coeff_k1(j, :) = coeff_sub';
end
%disp('Fitted values:');
%disp(fits_k1);

% k = 6
k = 6;
fits_k2 = zeros(size(y)); % array to store fitted value for k = 3
intv = linspace(-1,1,k+1); % so we have 4 pivot points that defines three equal spaced intervals
coeff_k2 = zeros(k, m + 1);
for j = 1:k
    
    %idx = x((x >= intv(1)) & (x <= intv(1)));
    x_min = intv(j);
    x_max = intv(j+1);
    idx = (x >= x_min) & (x <= x_max); % Logical indexing for points in the subinterval
    A_sub = A(idx, 1:(m+1));
    y_sub = y(idx);
    coeff_sub = A_sub \ y_sub(:);
    %coeff_k1(idx) = coeff_sub;
    fits_k2(idx) = A_sub * coeff_sub;
    coeff_k2(j, :) = coeff_sub';
end

% k = 6
k = 10;
fits_k3 = zeros(size(y)); % array to store fitted value for k = 3
intv = linspace(-1,1,k+1); % so we have 4 pivot points that defines three equal spaced intervals
coeff_k3 = zeros(k, m + 1);
for j = 1:k
    
    %idx = x((x >= intv(1)) & (x <= intv(1)));
    x_min = intv(j);
    x_max = intv(j+1);
    idx = (x >= x_min) & (x <= x_max); % Logical indexing for points in the subinterval
    A_sub = A(idx, 1:(m+1));
    y_sub = y(idx);
    coeff_sub = A_sub \ y_sub(:);
    %coeff_k1(idx) = coeff_sub;
    fits_k3(idx) = A_sub * coeff_sub;
    coeff_k3(j, :) = coeff_sub';
end

disp('RMSE values for piecewise fits:');
fprintf('k = 3: RMSE = %.4f\n', sqrt(mean((f(x) - fits_k1).^2, 'all')));
fprintf('k = 6: RMSE = %.4f\n', sqrt(mean((f(x) - fits_k2).^2, 'all')));
fprintf('k = 10: RMSE = %.4f\n', sqrt(mean((f(x) - fits_k3).^2, 'all')));


%% Analysis
% Plot Comparison of Global and Piecewise Fits
figure;
hold on;
plot(x, f(x), 'k-', 'LineWidth', 1.5, 'DisplayName', 'True Function');
plot(x, y, 'o', 'MarkerSize', 5, 'DisplayName', 'Noisy Data'); % noise data
% polynomial fits with whole dat
plot(x, p1_fit, 'b--', 'LineWidth', 1, 'DisplayName', 'Global Fit (m = 2)');
plot(x, p2_fit, 'r--', 'LineWidth', 1, 'DisplayName', 'Global Fit (m = 5)');
plot(x, p3_fit, 'g--', 'LineWidth', 1, 'DisplayName', 'Global Fit (m = 25)');
% piecewise fits
plot(x, fits_k1, 'c-', 'LineWidth', 1, 'DisplayName', 'Piecewise Fit (k = 3)');
plot(x, fits_k2, 'm-', 'LineWidth', 1, 'DisplayName', 'Piecewise Fit (k = 6)');
plot(x, fits_k3, 'y-', 'LineWidth', 1, 'DisplayName', 'Piecewise Fit (k = 10)');

legend('show');
title('Comparison of Global and Piecewise Polynomial Fits');
xlabel('x');
ylabel('y');
grid on;
set(gcf, 'PaperPositionMode', 'auto');
print('Fig1', '-dpdf', '-r300');
hold off;