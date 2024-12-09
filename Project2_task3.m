%% Task 3: 
rng(100);
% Matrix Setup
n = 25;%= 100;
density = 0.5; 
% Generate Sparse SPD Matrix A and random side vector b
%A = sprandsym(n, density, 0.1, 1); 
M = sprand (n, n, density ); % Sparse random matrix
A = M' * M + speye (n); % SPD matrix
b = rand(n, 1);
%disp(A)
% disp(b)

% Solve using Conjugate Gradient Method
tol = 1e-16;%1e-6; 
max_iter = 3000; 
[x_cg, num_iter, residuals] = conjugate_gradient(A, b, tol, max_iter);

% Solve using MATLAB's direct solver
x_direct = A \ b;
disp("C.G. Solution (first 10 entries):");
disp(x_cg(1:10));%disp(x_cg);
disp("Num iters:");
disp(num_iter);
disp("MatLab Solution (first 10 entries):");
disp(x_direct(1:10));%disp(x_direct); % only display part
disp("Norm difference:");
disp(norm(x_cg - x_direct));

% Plot convergence behavior
figure;
semilogy(1:numel(residuals), residuals, '-o', 'LineWidth', 2);
xlabel('Iteration Number');
ylabel('Residual Norm (log scale)');
title('Convergence History of C.G. Method, n = 25, p = 0.5');
grid on;
set(gcf, 'PaperPositionMode', 'auto');
print('Fig1-t3', '-dpdf', '-r300');
