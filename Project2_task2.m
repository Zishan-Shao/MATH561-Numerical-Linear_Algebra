%% Task 2: power iteration and QR algorithm without shifts for eigenvalue computation.
% Matrix Setup
% define a symmetric matrix
rng(100);
n = 10;
max_iter = 1000;
tol = 1e-10;

[Hn, eigen_mat, max_eigen, num_iter] = eigen_comp(n, max_iter, tol);


%disp(Hn);
%disp(eigen_mat);
%disp(max_eigen);
disp("Hessenberg Matrix (H):");
disp(Hn);
disp("Eigenvalues (QR Algorithm):");
disp(diag(eigen_mat));
disp("Largest Eigenvalue (Power Iteration):");
disp(max_eigen);
disp("Number of iterations (Power Iteration):");
disp(num_iter);

% Compare with MATLAB's eig function
eigen_true = sort(eig(Hn), 'descend'); 
%disp("Eigenvalues (MATLAB eig):");
%disp(eigen_true);
eigen_qr = diag(eigen_mat); 
diff_qr = norm(sort(eigen_qr, 'descend') - eigen_true);
disp("Difference between QR and MATLAB eig:");
disp(diff_qr);

diff_power = abs(max_eigen - eigen_true(1));
disp("Difference between Power Iteration and MATLAB eig:");
disp(diff_power);

%% Convergence Analysis
converge_qr = []; % Store successive differences (QR algorithm)
H = Hn; % Reinitialize Hessenberg matrix
for k = 1:num_iter
    [Q, R] = qr(H);
    H_new = R * Q; % Update Hessenberg matrix
    converge_qr = [converge_qr; norm(diag(H_new) - diag(H))]; % Record differences
    H = H_new;
end

figure;
plot(1:length(converge_qr), converge_qr, '-o', 'LineWidth', 2);
xlabel('Iteration Number');
ylabel('Difference Between Successive Eigenvalues');
title('Convergence Pattern of QR Algorithm');
grid on;
set(gcf, 'PaperPositionMode', 'auto'); 
print('Fig1-t2', '-dpdf', '-r300'); 

%% Matrix Size effects
matrix_sizes = [5, 10, 25, 50, 75, 100, 150, 250, 500]; % Varying sizes
converge_iter = []; 
for n = matrix_sizes
    [~, ~, ~, iter] = eigen_comp(n, max_iter, tol);
    converge_iter = [converge_iter; iter];
end

% Plot Effect of Matrix Size on Convergence
figure;
plot(matrix_sizes, converge_iter, '-o', 'LineWidth', 2);
xlabel('Matrix Size (n)');
ylabel('Number of Iterations for Convergence');
title('Effects of Matrix Size on Convergence');
grid on;
set(gcf, 'PaperPositionMode', 'auto'); 
print('Fig2-t2', '-dpdf', '-r300');


