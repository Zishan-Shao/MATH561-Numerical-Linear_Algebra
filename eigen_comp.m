function [Hn, eigen_mat, max_eigen, num_iter] = eigen_comp(n, max_iter, tol)

    % Matrix Setup
    B = rand(n, n);
    A = B * B'; % Symmetric positive definite matrix
    
    % Reduce to Hessenberg Form
    H = householder(A);
    H = round(H, 14);
    %disp("Verify if Hessenberg Matrix:");
    %disp(H); % Hessenberg matrix
    Hn = H;
    
    % Compute Eigenvalues with QR Algorithm
    %max_iter = 1000; 
    %tol = 1e-10; 
    k = 0; 
    
    while k < max_iter
        % QR Decomposition
        [Q, R] = qr(H);
        H_new = R * Q; % Update Hessenberg matrix
        if norm(H_new - diag(diag(H_new)), 'fro') < tol
            break; % Stop iteration if off-diagonal near zero
        end
        H = H_new;
        k = k + 1;
    end

    eigen_mat = H_new;
    
    % Power Iteration
    x = rand(n, 1); 
    x = x / norm(x,2); 
    eigen_old = 0; % Initialize eigenvalue
    
    for k = 1:max_iter
        
        x_new = H * x;
        x_new = x_new / norm(x_new);
        eigen_new = x_new' * H * x_new; % Rayleigh quotient
        
        if abs(eigen_new - eigen_old) < tol
            break;
        end
        x = x_new;
        eigen_old = eigen_new;
    end
    max_eigen = eigen_new;
    num_iter = k;
    
    % Display results
    %disp("Largest Eigenvalue (Power Iteration):");
    %disp(eigen_new);
    %disp("Number of iterations:");
    %disp(k);

   
end
