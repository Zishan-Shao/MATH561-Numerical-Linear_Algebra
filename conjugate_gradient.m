function [x, num_iter, res] = conjugate_gradient(A, b, tol, max_iter)

    n = length(b); 
    x = zeros(n, 1); 
    r = b - A * x; 
    p = r;
    res = norm(r); 
    num_iter = 0;

    for k = 1:max_iter
        Ap = A * p;
        alpha = (r' * r) / (p' * Ap); 
        x = x + alpha * p; 
        r_new = r - alpha * Ap; 
        residual_norm = norm(r_new); 
        res = [res; residual_norm]; % Store residual norm
        if residual_norm < tol
            num_iter = k;
            break;
        end
        % computes beta
        beta = (r_new' * r_new) / (r' * r);
        p = r_new + beta * p; 
        r = r_new; 
    end

    if num_iter == 0
        num_iter = max_iter;
    end
end
