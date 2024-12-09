function[R] = householder(A)
    [m,n] = size(A);
    for k = 1:m-2
        x = A(k+1:m,k);
        vk = sign(x(1)) * norm(x,2) * eye(m - k,1) + x;
        vk = vk / norm(vk, 2);
        A(k+1:m, k:m)= A(k+1:m, k:m) - 2*vk * (vk' * A(k+1:m, k:m));
        A(1:m, k+1:m) = A(1:m, k+1:m) - 2 * (A(1:m, k+1:m)*vk)*vk';
    end

    R = A;
end