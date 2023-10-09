function out = A_map(U, V, X)
    [m, ~] = size(U);
    out = zeros(m, 1);
    
    for i = 1:m
        out(i) = U(i,:)*X*V(i,:)';
    end