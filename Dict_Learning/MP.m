function Shat = MP(X, D, N0)
    [~, N] = size(D);
    [~, T] = size(X);
    Shat = zeros(N, T);
    
    % iterating over columns of S
    for t = 1:T
       x = X(:, t);
       for i = 1:N0
           ro = x'*D;
           [~, idx] = max(abs(ro));
           Shat(idx, t) = ro(idx);
           x = x - ro(idx)*D(:, idx);
       end
    end
end