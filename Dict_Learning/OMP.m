function Shat = OMP(X, D, N0)
    [~, N] = size(D);
    [~, T] = size(X);
    Shat = zeros(N, T);
    
    % iterating over columns of S
    for t = 1:T
       x = X(:, t);
       posOMP = zeros(1,N0);
       for i=1:N0
           ro = x'*D;
           [~, posOMP(i)] = max(abs(ro));
           if i>1
               Dsub = D(:,posOMP(1:i));
               Shat(posOMP(1:i), t) = pinv(Dsub)*X(:, t);
               x = X(:, t)-D*Shat(:, t);
           else
               Shat(posOMP(1), t) = ro(posOMP(1));
               x = x-Shat(posOMP(1), t)*D(:,posOMP(1));
           end
       end
    end