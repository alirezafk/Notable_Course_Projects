function states = find_states(M, L)
    if L == 1
        states = (1:M)';
    else
        U = find_states(M, L-1);
        S = [U, ones(M^(L-1), 1)];
        for i = 2:M
            S = [S; [U, i*ones(M^(L-1), 1)]];
        end
        states = S;
    end
    