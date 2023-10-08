% This function finds the index of the row a in the given matrix A
function idx = find_row(A, a)
    [m, n] = size(A);
    for i = 1:m
        if sum(a == A(i, :)) == n
            idx = i;
        end
    end