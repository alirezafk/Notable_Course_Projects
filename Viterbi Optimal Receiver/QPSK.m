% This function maps 1:4 to the QPSK constellation signals
function Y = QPSK(A)
    [m, n] = size(A);
    for i = 1:m
        for j = 1:n
            switch A(i, j)
                case 1
                    A(i, j) = (1+1i)/sqrt(2);
                case 2
                    A(i, j) = (1-1i)/sqrt(2);
                case 3
                    A(i, j) = (-1+1i)/sqrt(2);
                case 4
                    A(i, j) = (-1-1i)/sqrt(2);
                otherwise
                    A(i, j) = nan;
            end
        end
    end
    Y = A;