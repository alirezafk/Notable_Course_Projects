% This function finds all possible transtion states of each state from
% previous column in trellis matrices
function table = transition_states(states, M)
    [m, ~] = size(states);
    table = zeros(m, M);
    for i = 1:m
        for j = 1:M
            vec = [j, states(i, 1:(end-1))];
            table(i, j) = find_row(states, vec);
        end
    end