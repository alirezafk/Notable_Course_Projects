clc
clear
close all
% Import data
matrix_sco_data;
%% part c
number_of_runs = 10;
max_iter = 1e2;
epsilon = 1e-4;
threshold = epsilon^2;

figure
for num = 1:number_of_runs  
% initiate x, y
x = randn(n_1, 1);
y = randn(n_2, 1);
termination = inf;
obj_vals = zeros(1, max_iter+1);
obj_vals(1) = norm(b - A_map(U, V, x*y.'), 1);
iter = 1;

while((termination > threshold) && (iter <= max_iter))
    cvx_begin quiet
        variables delta_x(n_1) delta_y(n_2)
        minimize(norm(b - A_map(U, V, x*y.') - diag(V*y)*U*delta_x - diag(U*x)*V*delta_y , 1) + 0.5*delta_x.'*U.'*U*delta_x + 0.5*delta_y.'*V.'*V*delta_y);
    cvx_end
    
    % update x , y
    x = x + delta_x;
    y = y + delta_y;
    
    % calculate objective
    obj_vals(iter+1) = norm(b - A_map(U, V, x*y.'), 1);
    termination = norm(delta_x, 2)^2 + norm(delta_y, 2)^2;
    iter = iter + 1;
end

semilogy(1:iter, obj_vals(1:iter));
hold on
end
hold off
xlabel("iteration");
ylabel("error");
title("Error $f(x^k,y^k)$", "interpreter", "latex");
grid on