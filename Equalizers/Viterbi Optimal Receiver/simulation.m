clc
clear
close all
%% Defining g(t) , c(t) , SNR_range
% Define T=1 and time sampling interval
T = 1;
dt = 1e-1;  % time step

% define g(t)
g_range = 0:dt:T;
g = (1/sqrt(T))*ones(1, length(g_range));

% define c(t)
c_range = 0:dt:2*T;
c = (1/sqrt(2*T))*ones(1, length(c_range));

% Defining SNR values in range [0 20dB]
SNR_min = 0;
SNR_step = 2;
SNR_max = 20;
SNR_db = SNR_min:SNR_step:SNR_max;
SNR = 10.^(SNR_db/10);
%% Obtaining x(t) and h(t)
% h(t) = convolution(g(t), c(t))
[h, h_range] = cross_corr(g, g_range(1), c, c_range(1), dt);
figure
plot(h_range, h);
title("h(t)");

% x(t)
[x, x_range] = cross_corr([h, zeros(1, length(h))], h_range(1), [zeros(1, length(h)), conj(fliplr(h))], -h_range(end), dt);
x = x(1:(2*length(h_range)-1));
x_range = x_range(1:(2*length(h_range)-1));

% plot x(t)
figure
plot(x_range, x);
title("x(t)");
hold on
%% Finding L and xn
T_coeff = 1:(T/dt):length(x);
T_coeff = T_coeff(2:(end-1));
L = (length(T_coeff)-1)/2;
xn = x(T_coeff);
stem(x_range(T_coeff), x(T_coeff));
xn = xn(3:end); % kepping the samples in positive time domain
x0 = xn(1); % define x0 the sample at t=0
xvec = xn(2:end);
hold off
%% First we generate data (based on the modulation) and noise
prompt = "Please select the moulation (BPSK or QPSK):";
dlgtitle = 'Modulation';
answer = inputdlg(prompt, dlgtitle);
modulation = answer{1};

% According on the given modulation we obtain M and states
switch modulation
    case "BPSK"
        M = 2;  % for bpsk M=2
        states = 2*find_states(M, L)-3;
    case "QPSK"
        M = 4;  % for QPSK M=2
        states = QPSK(find_states(M, L));
    otherwise
        errordlg('Invalid Modulation! please run the code again and choose between QPSK or BPSK.','File Error');
        return
end

ber_vals = zeros(1, length(SNR));
depth = 6*L;    % depth of viterbi
num_states = M^L;   % Number of all states
states_table = find_states(M, L);   % This function returns all the states in (M^L*L) matrix
transition_table = transition_states(states_table, M);

N = 1e4*depth;  % Number of all bpsk data
estimated_seq = zeros(1, N);
nsmpl = T/dt;   % Number of samples in duration of T seconds
hlen = length(h);   % length of h
I = generate_data(modulation, N);
rlen = (depth-1)*nsmpl + hlen;

Initial_costs = zeros(num_states, 1);
for i = 1:length(SNR)
    start = 1;
    snr = SNR(i);
    snr_sym = snr*log2(M);
    N0 = 1/(2*x0*snr_sym);
    estimated_seq = zeros(1, N);
    while(start < N)
        I_curr = I(start: (start+depth-1));
        r = zeros(1, rlen);
        
        % In this for we obtain vector r(t)
        for j = 1:depth
            r((1+(j-1)*nsmpl) : ((j-1)*nsmpl+hlen)) = r((1+(j-1)*nsmpl) : ((j-1)*nsmpl+hlen)) + I_curr(j)*h;
        end
        
        % Adding normalized complex noise to r(t)
        noise = sqrt(N0)*(rand(1, rlen) + 1i*rand(1, rlen));
        r = r + noise;
        
        % Finding y[k] using riemann summation form of integral
        y = zeros(1, depth);
        for j = 1:depth
            y(j) = (((hlen-1)*dt)/hlen)*r((1+(j-1)*nsmpl) : ((j-1)*nsmpl+hlen))*h';
        end
        
        % Performing viterbi algorithm
        
        % First Define trellis matrices
        
        % This matrix keeps the costs in each state
        trellis_table_cost = [Initial_costs, zeros(num_states, depth)];
        
        % This matrix store the pathes (edges of graph)
        trellis_table_idx = zeros(num_states, depth);
        
        % Filling previous defined matrices using viterbi algorithm
        for col = 1:depth
            yn = y(col);
            for row = 1:num_states
                In = states(row, end);
                pre = zeros(1, M);
                % In this for we calculate cost of all M possible edges
                % drawn to current state
                for k = 1:M
                    cur_row = transition_table(row, k);
                    pre(k) = trellis_table_cost(cur_row, col) + real(In'*(2*yn-x0*In-2*xvec*fliplr(states(cur_row, :))'));
                end
                % Finding maximum of cost and corresponding index
                [max_vec, idx] = max(pre);
                % Filling trellis matrices
                trellis_table_cost(row, col+1) = max_vec;
                trellis_table_idx(row, col) = transition_table(row, idx);
            end
        end
        % Updating initial costs vector
        Initial_costs = trellis_table_cost(:, end);
        % estimate the sequence using the path that maximizes the cost
        % find index of maximum of last column
        [~, index] = max(trellis_table_cost(:, end));
        detection = zeros(1, depth);
        detection(depth) = states(index, end);
        for j = 1:(depth-1)
            index = trellis_table_idx(index, depth-j+1);
            detection(depth - j) = states(index, end);
        end
        estimated_seq(start: (start+depth-1)) = detection;
        
        % Updating start to iterate over the next block of length depth=6L
        start = start + depth;
    end
    
    % calculating ber of the current snr
    ber_vals(i) = sum(estimated_seq ~= I)/N;
end
%% Plotting BER
figure
if modulation == "BPSK"
    semilogy(SNR_db, ber_vals, 'bo-');
    title("BPSK BER using viterbi algorithm in range [0 20dB]");
else
    semilogy(SNR_db, ber_vals, 'ro-');
    title("QPSK BER using viterbi algorithm in range [0 20dB]");
end
xlabel("SNR(dB)");
ylabel("BER");
grid on