clc
clear
close all
%% Defining SNR values in range [0 20dB]
SNR_min = 0;
SNR_step = 10;
SNR_max = 100;
SNR_db = SNR_min:SNR_step:SNR_max;
SNR = 10.^(SNR_db/10);
%% Defining x[n] and innovation filter
x = [1/12, 0.5, 5/6, 0.5, 1/12];
x0 = 5/6;
L = 2;
Xz_roots = roots(x.');
syms z Xz
Xz = (1/12)*z^(-2)+0.5*z^(-1)+(5/6)+0.5*z + (1/12)*z^2;

% hence X(z) = F(z)F(1/z)
% where F(z) = (1+z^-1)(1+0.26792*z^-1)
f = [1, 1-Xz_roots(end), -Xz_roots(end)];
f = sqrt(x(1) /f(end))*f;
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

N = 1e6*depth;  % Number of all bpsk data
estimated_seq = zeros(1, N);
I = generate_data(modulation, N);

Initial_costs = zeros(num_states, 1);
for i = 1:length(SNR)
    start = 1;
    snr = SNR(i);
    snr_sym = snr*log2(M);
    N0 = 1/(2*snr_sym);
    estimated_seq = zeros(1, N);
    while(start < N)
        I_curr = I(start: (start+depth-1));
        noise = sqrt(N0)*(randn(1,depth+L)+1i*randn(1,depth+L));
        noise_out = conv(noise, f);
        nu = noise_out((L+1):(end-L));
        y = conv(x, I_curr);
        y = y((L+1):(end-L)) + nu;
        
        
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
                    pre(k) = trellis_table_cost(cur_row, col) + real(In'*(2*yn-x0*In-2*x((end-L+1):end)*fliplr(states(cur_row, :))'));
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
    title("BPSK BER using viterbi algorithm in range [0 15dB]");
else
    semilogy(SNR_db, ber_vals, 'ro-');
    title("QPSK BER using viterbi algorithm in range [0 20dB]");
end
xlabel("SNR(dB)");
ylabel("BER");
grid on