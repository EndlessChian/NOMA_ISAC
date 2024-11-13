
% MATLAB Code for Enhanced NOMA-ISAC and FDSAC Simulation
% Including Dynamic Power Allocation and Sensing Rate vs. Sum Communication Rate Plot

% Clear workspace and command window
clear; clc;

%% 1. Define System Parameters
% Power Allocation Factors (Dynamic for main simulation)
% Initially set to fixed values; will vary dynamically in trials
alpha_N_fixed = 0.2;             % Base power allocation for Near User (NU) in ISAC
alpha_F_fixed = 0.8;             % Base power allocation for Far User (FU) in ISAC

% Noise Variances
sigma_c = 0.01;                  % Variance of noise for communication
sigma_s = 0.01;                  % Variance of noise for sensing

% Target Rates
target_rate_N = 0.8;             % Target rate for NU in bps/Hz
target_rate_F = 0.8;             % Target rate for FU in bps/Hz

% Monte Carlo Simulation Parameters
num_trials = 10000;              % Number of trials for Monte Carlo simulation

% SNR Ranges
SNR_dB_outage = 0:5:60;          % SNR range in dB for outage probability and ECR analysis
SNR_dB_sensing = -20:5:20;        % SNR range in dB for sensing rate analysis

% FDSAC Parameters
kappa_fixed = 0.5;                % Fixed fraction of bandwidth allocated to communications
mu_fixed = 0.5;                   % Fixed fraction of power allocated to communications

% Sensing Parameters
L = 30;                           % Length of communication frame
M = 8;                            % Number of range cells for sensing
lambda_mean = [5, 3, 3.5, 2.5, 1.5, 2, 1, 0.5];  % Mean eigenvalues for R
num_eigen = length(lambda_mean);  % Number of eigenvalues

% Convert Target Rates to SINR Thresholds
gamma_target_N = 2^target_rate_N - 1;  % SINR threshold for NU
gamma_target_F = 2^target_rate_F - 1;  % SINR threshold for FU

%% 2. Initialize Arrays to Store Results
% Outage Probabilities
OP_NU_ISAC = zeros(1, length(SNR_dB_outage));
OP_FU_ISAC = zeros(1, length(SNR_dB_outage));
OP_NU_FDSAC = zeros(1, length(SNR_dB_outage));
OP_FU_FDSAC = zeros(1, length(SNR_dB_outage));

% Sensing Rates
sensing_rate_ISAC_sim = zeros(1, length(SNR_dB_sensing));
sensing_rate_FDSAC_sim = zeros(1, length(SNR_dB_sensing));

% Ergodic Communication Rates (ECR)
ECR_NU_ISAC = zeros(1, length(SNR_dB_outage));
ECR_FU_ISAC = zeros(1, length(SNR_dB_outage));
ECR_NU_FDSAC = zeros(1, length(SNR_dB_outage));
ECR_FU_FDSAC = zeros(1, length(SNR_dB_outage));

%% 3. Outage Probability and ECR Simulation with Dynamic Power Allocation
for idx = 1:length(SNR_dB_outage)
    SNR_linear = 10^(SNR_dB_outage(idx) / 10);  % Convert SNR from dB to linear scale
    p = SNR_linear * sigma_c^2;                 % Adjust power budget for each SNR value
    
    % Initialize Counters for Outage Probabilities
    outage_count_NU_ISAC = 0;
    outage_count_FU_ISAC = 0;
    outage_count_NU_FDSAC = 0;
    outage_count_FU_FDSAC = 0;
    
    % Initialize Accumulators for Ergodic Communication Rates
    sum_RN_ISAC = 0;
    sum_RF_ISAC = 0;
    sum_RN_FDSAC = 0;
    sum_RF_FDSAC = 0;
    
    for trial = 1:num_trials
        % Dynamic Power Allocation based on Real-Time Channel Feedback
        alpha_N = alpha_N_fixed + 0.1 * rand();     % alpha_N ∈ [0.2, 0.3]
        alpha_F = alpha_F_fixed - 0.1 * rand();     % alpha_F ∈ [0.5, 0.8]
        
        % Ensure that alpha_N + alpha_F <= 1 for valid power allocation
        if (alpha_N + alpha_F) > 1
            alpha_F = 1 - alpha_N;
        end
        
        % Generate Rayleigh Fading Channels for NU and FU
        hN = sqrt(0.9) * (randn(1,1) + 1i * randn(1,1)) / sqrt(2);  % Channel for NU
        hF = sqrt(0.2) * (randn(1,1) + 1i * randn(1,1)) / sqrt(2);  % Channel for FU
        
        % Calculate SINRs for ISAC
        gamma_sic_NU_ISAC = (p * abs(hN)^2 * alpha_F) / (sigma_c^2 + p * abs(hN)^2 * alpha_N);
        gamma_NU_ISAC = (p * abs(hN)^2 * alpha_N) / sigma_c^2;
        gamma_FU_ISAC = (p * abs(hF)^2 * alpha_F) / (sigma_c^2 + p * abs(hF)^2 * alpha_N);
        
        % Determine Outage for ISAC
        if ~(gamma_sic_NU_ISAC > gamma_target_F && gamma_NU_ISAC > gamma_target_N)
            outage_count_NU_ISAC = outage_count_NU_ISAC + 1;
        end
        if gamma_FU_ISAC < gamma_target_F
            outage_count_FU_ISAC = outage_count_FU_ISAC + 1;
        end
        
        % Calculate SINRs for FDSAC (Corrected: Removed division by kappa)
        gamma_NU_FDSAC = (mu_fixed * p * abs(hN)^2) / sigma_c^2;
        gamma_FU_FDSAC = (mu_fixed * p * abs(hF)^2) / sigma_c^2;
        
        % Determine Outage for FDSAC
        if gamma_NU_FDSAC < gamma_target_N
            outage_count_NU_FDSAC = outage_count_NU_FDSAC + 1;
        end
        if gamma_FU_FDSAC < gamma_target_F
            outage_count_FU_FDSAC = outage_count_FU_FDSAC + 1;
        end
        
        % Accumulate Ergodic Communication Rates for ISAC
        sum_RN_ISAC = sum_RN_ISAC + log2(1 + gamma_NU_ISAC);
        sum_RF_ISAC = sum_RF_ISAC + log2(1 + gamma_FU_ISAC);
        
        % Accumulate Ergodic Communication Rates for FDSAC (scaled by kappa)
        sum_RN_FDSAC = sum_RN_FDSAC + kappa_fixed * log2(1 + gamma_NU_FDSAC);
        sum_RF_FDSAC = sum_RF_FDSAC + kappa_fixed * log2(1 + gamma_FU_FDSAC);
    end
    
    % Compute Outage Probabilities
    OP_NU_ISAC(idx) = outage_count_NU_ISAC / num_trials;
    OP_FU_ISAC(idx) = outage_count_FU_ISAC / num_trials;
    OP_NU_FDSAC(idx) = outage_count_NU_FDSAC / num_trials;
    OP_FU_FDSAC(idx) = outage_count_FU_FDSAC / num_trials;
    
    % Compute Ergodic Communication Rates (Average over Trials)
    ECR_NU_ISAC(idx) = sum_RN_ISAC / num_trials;
    ECR_FU_ISAC(idx) = sum_RF_ISAC / num_trials;
    ECR_NU_FDSAC(idx) = sum_RN_FDSAC / num_trials;
    ECR_FU_FDSAC(idx) = sum_RF_FDSAC / num_trials;
end

%% 4. Sensing Rate Simulation
for idx = 1:length(SNR_dB_sensing)
    SNR_linear = 10^(SNR_dB_sensing(idx) / 10);  % Convert SNR from dB to linear scale
    p = SNR_linear * sigma_s^2;                 % Adjust power budget for each SNR value
    
    % Initialize Accumulators for Sensing Rates
    sum_Rs_ISAC = 0;
    sum_Rs_FDSAC = 0;
    
    for trial = 1:num_trials
        % Introduce Randomness in Correlation Matrix R by Perturbing Eigenvalues
        perturbation_std = 0.1;                                 % Standard deviation for perturbation
        lambda_perturbed = lambda_mean + perturbation_std * randn(1, num_eigen);  % Perturbed eigenvalues
        lambda_perturbed(lambda_perturbed < 0) = 0;              % Ensure non-negative eigenvalues
        
        % Construct the Random Correlation Matrix R
        R_random = diag(lambda_perturbed);
        
        %% ISAC Sensing Rate Calculation
        Rs_ISAC = (1 / L) * log2(det(eye(M) + (p * L / sigma_s^2) * R_random));
        sum_Rs_ISAC = sum_Rs_ISAC + Rs_ISAC;
        
        %% FDSAC Sensing Rate Calculation
        Rs_FDSAC = ((1 - kappa_fixed) / L) * log2(det(eye(M) + ((1 - mu_fixed) * p * L / ((1 - kappa_fixed) * sigma_s^2)) * R_random));
        sum_Rs_FDSAC = sum_Rs_FDSAC + Rs_FDSAC;
    end
    
    % Average Sensing Rates over Trials
    sensing_rate_ISAC_sim(idx) = sum_Rs_ISAC / num_trials;
    sensing_rate_FDSAC_sim(idx) = sum_Rs_FDSAC / num_trials;
end

%% 5. Sensing Rate vs. Sum Communication Rate Plot at p = 5 dB
% Fixed SNR for this plot
p_fixed_dB = 5;
p_fixed = 10^(p_fixed_dB / 10) * sigma_c^2;  % Convert to linear scale and adjust power
    
% Define Range for mu (Power Allocation Factor) from 0 to 1
mu_range = 0:0.05:1;            % Vary mu from 0 to 1 in steps of 0.05
num_mu = length(mu_range);

% Initialize Arrays to Store Sensing Rate and Sum ECR for FDSAC
sensing_rate_FDSAC_mu = zeros(1, num_mu);
sum_ECR_FDSAC_mu = zeros(1, num_mu);

% Initialize Arrays to Store Sensing Rate and Sum ECR for ISAC (mu = 1)
sensing_rate_ISAC_mu = zeros(1, num_mu);
sum_ECR_ISAC_mu = zeros(1, num_mu);

% Perform Simulations for Varying mu
for mu_idx = 1:num_mu
    mu = mu_range(mu_idx);         % Current mu value
    
    % Initialize Accumulators
    sum_RN_FDSAC_mu = 0;
    sum_RF_FDSAC_mu = 0;
    sum_RN_ISAC_mu = 0;
    sum_RF_ISAC_mu = 0;
    sum_Rs_FDSAC_mu_val = 0;
    sum_Rs_ISAC_mu_val = 0;
    
    for trial = 1:num_trials
        % Dynamic Power Allocation for ISAC
        alpha_N = alpha_N_fixed + 0.1 * rand();     % alpha_N ∈ [0.2, 0.3]
        alpha_F = alpha_F_fixed - 0.1 * rand();     % alpha_F ∈ [0.5, 0.8]
        
        % Ensure that alpha_N + alpha_F <= 1 for valid power allocation
        if (alpha_N + alpha_F) > 1
            alpha_F = 1 - alpha_N;
        end
        
        % Generate Rayleigh Fading Channels for NU and FU
        hN = sqrt(0.9) * (randn(1,1) + 1i * randn(1,1)) / sqrt(2);  % Channel for NU
        hF = sqrt(0.2) * (randn(1,1) + 1i * randn(1,1)) / sqrt(2);  % Channel for FU
        
        %% ISAC Ergodic Communication Rate Calculation
        % SINR Calculation for NU and FU in ISAC
        gamma_sic_NU_ISAC = (p_fixed * abs(hN)^2 * alpha_F) / (sigma_c^2 + p_fixed * abs(hN)^2 * alpha_N);
        gamma_NU_ISAC = (p_fixed * abs(hN)^2 * alpha_N) / sigma_c^2;
        gamma_FU_ISAC = (p_fixed * abs(hF)^2 * alpha_F) / (sigma_c^2 + p_fixed * abs(hF)^2 * alpha_N);
        
        % Accumulate Ergodic Communication Rates for ISAC
        sum_RN_ISAC_mu = sum_RN_ISAC_mu + log2(1 + gamma_NU_ISAC);
        sum_RF_ISAC_mu = sum_RF_ISAC_mu + log2(1 + gamma_FU_ISAC);
        
        %% FDSAC Ergodic Communication Rate Calculation
        % SINR Calculation for NU and FU in FDSAC
        gamma_NU_FDSAC_mu = (mu * p_fixed * abs(hN)^2) / sigma_c^2;
        gamma_FU_FDSAC_mu = (mu * p_fixed * abs(hF)^2) / sigma_c^2;
        
        % Accumulate Ergodic Communication Rates for FDSAC
        sum_RN_FDSAC_mu = sum_RN_FDSAC_mu + kappa_fixed * log2(1 + gamma_NU_FDSAC_mu);
        sum_RF_FDSAC_mu = sum_RF_FDSAC_mu + kappa_fixed * log2(1 + gamma_FU_FDSAC_mu);
        
        %% Sensing Rate Calculation at p = 5 dB for Current mu
        % Generate Perturbed Correlation Matrix R
        perturbation_std = 0.1;                             % Standard deviation for perturbation
        lambda_perturbed = lambda_mean + perturbation_std * randn(1, num_eigen);  % Perturbed eigenvalues
        lambda_perturbed(lambda_perturbed < 0) = 0;          % Ensure non-negative eigenvalues
        R_random = diag(lambda_perturbed);                   % Correlation matrix
        
        % ISAC Sensing Rate Calculation
        Rs_ISAC_mu = (1 / L) * log2(det(eye(M) + (p_fixed * L / sigma_s^2) * R_random));
        sum_Rs_ISAC_mu_val = sum_Rs_ISAC_mu_val + Rs_ISAC_mu;
        
        % FDSAC Sensing Rate Calculation
        Rs_FDSAC_mu = ((1 - kappa_fixed) / L) * log2(det(eye(M) + ((1 - mu) * p_fixed * L / ((1 - kappa_fixed) * sigma_s^2)) * R_random));
        sum_Rs_FDSAC_mu_val = sum_Rs_FDSAC_mu_val + Rs_FDSAC_mu;
    end
    
    % Compute Average Ergodic Communication Rates
    ECR_NU_FDSAC_mu = sum_RN_FDSAC_mu / num_trials;
    ECR_FU_FDSAC_mu = sum_RF_FDSAC_mu / num_trials;
    sum_ECR_FDSAC_mu(mu_idx) = ECR_NU_FDSAC_mu + ECR_FU_FDSAC_mu;
    
    ECR_NU_ISAC_mu_val = sum_RN_ISAC_mu / num_trials;
    ECR_FU_ISAC_mu_val = sum_RF_ISAC_mu / num_trials;
    sum_ECR_ISAC_mu(mu_idx) = ECR_NU_ISAC_mu_val + ECR_FU_ISAC_mu_val;
    
    % Compute Average Sensing Rates
    sensing_rate_ISAC_mu(mu_idx) = sum_Rs_ISAC_mu_val / num_trials;
    sensing_rate_FDSAC_mu(mu_idx) = sum_Rs_FDSAC_mu_val / num_trials;
end

%% 6. Plot Outage Probability Results
figure;
semilogy(SNR_dB_outage, OP_NU_ISAC, '-o', 'LineWidth', 1.5, 'DisplayName', 'Outage Probability NU (ISAC)');
hold on;
semilogy(SNR_dB_outage, OP_FU_ISAC, '-x', 'LineWidth', 1.5, 'DisplayName', 'Outage Probability FU (ISAC)');
semilogy(SNR_dB_outage, OP_NU_FDSAC, '--o', 'LineWidth', 1.5, 'DisplayName', 'Outage Probability NU (FDSAC)');
semilogy(SNR_dB_outage, OP_FU_FDSAC, '--x', 'LineWidth', 1.5, 'DisplayName', 'Outage Probability FU (FDSAC)');
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Outage Probability', 'FontSize', 12);
legend('Location', 'southwest', 'FontSize', 10);
title('Outage Probability for NU and FU in ISAC and FDSAC', 'FontSize', 14);
grid on;
hold off;

%% 7. Plot Sensing Rate Results
figure;
plot(SNR_dB_sensing, sensing_rate_ISAC_sim, '-s', 'LineWidth', 1.5, 'DisplayName', 'Sensing Rate (ISAC) - Simulated');
hold on;
plot(SNR_dB_sensing, sensing_rate_FDSAC_sim, '--o', 'LineWidth', 1.5, 'DisplayName', 'Sensing Rate (FDSAC) - Simulated');
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Sensing Rate (bps/Hz)', 'FontSize', 12);
legend('Location', 'northwest', 'FontSize', 10);
title('Sensing Rate vs SNR for ISAC and FDSAC', 'FontSize', 14);
grid on;
hold off;

%% 8. Plot Ergodic Communication Rate Results
figure;
plot(SNR_dB_outage, ECR_NU_ISAC, '-^', 'LineWidth', 1.5, 'DisplayName', 'ECR NU (ISAC)');
hold on;
plot(SNR_dB_outage, ECR_FU_ISAC, '-v', 'LineWidth', 1.5, 'DisplayName', 'ECR FU (ISAC)');
plot(SNR_dB_outage, ECR_NU_FDSAC, '--^', 'LineWidth', 1.5, 'DisplayName', 'ECR NU (FDSAC)');
plot(SNR_dB_outage, ECR_FU_FDSAC, '--v', 'LineWidth', 1.5, 'DisplayName', 'ECR FU (FDSAC)');
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Ergodic Communication Rate (bps/Hz)', 'FontSize', 12);
legend('Location', 'northwest', 'FontSize', 10);
title('Ergodic Communication Rate vs SNR for NU and FU in ISAC and FDSAC', 'FontSize', 14);
grid on;
hold off;

%% 9. Plot Sensing Rate vs Sum Communication Rate at p = 5 dB
figure;
plot(sum_ECR_ISAC_mu, sensing_rate_ISAC_mu, '-^', 'LineWidth', 2, 'DisplayName', 'ISAC');
hold on;
plot(sum_ECR_FDSAC_mu, sensing_rate_FDSAC_mu, '--o', 'LineWidth', 2, 'DisplayName', 'FDSAC');
xlabel('Sum Communication Rate (bps/Hz)', 'FontSize', 12);
ylabel('Sensing Rate (bps/Hz)', 'FontSize', 12);
legend('Location', 'northwest', 'FontSize', 10);
title(['Sensing Rate vs Sum Communication Rate at p = ', num2str(p_fixed_dB), ' dB'], 'FontSize', 14);
grid on;
hold off;