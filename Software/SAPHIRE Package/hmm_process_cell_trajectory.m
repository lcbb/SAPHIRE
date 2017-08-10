%%%%%%%%%%%%%%%%%%%%
% Hidden Markov modeling and state annotation of cellular image time series
% with Bayesian modeling selection.
%
% Inputs:   traj - a dxT matrix, where d are the independent variables that
%           comprise the trajectory and T is the trajectory length.
%           Kmax - maximum number of hidden states to test
%           mcmc_params.parallel = 'on' or 'off' for parallel model
%           procesing         
% Outputs:  PrM - model probabilities
%           ML_states - Maximum likelihood state annotations
%           ML_params - Maximum likelihood model parameters
%           full_results - Results for all models tested
%           full_fitting - MCMC samples, used for error calculations
%           logI - log-likelihoods of models tested
%               
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function [PrM, ML_states, ML_params, full_results, full_fitting, logI] = hmm_process_cell_trajectory(traj,Kmax,mcmc_params)

% Ensure that, if there are multiple trajectories, has dimensions 1xn where
% n is the number of trajectories.
if iscell(traj)
    traj = reshape(traj,1,numel(traj));
end

% Choose ranges for MCMC parameter initializations
if iscell(traj)
    traj_matrix = cell2mat(traj);
else
    traj_matrix = traj;
end
mu_max = max(traj_matrix,[],2); 
mu_min = min(traj_matrix,[],2); 
mu_range = [mu_min mu_max];
sigma_max = max(mu_max-mu_min); 

%Set number of states in models
nModels = Kmax;
Klist = cell(0);
states_list = cell(0);
for i = 1:nModels
    Klist{i} = i;
    states_list{i} = ones(1,i); 
end

logI = zeros(1,nModels);
full_results = struct('PrM',[],'K',Klist,'states',states_list,'ML_states',[],'ML_params',[]);
full_fitting = struct('mcmc_params_final',[],'logI',cell(1,nModels),'logprobs',[],'samples',[]);

%Set additional MCMC parameters
if ~isfield(mcmc_params,'proposaltype')
    mcmc_params.proposaltype = 'gaussian';
end
if ~isfield(mcmc_params,'move')
    mcmc_params.move = 'block';
end

% Find the marginal likelihood of the data for each number of states
if isfield(mcmc_params,'parallel') && strcmp(mcmc_params.parallel,'on')
    
    parpool()
    
    parfor m = 1:nModels
        
        K = Klist{m};
        full_fitting(m).mcmc_params_final = mcmc_params;
        full_fitting(m).mcmc_params_final.states =states_list{m};

        % Run initializations with different starting parameters
        [best_sample, full_fitting(m).mcmc_params_final] = hmm_mcmc_initialization(traj,K,mu_range,sigma_max,full_fitting(m).mcmc_params_final);

        % Run a longer MCMC trajectory with the best starting parameters
        full_fitting(m).mcmc_params_final.nIter = 100000;
        [full_fitting(m).samples, full_fitting(m).logprobs] = hmm_mcmc(traj,best_sample.p_start,best_sample.p_trans,best_sample.mu_emit,best_sample.sigma_emit,full_fitting(m).mcmc_params_final);

        % Perform integration to get marginal likelihood
        nIntegration = 200000;
        logI(m) = hmm_integration_gaussian(traj,full_fitting(m).samples(ceil(full_fitting(m).mcmc_params_final.nIter/2):end),nIntegration,full_fitting(m).mcmc_params_final.states);
        full_fitting(m).logI = logI(m);
        
    end
    
    %Close parallel pool
    delete(gcp)
    
else
    
    for m = 1:nModels

        K = Klist{m};
        full_fitting(m).mcmc_params_final = mcmc_params;
        full_fitting(m).mcmc_params_final.states = states_list{m};

        % Run initializations with different starting parameters
        [best_sample, full_fitting(m).mcmc_params_final] = hmm_mcmc_initialization(traj,K,mu_range,sigma_max,full_fitting(m).mcmc_params_final);

        % Run a longer MCMC trajectory with the best starting parameters
        full_fitting(m).mcmc_params_final.nIter = 100000;
        [full_fitting(m).samples, full_fitting(m).logprobs] = hmm_mcmc(traj,best_sample.p_start,best_sample.p_trans,best_sample.mu_emit,best_sample.sigma_emit,full_fitting(m).mcmc_params_final);

        % Perform integration to get marginal likelihood
        nIntegration = 200000;
        logI(m) = hmm_integration_gaussian(traj,full_fitting(m).samples(ceil(full_fitting(m).mcmc_params_final.nIter/2):end),nIntegration,full_fitting(m).mcmc_params_final.states);
        full_fitting(m).logI = logI(m);
        
    end
end

%Calculate model probabilities
PrM = exp(logI-max(logI))/sum(exp(logI-max(logI)));  %subtract max to avoid underflow

for m=1:nModels
    
    full_results(m).PrM = PrM(m);
    
    % Maximum likelihood parameters for each model
    [~, ML_idx] = max(full_fitting(m).logprobs);
    full_results(m).ML_params = full_fitting(m).samples(ML_idx);
    
    % Uncertainty in ML parameters
    full_results(m).ML_params_error = hmm_param_errors(full_fitting(m).samples(ceil(full_fitting(m).mcmc_params_final.nIter/2):end),full_fitting(m).mcmc_params_final.states);
    
    % Viterbi algorithm to get most probable state sequence
    if iscell(traj)
        full_results(m).ML_states = cell(size(traj));
        for i=1:length(traj)
            full_results(m).ML_states{i} = hmm_viterbi(traj{i},full_results(m).ML_params.p_start,full_results(m).ML_params.p_trans,full_results(m).ML_params.mu_emit,full_results(m).ML_params.sigma_emit);
        end
    else
        full_results(m).ML_states = hmm_viterbi(traj,full_results(m).ML_params.p_start,full_results(m).ML_params.p_trans,full_results(m).ML_params.mu_emit,full_results(m).ML_params.sigma_emit);
    end
    
    % Order states from lowest to largest distance from the origin
    [full_results(m).ML_states, full_results(m).ML_params] = order_states(full_results(m).ML_states, full_results(m).ML_params);
    
end

% Parameters and states of the best model
[~, ML_M] = max(PrM);
ML_states = full_results(ML_M).ML_states;
ML_params = full_results(ML_M).ML_params;

end

