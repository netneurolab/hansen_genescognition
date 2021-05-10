% This script relates the PLS-derived score pattern to individual
% differences in behaviour using the Human Connectome Project dataset. 
% This section got cut in the review process, but I left it because it's
% in the preprint.

%% load

load('hcp_smyl_all_125.mat') % t1w/t2w ratios
load('hcp_thi_all_125.mat')  % cortical thickness
load('behav_idx.mat')        % indeces of relevant behavioural measures
load('subject_idx.mat')      % indeces of subjects
load('result.mat')           % PLS result (from scpt_genes_cog_pls.m)

%% get individual correlations

% individual manifestation of gene score pattern via cortical thickness or
% T1w/T2w ratio
subject_corrs = zeros(size(hcp_smyl,1),2);                        % this is a subjects x 2 array of correlations
subject_corrs(:,1) = corr(hcp_thi(:,109:end)',result.usc(:,1));   % correlate cortical thickness to gene scores for each subject
subject_corrs(:,2) = corr(hcp_smyl(:,109:end)', result.usc(:,1)); % correlate T1w/T2w ratio to gene scores for each subject

% correlate subjects' manifestation against behavioural performance
y = table2array(HCP_unrestricted(subject_idx, behav_idx)); % get behavioural performance for the behavioural measures in behav_idx and the subjects in subject_idx

behav_corr = zeros(size(y,2),2);
for j = 1:2                                                                  % for cortical thickness and T1w/T2w ratios
    for k = 1:size(behav_corr,1)                                             % for each behavioural measure
        behav_corr(k,j) = corr(subject_corrs(:,j),y(:,k),'rows','complete'); % correlate subjects' correlations to behavioural performance
    end
end

%% permuted null model

nperm = 10000;
null = zeros(nperm,length(behav_idx),2);

% get permuted null models
for j = 1:2                                                                % for cortical thickness and T1w/T2w ratio         
    x = subject_corrs(:,j);                                                % get the vector of physical brain measure to gene score correlation
    for i = 1:size(y,2)                                                    % for each behavioural measure
        for k = 1:nperm                                                    % for each permutation
            idx = randperm(length(subject_idx));                           % permute subjects randomly
            null(k,i,j) = corr(x(idx),y(:,i),'rows','complete');           % correlate permuted subjects with nonpermuted behavioural performance
        end
    end
end

% two-tailed fdr corrected significance test
pval = zeros(size(y,2),2); 
behav_corr = behav_corr - squeeze(mean(null,1)); % mean centre
null = null - mean(null,1);                      % mean centre
for j = 1:2
    for i = 1:size(y,2)
        pval(i,j) = (1+(nnz(find(abs(null(:,i,j)) >= abs(behav_corr(i,j))))))/(nperm+1);
    end
end

for j = 1:2
    pval(:,j) = mafdr(pval(:,j),'BHFDR',true); % FDR correction
end

%% visualize
% this plots the behavioural measures that significantly correlate with
% subject_corrs (manifestation of PLS-derived gene score via cortical
% thickness or T1w/T2w ratio). Null models are shown as boxplots and
% empirical values as points.

t = {'gene score - cortical thickness','gene score - T1w/T2w ratio'};

figure;
for k = 1:2
    subplot(2,2,k)
    sig = find(fdr(:,k) < 0.05);
    boxplot(null_behav_corr(:,sig,k))
    hold on;
    scatter([1:length(sig)],behav_corr(sig,k),20,'filled')
    xticklabels(HCPmeasures(behav_idx(sig)))
    xtickangle(90)
    ylim([-0.15,0.16])
    title(t{k})
end