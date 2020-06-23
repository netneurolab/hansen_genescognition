%% load

load('hcp_smyl_all_125.mat') % t1w/t2w ratios
load('hcp_thi_all_125.mat')  % cortical thickness
load('behav_idx.mat')        % behavioural measures
load('subject_idx.mat')      % relevant subjects
load('result.mat')           % PLS result

%% get individual correlations

% individual manifestation of gene score pattern via cortical thickness or
% T1w/T2w ratio
subject_corrs = zeros(size(hcp_smyl,1),2);
subject_corrs(:,1) = corr(hcp_thi(:,109:end)',result.usc(:,1));
subject_corrs(:,2) = corr(hcp_smyl(:,109:end)', result.usc(:,1));

% correlate subjects' manifestation against behavioural performance
y = table2array(HCP_unrestricted(subject_idx, behav_idx));

behav_corr = zeros(size(y,2),2);
for j = 1:2
    for k = 1:size(behav_corr,1)
        behav_corr(k,j) = corr(subject_corrs(:,j),y(:,k),'rows','complete');
    end
end

%% permuted null model

nperm = 10000;
null = zeros(nperm,length(behav_idx),2);

for j = 1:2 % get permuted null models
    x = subject_corrs(:,j);
    for i = 1:size(y,2)
        for k = 1:nperm
            idx = randperm(length(subject_idx));
            null(k,i,j) = corr(x(idx),y(:,i),'rows','complete');
        end
    end
end

% two-tailed fdr corrected significance test
pval = zeros(size(y,2),2); 
behav_corr = behav_corr - squeeze(mean(null,1));
null = null - mean(null,1);
for j = 1:2
    for i = 1:size(y,2)
        pval(i,j) = nnz(find(abs(null(:,i,j)) >= abs(behav_corr(i,j))))/nperm;
    end
end

for j = 1:2
    pval(:,j) = mafdr(pval(:,j),'BHFDR',true);
end

%% visualize

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