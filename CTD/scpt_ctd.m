
%% load

celltypes = table2cell(readtable('celltypes_PSP.csv'));
load('label.mat')  % gene names
load('result.mat') % PLS result
load('genes.mat')  % gene idx

genenames = cellstr(celltypes(:,1));

%% get index of genes with specific cell type expression

for k = 1:length(genenames)
    if ismember(genenames(k),label)
        celltypes(k,3) = num2cell(find(strcmp(label,genenames(k))));
    else 
        celltypes(k,3) = {0};
    end
end
bad_idx = cell2mat(celltypes(:,3))==0;
celltypes(bad_idx,:) = [];
genenames(bad_idx) = [];

% find index of genes in each cell type
[C,~,i] = unique(cellstr(celltypes(:,2)));
ctd_idx = cell(length(C),1);

% number of cell types
ntypes = length(C);

for k = 1:ntypes
    ctd_idx{k} = celltypes(i==k,3);
end

%% find cell type ratio

% get gene sets
gpos_idx = find(result.boot_result.compare_u(:,1) > 2.58);
gneg_idx = find(result.boot_result.compare_u(:,1) < -2.58);
g = genes.scale125.stable;

% find empirical cell type ratio
ctd_ratios = zeros(ntypes,2);
for k = 1:length(ctd_ratios)
    ctd_ratios(k,1) = length(intersect(g(gpos_idx),cell2mat(ctd_idx{k})))/length(gpos_idx);
    ctd_ratios(k,2) = length(intersect(g(gneg_idx),cell2mat(ctd_idx{k})))/length(gneg_idx);
end

%% null model

nperm = 10000;
ctd_null = zeros(ntypes,2,nperm);
for k = 1:nperm
    
    % positive nulls
    y = datasample([1:length(g)],length(gpos_idx),'Replace',false);
    for j = 1:ntypes
        ctd_null(j,1,k) = length(intersect(g(y),cell2mat(ctd_idx{j})))/length(gpos_idx);
    end
    
    % negative nulls
    y = datasample([1:length(g)],length(gneg_idx),'Replace',false);
    for j = 1:ntypes
        ctd_null(j,2,k) = length(intersect(g(y),cell2mat(ctd_idx{j})))/length(gneg_idx);
    end
end

% get p-values with two-tailed permutation test
p_ctd = ctd_ratios - mean(ctd_null,3);
p_null = ctd_null - mean(ctd_null,3);
pval = zeros(ntypes,2);
for k = 1:ntypes
    pval(k,1) = nnz(find(abs(p_null(k,1,:)) >= abs(p_ctd(k,1))))/nperm; % positive genes
    pval(k,2) = nnz(find(abs(p_null(k,2,:)) >= abs(p_ctd(k,2))))/nperm; % negative genes
end

%% visualize

figure;

subplot(1,2,1) % positive genes
scatter(1:ntypes,ctd_ratios(:,1),50,'filled')
hold on
boxplot(squeeze(ctd_null(:,1,:))')
set(gca,'xticklabel',C)
ylim([min(ctd_ratios,[],'all') max(ctd_ratios,[],'all')])
title('pos')

subplot(1,2,2) % negative genes
scatter(1:ntypes,ctd_ratios(:,2),50,'filled')
hold on
boxplot(squeeze(ctd_null(:,2,:))')
set(gca,'xticklabel',C)
ylim([min(ctd_ratios,[],'all') max(ctd_ratios,[],'all')])
title('neg')    