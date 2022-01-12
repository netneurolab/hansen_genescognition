% This script calculates the ratio of genes in two gene sets that are
% preferentially expressed in 7 cell types:
% astrocytes (astro), microglia (micro), endothelial cells (endo), 
% excitatory neurons (neuro-ex), inhibitory neurons (neuro-in), 
% oligodendrocyte precursors (opc), and oligodendrocytes (oligo).
% Significance is assessed using a null model constructed from random gene
% sets.
%% load

celltypes = table2cell(readtable('celltypes_PSP.csv')); % load specific cell type expression
load('label.mat')           % gene names
load('result.mat')          % PLS result from scpt_genes_cog_pls.m
load('genes.mat')           % gene info
load('gene_expression.mat') % node x gene expression matrix
g = genes.scale125.stable;  % index of genes

genenames = cellstr(celltypes(:,1));

%% get index of genes with specific cell type expression

for k = 1:length(genenames)                                           % for each gene with specific cell type expression 
    if ismember(genenames(k),label)                                   % if gene overlaps with AHBA genes
        celltypes(k,3) = num2cell(find(strcmp(label,genenames(k))));  % add index of gene
    else 
        celltypes(k,3) = {0};                                         % otherwise, add 0
    end
end
bad_idx = cell2mat(celltypes(:,3))==0; % find indices of genes not in AHBA
celltypes(bad_idx,:) = [];             % remove these genes

% find index of genes in each cell type
[C,~,i] = unique(cellstr(celltypes(:,2))); % index genes by which cell type (in C) they're expressed in
ctd_idx = cell(length(C),1);               % ctd_idx will store the gene idx associated to each cell type

% number of cell types
ntypes = length(C);

for k = 1:ntypes                    % for each cell type
    ctd_idx{k} = celltypes(i==k,3); % store gene indices related to the cell type
end

%% find cell type ratio

% get gene sets
% compute the loading of each gene as the correlation between the original
% data and the gene scores
gload = corr(expression125(:,g),result.usc(:,1));

ipos = find(gload > 0); % index of genes with positive loading
ineg = find(gload < 0); % index of genes with negative loading
gload_pos = gload(gload > 0); % loading of genes with positive loading
gload_neg = gload(gload < 0); % loading of genes with negative loading
[~,Ipos] = sort(gload_pos); % sorted
[~,Ineg] = sort(gload_neg); % sorted 

threshold = 0.5; % top 50% of pos/neg genes constitute each gene set

gpos_idx = Ipos(end-floor(threshold*length(gload_pos)):end); % top 50% of genes with positive loading
gneg_idx = Ineg(1:floor(threshold*length(gload_neg)));       % top 50% of genes with negative loading

% find empirical cell type ratio
ctd_ratios = zeros(ntypes,2);
for k = 1:length(ctd_ratios)
    ctd_ratios(k,1) = length(intersect(g(ipos(gpos_idx)),cell2mat(ctd_idx{k})))/length(gpos_idx);
    ctd_ratios(k,2) = length(intersect(g(ineg(gneg_idx)),cell2mat(ctd_idx{k})))/length(gneg_idx);
end

%% null model from random gene set

n = 10000;
ctd_null = zeros(ntypes,2,n);
for k = 1:n % for each repetition
    
    % positive nulls
    y = datasample([1:length(g)],length(gpos_idx),'Replace',false);                      % get random gene set the size of the positive gene set
    for j = 1:ntypes                                                                     % for each cell type
        ctd_null(j,1,k) = length(intersect(g(y),cell2mat(ctd_idx{j})))/length(gpos_idx); % find ratio of genes expressed in cell type to all genes
    end
    
    % negative nulls
    y = datasample([1:length(g)],length(gneg_idx),'Replace',false);                      % get random gene set the size of the negative gene set
    for j = 1:ntypes                                                                     % for each cell type
        ctd_null(j,2,k) = length(intersect(g(y),cell2mat(ctd_idx{j})))/length(gneg_idx); % find ratio of genes expressed in cell type to all genes 
    end
end

% get p-values with two-tailed significance test
p_ctd = ctd_ratios - mean(ctd_null,3); % mean centre
p_null = ctd_null - mean(ctd_null,3);
pval = zeros(ntypes,2);
for k = 1:ntypes
    pval(k,1) = (1+(nnz(find(abs(p_null(k,1,:)) >= abs(p_ctd(k,1))))))/(n+1); % pval for positive gene set
    pval(k,2) = (1+(nnz(find(abs(p_null(k,2,:)) >= abs(p_ctd(k,2))))))/(n+1); % pval for negative gene set
end

for k = 1:2
    pval(:,k) = mafdr(pval(:,k),'BHFDR',true); % FDR correction
end

%% visualize
% plot empirical ratios as points and null models as boxplots

o = [1,3,6,7,2,4,5]; % order in which cell types appear

figure;

subplot(1,2,1) % specific cell type expression for positive gene set
scatter(1:ntypes,ctd_ratios(o,1),30,'filled')
hold on
boxplot(squeeze(ctd_null(o,1,:))')
set(gca,'xticklabel',C(o))
xtickangle(90)
% ylim([min(ctd_ratios,[],'all') max(ctd_ratios,[],'all')])
ylim([-0.01 0.12])
title('neg')

subplot(1,2,2) % specific cell type expression for negative gene set
scatter(1:ntypes,ctd_ratios(o,2),30,'filled')
hold on
boxplot(squeeze(ctd_null(o,2,:))')
set(gca,'xticklabel',C(o))
xtickangle(90)
ylim([-0.01 0.12])
% ylim([min(ctd_ratios,[],'all') max(ctd_ratios,[],'all')])
title('pos')    