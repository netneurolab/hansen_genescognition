% This script finds the biological processes with which two gene sets are
% most involved. For each biological process (category), the mean loading
% of genes is computed (called the category score). Significance is
% assessed against a null model of category scores. This null model is
% constructed using the loadings of the same genes of interest,
% where null loadings come from PLS analysis on a
% spatial-autocorrelation preserving permutation of one matrix.

% see https://github.com/benfulcher/GeneSetEnrichmentAnalysis for more
% details on GO annotations

%% load
load('gene_expression.mat') % node by gene expression matrix
load('neurosynth.mat')      % node by term probability matrix
load('nodes.mat')           % relevant node indices
load('genes.mat')           % relevant gene indices
load('terms.mat')           % relevant terms
load('label.mat')           % all gene names
load('spins.mat')           % spin test indices

load('GeneSetEnrichmentAnalysis-master\ProcessedData\GOAnnotationDirect-human.mat')
load('GeneSetEnrichmentAnalysis-master\GOTerms_BP.mat')

%% get genes with entrezID

T = table2cell(readtable('gene_entrez_ids')); % load entrezID of genes

gene_name = label(genes.scale125.stable);     % get relevant gene names
entrezIDs = zeros(size(gene_name));

idx = [];
for k = 1:length(gene_name)                                                % for each gene
    if ismember(gene_name(k), T(:,1))                                      % if the gene has an entrezID
        entrezIDs(k) = cell2mat(T(find(strcmp(gene_name(k), T(:,1))),2));  % store the entrezID
        idx = [idx;k];                                                     % also store the index of the gene
    end
end
entrezIDs = entrezIDs(entrezIDs ~= 0);                                     % remove all genes without entrezID


%% get category scores

% set up PLS analysis
X = zscore(expression125);
Y = cogact125;

n = nodes.scale125.lefthem;
g = genes.scale125.stable(idx); % stable genes with entrez ID
t = terms.all;

nnodes = length(n);
nterms = length(t);
ngenes = length(g);

% behav pls
option.method = 3;
option.stacked_behavdata = zscore(Y(n,t));

exp{1} = X(:,g);

result = pls_analysis(exp, nnodes, 1, option); % empirical category scores come from here

% get gene sets
% compute the loading of each gene as the correlation between the original
% data and the gene scores
gload = zeros(ngenes,1);
for k = 1:ngenes
    gload(k) = corr(expression125(:,g(k)),result.usc(:,1));
end

ipos = find(gload > 0); % index of genes with positive loading
ineg = find(gload < 0); % index of genes with negative loading
gload_pos = gload(gload > 0); % loading of genes with positive loading
gload_neg = gload(gload < 0); % loading of genes with negative loading
[~,Ipos] = sort(gload_pos); % sorted
[~,Ineg] = sort(gload_neg); % sorted 

threshold = 0.5; % top 50% of pos/neg genes constitute each gene set

gpos_idx = Ipos(end-floor(threshold*length(gload_pos)):end); % top 50% of genes with positive loading
gneg_idx = Ineg(1:floor(threshold*length(gload_neg)));       % top 50% of genes with negative loading

gpos_ID = entrezIDs(ipos(gpos_idx));   % these are the entrezIDs of the genes in the positive set
gneg_ID = entrezIDs(ineg(gneg_idx));  % these are the entrezIDs of the genes in the negative set

% get category scores
categoryScores = zeros(length(allGOCategories),2);
for k = 1:length(categoryScores)           % for each category
    
    GOcategory = geneEntrezAnnotations{k}; % get genes in this category
    
    category_pos = [];
    category_neg = [];
    
    for j = 1:length(GOcategory)           % for each gene in the category
        
        % check if the gene is in the positive/negative gene set
        % if so, save the loading
        if ismember(GOcategory(j),gpos_ID)
            category_pos = [category_pos; gload(find(entrezIDs == GOcategory(j)))];
        elseif ismember(GOcategory(j),gneg_ID)
            category_neg = [category_neg; gload(find(entrezIDs == GOcategory(j)))];
        end
        
    end
    
    % category score = mean loading of relevant genes in category
    categoryScores(k,1) = mean(category_pos);
    categoryScores(k,2) = mean(category_neg);   
end

%% get null category scores

nperm = 10000;
categoryScores_null = zeros(length(allGOCategories),2,nperm);

for m = 1:nperm % for each permutation
    
    % PLS analysis
    option.stacked_behavdata = zscore(Y(scale125spins(:,m),t)); % rows of neurosynth matrix has been permuted while preserving spatial autocorrelation
    result = pls_analysis(exp, nnodes, 1, option);              % this is the null PLS result
    
    gload = zeros(ngenes,1);                                    % these are the null gene loadings
    for k = 1:ngenes
        gload(k) = corr(expression125(:,g(k)),result.usc(:,1));
    end
    
    for k = 1:length(allGOCategories)          % for each gene in category
        
        GOcategory = geneEntrezAnnotations{k}; % get genes in the category
        
        category_pos = [];
        category_neg = [];
        
        for j = 1:length(GOcategory)           % for each gene in the category
            
            % check if the gene is in the positive/negative gene set
            % if so, save the loading
            if ismember(GOcategory(j),gpos_ID)
                category_pos = [category_pos; gload(find(entrezIDs == GOcategory(j)))];
            elseif ismember(GOcategory(j),gneg_ID)
                category_neg = [category_neg; gload(find(entrezIDs == GOcategory(j)))];
            end
        end
        
        % category score = mean loading of relevant genes in category        
        categoryScores_null(k,1,m) = mean(category_pos);
        categoryScores_null(k,2,m) = mean(category_neg);
    end
end

%% get biological processes

% get categories with a score
pos_cat = find(~isnan(categoryScores(:,1)));
neg_cat = find(~isnan(categoryScores(:,2)));

% set up p-value template
pos_pvals = zeros(length(pos_cat),1);
neg_pvals = zeros(length(neg_cat),1);

% calculate permuted p-value
for k = 1:length(pos_pvals)
    pos_pvals(k) = (1+(nnz(find(categoryScores_null(pos_cat(k),1,:) > categoryScores(pos_cat(k),1)))))/(nperm+1);
end
for k = 1:length(neg_pvals)
    neg_pvals(k) = (1+(nnz(find(categoryScores_null(neg_cat(k),2,:) < categoryScores(neg_cat(k),2)))))/(nperm+1);
end

% get biological processeses
GO_categories = table2cell(GOTable(:,2));    % these are all the cateogy names
GO_categoriesID = table2array(GOTable(:,3)); % this is the ID of the category

positive_terms = cell(length(pos_cat),1);
negative_terms = cell(length(neg_cat),1);

% get biological process category if it exists and is available

% for the positive gene set:
for k = 1:length(pos_cat)
    positive_terms{k} = GO_categories(GO_categoriesID == allGOCategories(pos_cat(k)));
end
p_idx = find(~cellfun('isempty',positive_terms));
positive_terms = positive_terms(p_idx);           % remove empty indices

% for the negative gene set:
for k = 1:length(neg_cat)
    negative_terms{k} = GO_categories(GO_categoriesID == allGOCategories(neg_cat(k)));
end
n_idx = find(~cellfun('isempty',negative_terms));
negative_terms = negative_terms(n_idx);           % remove empty indices