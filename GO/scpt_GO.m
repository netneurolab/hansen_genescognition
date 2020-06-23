% see https://github.com/benfulcher/GeneSetEnrichmentAnalysis for more
% details on GO annotations

%% load
load('gene_expression.mat')
load('neurosynth.mat')
load('nodes.mat')
load('genes.mat')
load('terms.mat')
load('label.mat')
load('spins.mat')

% from https://figshare.com/s/71fe1d9b2386ec05f421
load('GeneSetEnrichmentAnalysis-master\ProcessedData\GOAnnotationDirect-human.mat')
load('GeneSetEnrichmentAnalysis-master\GOTerms_BP.mat')

%% get genes with entrezID

T = table2cell(readtable('gene_entrez_ids'));

gene_name = label(genes.scale125.stable);
entrezIDs = zeros(size(gene_name));

idx = [];
for k = 1:length(gene_name)
    if ismember(gene_name(k), T(:,1))
        entrezIDs(k) = cell2mat(T(find(strcmp(gene_name(k), T(:,1))),2));
        idx = [idx;k];
    end
end
entrezIDs = entrezIDs(entrezIDs ~= 0);


%% get category scores

% set up PLS analysis
X = zscore(expression125);
Y = cogact125;

n = nodes.scale125.lefthem;
g = genes.scale125.stable(idx);
t = terms.all;

nnodes = length(n);
nterms = length(t);
ngenes = length(g);

% behav pls
option.method = 3;
option.num_boot = 1000;
option.num_perm = 0;
option.stacked_behavdata = zscore(Y(n,t));

exp{1} = X(:,g);

result = pls_analysis(exp, nnodes, 1, option);

% get positive/negative gene sets
btsp = result.boot_result.compare_u(:,1);
gpos_btsp = btsp(btsp > 2.58); % positive gene set
gneg_btsp = btsp(btsp < -2.58); % negative gene set

gpos_ID = entrezIDs(find(btsp > 2.58));
gneg_ID = entrezIDs(find(btsp < -2.58));

% get category scores
categoryScores = zeros(length(allGOCategories),2);
for k = 1:length(categoryScores) % for each category
    
    GOcategory = geneEntrezAnnotations{k}; % get genes in this category
    
    category_pos = [];
    category_neg = [];
    
    for j = 1:length(GOcategory) % for each gene in the category
        
        % check if the gene is in the positive/negative gene set
        % if so, save the bootstrap ratio
        if ismember(GOcategory(j),gpos_ID)
            category_pos = [category_pos; btsp(find(entrezIDs == GOcategory(j)))];
        elseif ismember(GOcategory(j),gneg_ID)
            category_neg = [category_neg; btsp(find(entrezIDs == GOcategory(j)))];
        end
        
    end
    
    % category score = mean bootstrap ratio of relevant genes in category
    categoryScores(k,1) = mean(category_pos);
    categoryScores(k,2) = mean(category_neg);
    
end

%% get null category scores

nperm = 1000;
categoryScores_null = zeros(length(allGOCategories),2,nperm);

for m = 1:nperm
    
    % behav PLS with spatially permuted rows
    option.stacked_behavdata = zscore(Y(scale125spins(:,m),t));
    result = pls_analysis(exp, nnodes, 1, option);
    
    btsp = result.boot_result.compare_u(:,1);
    
    for k = 1:length(allGOCategories) % for each gene in category
        
        GOcategory = geneEntrezAnnotations{k}; % get genes in the category
        
        category_pos = [];
        category_neg = [];
        
        for j = 1:length(GOcategory) % for each gene in the category
            
            % check if the gene is in the positive/negative gene set
            % if so, save the bootstrap ratio
            if ismember(GOcategory(j),gpos_ID)
                category_pos = [category_pos; btsp(find(entrezIDs == GOcategory(j)))];
            elseif ismember(GOcategory(j),gneg_ID)
                category_neg = [category_neg; btsp(find(entrezIDs == GOcategory(j)))];
            end
        end
        
        % category score = mean bootstrap ratio of relevant genes in category        
        categoryScores_null(k,1,m) = mean(category_pos);
        categoryScores_null(k,2,m) = mean(category_neg);
    end
end

%% get biological processes

% get categories with a score
pos_cat = find(~isnan(categoryScores(:,1)));
neg_cat = find(~isnan(categoryScores(:,2)));

pos_pvals = zeros(length(pos_cat),1);
neg_pvals = zeros(length(neg_cat),1);

% calculate p-value
for k = 1:length(pos_pvals)
    pos_pvals(k) = nnz(find(categoryScores_null(pos_cat(k),1,:) > categoryScores(pos_cat(k),1)))/nperm;
end
for k = 1:length(neg_pvals)
    neg_pvals(k) = nnz(find(categoryScores_null(neg_cat(k),2,:) < categoryScores(neg_cat(k),2)))/nperm;
end

% get biological processeses
GO_categories = table2cell(GOTable(:,2));
GO_categoriesID = table2array(GOTable(:,3));

positive_terms = cell(length(pos_cat),1);
negative_terms = cell(length(neg_cat),1);

% get available biological processes for each category
for k = 1:length(pos_cat)
    positive_terms{k} = GO_categories(GO_categoriesID == allGOCategories(pos_cat(k)));
end
p_idx = find(~cellfun('isempty',positive_terms));
positive_terms = positive_terms(p_idx);

for k = 1:length(neg_cat)
    negative_terms{k} = GO_categories(GO_categoriesID == allGOCategories(neg_cat(k)));
end
n_idx = find(~cellfun('isempty',negative_terms));
negative_terms = negative_terms(n_idx);

%% make table

pGeneCats = cell(length(p_idx),3);
nGeneCats = cell(length(n_idx),3);

% positive biological processes
pGeneCats(:,1) = positive_terms(p_idx);
pGeneCats(:,2) = num2cell(categoryScores(pos_cat(p_idx),1));
pGeneCats(:,3) = num2cell(pos_pvals(p_idx));
[~,I] = sort(pos_pvals(p_idx));
pGeneCats = pGeneCats(I,:);
pGeneCats = cell2table(pGeneCats);
pGeneCats.Properties.VariableNames{'pGeneCats1'} = 'Biological Process';
pGeneCats.Properties.VariableNames{'pGeneCats2'} = 'Category Score';
pGeneCats.Properties.VariableNames{'pGeneCats3'} = 'p-value';
writetable(pGeneCats,'positive_gene_categories.csv')

% negative biological processes
nGeneCats(:,1) = negative_terms(n_idx);
nGeneCats(:,2) = num2cell(categoryScores(neg_cat(n_idx),2));
nGeneCats(:,3) = num2cell(neg_pvals(n_idx));
[~,I] = sort(neg_pvals(n_idx));
nGeneCats = nGeneCats(I,:);
nGeneCats = cell2table(nGeneCats);
nGeneCats.Properties.VariableNames{'nGeneCats1'} = 'Biological Process';
nGeneCats.Properties.VariableNames{'nGeneCats2'} = 'Category Score';
nGeneCats.Properties.VariableNames{'nGeneCats3'} = 'p-value';
writetable(nGeneCats,'negative_gene_categories.csv')