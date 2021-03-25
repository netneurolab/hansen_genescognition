% This script organizes and parcellates data from BrainSpan, then computes
% estimated gene scores based on the gene weights computed in the original
% analyses (see scpt_genes_cog_pls.m). Note that the figure generation will
% run into an error if you don't have cbrewer added to your path
% (https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab)

% BrainSpan data can be downloaded in its original form from
% https://www.brainspan.org/static/download.html

%% load

load('genes.mat')    % gene idx
load('label.mat')    % gene names
load('result34.mat') % PLS result from scpt_genes_cog_pls.m on 34-node parcellation
load('mapping.mat')  % mapping from 34 node parcellation to the 16 unique cortical regions included in BrainSpan


% AHBA harmonized files have been organized to include only genes included
% in AHBA
brainspan = readtable('gene_expression_AHBA_harmonized.csv'); % gene by sample matrix of gene expression
gene_info = readtable('gene_metadata_AHBA_harmonized.csv');   % metadata on genes
sample_info = readtable('samples_metadata.csv');              % metadata on tissue samples
brainspan = table2array(brainspan(2:end,2:end));              % remove column of gene names and row of sample IDs

%% remove non-cortical samples
% remove samples labeled 'not cortex' and amygdala samples

notcortex_idx = [];
sensorifugal = table2cell(sample_info(:,15)); % this column labels noncortical regions as 'Not_Cortex' 
sensorifugal = string(sensorifugal);

% find noncortical defined by sensorifugal
for k = 1:length(sensorifugal)
    if strcmp(sensorifugal(k),'Not_Cortex')
        notcortex_idx = [notcortex_idx; k];
    end
end

% find amygdala indices
amygdala = find(contains(table2cell(sample_info(:,8)),'amygdaloid complex'));
notcortex_idx = [notcortex_idx; amygdala];
notcortex_idx = sort(notcortex_idx);

% get all relevant (cortical) indices
cortex_idx = setdiff([1:length(sensorifugal)],notcortex_idx); 

% remove noncortical indices
brainspan(:,notcortex_idx) = [];

%% remove unstable genes

% genes included in original analyses are stable (differential stability >
% 0.1). For comparability, we only include (available) stable genes as defined on the
% 34-node parcellation.

stable_label = label(genes.scale033.stable); % names of stable genes as defined on 34-node parcellation
gene_label = table2cell(gene_info(:,4));     % names of genes in BrainSpan

notstable_idx = [];
bspan_gidx033 = [];
for k = 1:length(gene_label)                            % for each gene in BrainSpan
    if ~ismember(stable_label,gene_label(k))            % if gene isn't stable
        notstable_idx = [notstable_idx; k];             % keep its index
    else                                                % if gene is stable
        i = find(ismember(stable_label,gene_label(k))); % get index of stable gene
        bspan_gidx033 = [bspan_gidx033; i];             % store it
    end
end

% remove nonstable genes
brainspan(notstable_idx,:) = [];

%% organize data by life stage
% samples are organized into 5 life stages: fetal, infant, child,
% adolescent, and adult

% get indices
fetal_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'fetal'));
infant_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'infant'));
child_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'child'));
adolescent_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'adolescent'));
adult_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'adult'));

% make gene expression matrices for each life stage
fetal = brainspan(:,fetal_idx);
infant = brainspan(:,infant_idx);
child = brainspan(:,child_idx);
adolescent = brainspan(:,adolescent_idx);
adult = brainspan(:,adult_idx);

% organize matrices and indices
M_idx = {fetal_idx, infant_idx, child_idx, adolescent_idx, adult_idx};
M = {fetal, infant, child, adolescent, adult};

%% parcellate brain regions

sample_regions = table2cell(sample_info(cortex_idx,8)); % get region name of each sample
regions = unique(string(sample_regions));               % get unique regions - this is the maximum number of brain regions in each gene expression matrix

% make value-based region mapping instead of string-based (for
% simplicity)
region_mapping = zeros(length(sample_regions),1);

for k = 1:length(regions)                               % for each unique region
    i = find(contains(sample_regions,regions(k)));      % find all samples from that region
    region_mapping(i) = k;                              % make index-based map of regions
end

% average expression of each gene from identical regions
% done separately for each life stage, where the number of unique regions
% with gene expression estimates varies across life stage
regionsIncluded = cell(5,1);
for k = 1:length(M)                                          % for each life stage
    mat = M{k};                                              % get original gene expression matrix (genes x samples)
    r = region_mapping(M_idx{k});                            % get region mapping of this life stage
    regionsIncluded{k} = unique(r);                          % this is how many regions have gene expression estimates
    mat_tmp = zeros(size(mat,1),length(regionsIncluded{k})); % make parcellated matrix template
    i = 1;
    for j = regionsIncluded{k}'                              % for each region that has gene expression estimates
        aRegion = mat(:,find(r==j));                         % find all columns corresponding to same region
        mat_tmp(:,i) = sum(aRegion,2) ./ sum(aRegion~=0,2);  % fill parcellated matrix with mean exp ignoring 0 values
        i = i+1;
    end
    M{k} = mat_tmp';                                         % reassign with new parcellated matrix (genes x regions)
end

%% get gene and term scores

% estimate gene scores by multiplying PLS-derived gene weights with gene
% expression matrices

gscore = cell(length(M),1);

% for k = 1:length(M)                             % for each life stage
%     m = M{k};                                   % get gene expression matrix (genes x brain regions)
%     med = repmat(nanmedian(m,2), 1, size(m,2)); % fill missing data with median expression of gene across regions
%     m(isnan(m)) = med(isnan(m));
%     gscore{k} = m*result.u(bspan_gidx033);      % estimate gene score
% end

% as per reviewer's request, use only genes with complete data instead of
% imputing with median expression (commented above, method for the preprint)
missing_genes = []; 
for k = 1:length(M) % for each life stage
    m = M{k}; % get the gene expression matrix (genes x brain regions)
    [~,i] = find(isnan(m)); % get 
    missing_genes = union(missing_genes, i);
end

bspan_gidx033(missing_genes) = [];

for k = 1:length(M)
    m = M{k};
    m(:,missing_genes) = [];
    M{k} = m;
    gscore{k} = m * result.u(bspan_gidx033);
end

% convert PLS-derived term scores into 16-node parcellation, averaging
% sibling nodes
tscore = zeros(length(regions),1);
for k = 1:length(regions)                       % for each unique region
    tscore(k) = mean(result.vsc(mapping==k,1)); % average term scores (from result.vsc) of all sibling nodes corresponding to current region
end

%% track gene scores across development
% see how gene scores at the brain regions with gene expression estimates
% across all 5 life stages change with development

% get regions with available gene expression estimates across all life stages
reg = intersect(regionsIncluded{1},[regionsIncluded{2};regionsIncluded{3};regionsIncluded{4};regionsIncluded{5}]);
tmp = gscore;
tmp{1} = tmp{1}(reg); % only need to change first life stage (16 regions)
gscore_mat = [tmp{1} tmp{2} tmp{3} tmp{4} tmp{5}]; % organize gene scores

%% visualize

% change the colourmap if you don't want to download cbrewer
cm=cbrewer('qual', 'Paired', 16, 'PCHIP');

% gene scores across development

figure;
for k = 1:length(gscore_mat)                                    % for each brain region
    hold on
    plot(gscore_mat(k,:),'LineWidth',1.3,'Color',cm(reg(k),:))  % plot a curve of gene score in that brain region across all five life stages
end
legend(regions(reg),'Location','northwest');
xticks(1:5)
xticklabels({'fetus','infant','child','adolescent','adult'})
ylabel('estimated gene scores')

% gene-term score correlation across development

figure;
for k = 1:length(gscore)                                         % for each life stage
    subplot(2,3,k)
    t = tscore;                                                  % get term scores
    g = gscore{k};                                               % get gene scores
    notIncluded = setdiff(1:length(regions),regionsIncluded{k}); % get regions without a gene score (varies across life stage)
    if size(notIncluded,2) > 0                                   % if there are regions without a gene score
        t(notIncluded) = [];                                     % remove the corresponding term scores
    end
    [rho,pval] = corr(g,t);                                      % correlate gene and term scores in this life stage
    scatter(g,t,50,regionsIncluded{k},'filled')                  % plot and colour each point by brain region
    xlabel('gene score')
    ylabel('term score')
    title([' rho=', num2str(rho), ' p=', num2str(pval)])
end
colormap(cm);
