
%% load

load('genes.mat')
load('label.mat')
load('result.mat')
load('mapping.mat')

brainspan = readtable('gene_expression_AHBA_harmonized.csv');
gene_info = readtable('gene_metadata_AHBA_harmonized.csv');
sample_info = readtable('samples_metadata.csv');

brainspan = table2array(brainspan(2:end,2:end));

%% remove non-cortical samples

% find noncortical indices including amygdala
notcortex_idx = [];
sensorifugal = table2cell(sample_info(:,15));
sensorifugal = string(sensorifugal);

for k = 1:length(sensorifugal)
    if strcmp(sensorifugal(k),'Not_Cortex')
        notcortex_idx = [notcortex_idx; k];
    end
end

amygdala = find(contains(table2cell(sample_info(:,8)),'amygdaloid complex'));
notcortex_idx = [notcortex_idx; amygdala];
notcortex_idx = sort(notcortex_idx);

cortex_idx = setdiff([1:length(sensorifugal)],notcortex_idx);

% remove noncortical indices
brainspan(:,notcortex_idx) = [];

%% remove unstable genes

% find stable genes as defined on the 34-node parcellation
% keep (available) stable genes

stable_label = label(genes.scale033.stable);
gene_label = table2cell(gene_info(:,4));
notstable_idx = [];
bspan_gidx033 = [];
for k = 1:length(gene_label)
    if ~ismember(stable_label,gene_label(k)) % if gene isn't stable
        notstable_idx = [notstable_idx; k];
    else
        i = find(ismember(stable_label,gene_label(k))); % store indices of stable genes
        bspan_gidx033 = [bspan_gidx033; i];
    end
end

% remove nonstable genes
brainspan(notstable_idx,:) = [];

%% organize data by life stage

fetal_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'fetal'));
infant_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'infant'));
child_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'child'));
adolescent_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'adolescent'));
adult_idx = find(contains(table2cell(sample_info(cortex_idx,9)),'adult'));

fetal = brainspan(:,fetal_idx);
infant = brainspan(:,infant_idx);
child = brainspan(:,child_idx);
adolescent = brainspan(:,adolescent_idx);
adult = brainspan(:,adult_idx);

M = {fetal, infant, child, adolescent, adult};
M_idx = {fetal_idx, infant_idx, child_idx, adolescent_idx, adult_idx};

%% parcellate brain regions

sample_regions = table2cell(sample_info(cortex_idx,8));
regions = unique(string(sample_regions));
region_mapping = zeros(length(sample_regions),1);
for k = 1:length(regions)
    i = find(contains(sample_regions,regions(k)));
    region_mapping(i) = k; % make index-based map of regions
end

regionsIncluded = cell(5,1);
for k = 1:length(M)                                          % for each life stage
    mat = M{k};                                              % get original matrix
    r = region_mapping(M_idx{k});                            % get region mapping
    regionsIncluded{k} = unique(r);                          % this is how many regions there will be
    mat_tmp = zeros(size(mat,1),length(regionsIncluded{k})); % parcellated matrix
    i = 1;
    for j = regionsIncluded{k}'
        aRegion = mat(:,find(r==j));                         % find all columns corresponding to same region
        mat_tmp(:,i) = sum(aRegion,2) ./ sum(aRegion~=0,2);  % fill parcellated matrix with mean exp ignoring 0 
        i = i+1;
    end
    M{k} = mat_tmp';                                         % reassign
end

%% get gene and term scores

% estimate gene scores by multiplying PLS-derived gene weights with gene
% expression matrices

gscore = cell(length(M),1);

for k = 1:length(M)                             % for each life stage
    m = M{k};                                   % get gene expression matrix
    med = repmat(nanmedian(m,2), 1, size(m,2)); % fill missing data with median expression
    m(isnan(m)) = med(isnan(m));
    gscore{k} = m*result.u(bspan_gidx033);      % estimate gene score
end

% convert PLS-derived term scores into 16-node parcellation, averaging
% sibling nodes
tscore = zeros(length(regions),1);
for k = 1:length(regions)
    tscore(k) = mean(result.vsc(mapping==k,1));
end

%% track gene scores across development

% get regions available across all life stages
reg = intersect(regionsIncluded{1},[regionsIncluded{2};regionsIncluded{3};regionsIncluded{4};regionsIncluded{5}]);

for k = 1:5 % keep gene scores in regions available across all life stages
    gscore{k} = gscore{k}(reg);
end
gscore_mat = [gscore{1} gscore{2} gscore{3} gscore{4} gscore{5}];

%% visualize

cm=cbrewer('qual', 'Paired', 16, 'PCHIP');

% gene scores across development

figure;
for k = 1:size(gscore_mat,2)
    hold on
    plot(gscore_mat(:,k),'LineWidth',1.3,'Color',cm(reg(k),:))
end
legend(regions(reg),'Location','northwest');
xticks(1:5)
xticklabels({'fetus','infant','child','adolescent','adult'})
ylabel('estimated gene scores')

% gene-term score correlation across development

figure;
for k = 1:length(gscore)
    subplot(2,3,k)
    t = tscore;
    g = gscore{k};
    notIncluded = setdiff(1:length(regions),regionsIncluded{k});
    if size(notIncluded,2) > 0
        t(notIncluded) = [];
    end
    [rho,pval] = corr(g,t);
    scatter(g,t,50,regionsIncluded{k},'filled')
    xlabel('gene score')
    ylabel('term score')
    title([' rho=', num2str(rho), ' p=', num2str(pval)])
end
colormap(cm);
