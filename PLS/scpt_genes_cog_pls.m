% This script runs PLS on gene expression and neurosynth probability maps.
% Significance of latent variables is assessed using a
% spatial autocorrelation-preserving permutation test. Correlation between
% gene and term scores are cross-validated using a distance-based set
% assignment (see fcn_crossval_pls_brain_obvs.m). Most contributing terms
% are retained and scores are distributed among structural and functional
% networks.

%% load

load('gene_expression.mat') % node by gene expression matrix
load('neurosynth.mat')      % node by term probability matrix
load('nodes.mat')           % relevant node indices
load('genes.mat')           % relevant gene indices
load('terms.mat')           % relevant term indices and names
load('spins.mat')           % spin test indices
load('coords.mat')          % (x,y,z) coordinates for brain regions

%% PLS analysis

% shown here for neurosynth - replicated for brainmap

% set up PLS analysis
n = nodes.scale125.lefthem;
g = genes.scale125.stable;
t = terms.all; 

X = zscore(expression125(:,g));
Y = zscore(cogact125(n,t));

nnodes = length(n);
nterms = length(t);
ngenes = length(g);

% behav pls
option.method = 3;
option.num_boot = 10000;
option.num_perm = 0;               % zero permutations because they will be run manually later to account for spatial autocorrelation
option.stacked_behavdata = Y;

exp{1} = X;

result = pls_analysis(exp, nnodes, 1, option); % this is the PLS result that is used in all other analyses
% save('result.mat','result')

%% spin test
% this code comes from pls_analysis.m and is modified to account for a
% spatial autocorrelation-preserving permutation test

spins = scale125spins;          % spatial autocorrelation-preserving permutation assignments
nspins = 10000;                 % number of permutations ("spins")
s_spins = zeros(nterms,nspins); % singular values
option.method = 3;              % set up PLS
option.num_boot = 0;
option.num_perm = 0;
exp{1} = X;
for k = 1:nspins    
    option.stacked_behavdata = Y(spins(:,k),:);  % permute neurosynth matrix
    
    datamatsvd=rri_xcor(option.stacked_behavdata,exp{1},0); % refer to pls_analysis.m
    [r,c] = size(datamatsvd);
    if r <= c
        [pu, sperm, pv] = svd(datamatsvd',0);
    else
        [pv, sperm, pu] = svd(datamatsvd,0);
    end
    
    %  rotate pv to align with the original v
    rotatemat = rri_bootprocrust(result.v,pv);
 
    %  rescale the vectors
    pv = pv * sperm * rotatemat;

    sperm = sqrt(sum(pv.^2));
    
    s_spins(:,k) = sperm;
end

sprob = zeros(nterms,1); % p-value for each latent variable

for k = 1:nterms % get permuted (via spin test) p-values
    sprob(k) = (1+(nnz(find(s_spins(k,:)>=result.s(k)))))/(1+nspins);
end  

%% cross validation

% cross-validate the correlation between gene and term scores using
% distance-based set assignment
[rtrain,rtest] = fcn_crossval_pls_brain_obvs(exp,Y,100,0.75,1,coords125(116:226,:));

% visualize boxplots as distributions of the correlations computed using
% the training and testing set
figure;
boxplot([rtrain rtest])
set(gca,'xticklabel',{'train','test'})
ylabel('correlation')
title('cross validation')

%% get terms

[B,I] = sort(result.boot_result.orig_corr(:,lv),'descend');  % sort loadings
npos = length(find(result.boot_result.orig_corr(:,lv) > 0)); % number of positive loadings
nneg = length(find(result.boot_result.orig_corr(:,lv) < 0)); % number of negative loadings
tpos = I(1:floor(0.25*npos));                                % get top 25% of positive loadings
tneg = I(end-floor(0.25*nneg):end);                          % get top 25% of negative loadings
pos_terms = t(tpos);                                         % these are the positive terms contributing most
neg_terms = t(tneg);                                         % these are the negative terms contributing most

%% distribute PLS-derived gene and term scores

% intrinsic (resting-state) networks from Yeo et al., 2011
load('rsn.mat')       % load rsn network assignment
load('rsn_names.mat') % load rsn names
rsn = rsn(109:end);   % keep left hemisphere only

figure;               % distribute scores (as boxplots) in 7 networks
subplot(2,1,1)
boxplot(result.usc(:,1),rsn) % gene scores
xticklabel(rsn_names)
title('gene scores')
subplot(2,1,2)
boxplot(result.vsc(:,1),rsn) % term scores
xticklabel(rsn_names)
title('term scores')

% von economo cytoarchitectonic classes

load('voneconomo.mat')            % load von economo network assignment
load('ve_names.mat')              % load von economo names
voneconomo = voneconomo(109:end); % keep left hemisphere only

figure;                           % distribute scores (as boxplots) in 7 networks
subplot(2,1,1)
boxplot(result.usc(:,1),voneconomo)  % gene scores
xticklabel(ve_names)
title('gene scores')
subplot(2,1,2)
boxplot(-result.vsc(:,1),voneconomo) % term scores
xticklabel(ve_names)
title('term scores')

% mesulam levels of laminar differentiation
% level names from https://github.com/MICA-MNI/micaopen/tree/master/MPC

mesulam = table2array(readtable('mesulam_mapping.csv'));              % get mesulam class assignments
mesulam = mesulam(109:end);                                           % keep left hemisphere only
mesulam_names = {'paralimbic','heteromodal','unimodal','idiotypic'};  % class names

figure;                           % distribute scores (as boxplots) in 4 classes
subplot(2,1,1)
boxplot(result.usc(:,1),mesulam)  % gene scores
xticklabel(mesulam_names)
title('gene scores')
subplot(2,1,2)
boxplot(result.vsc(:,1),mesulam)  % term scores
xticklabel(mesulam_names)
title('term scores')
