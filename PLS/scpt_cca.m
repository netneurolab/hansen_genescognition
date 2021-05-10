% This script runs CCA on the gene expression and functional activation
% matrices to compare whether results are similar between CCA and PLS
% methods.

%% load

load('gene_expression.mat') % node by gene expression matrix
load('neurosynth.mat')      % node by term probability matrix
load('nodes.mat')           % relevant node indices
load('genes.mat')           % relevant gene indices
load('terms.mat')           % relevant term indices and names
load('result.mat')          % PLS result

X = zscore(expression125(:,g)); % gene expression matrix
Y = zscore(cogact125(n,t));     % neurosynth matrix
 
[g_coeff,g_score,g_latent] = pca(X);
[t_coeff,t_score,t_latent] = pca(Y);

ncmp = 50; % number of components
[A,B,r,U,V] = canoncorr(zscore(g_score(:,1:ncmp)),zscore(t_score(:,1:ncmp)));

% go backwards to get weights (first component only)
% method from Smith et al 2015
gweights = corr(U(:,1),X);
tweights = corr(V(:,1),Y);

% project onto original data to get scores
gscores = X*gweights';
tscores = Y*tweights';

% check if CCA and PLS give similar results
corr(gscores,result.usc(:,1))
corr(tscores,result.vsc(:,1))

corr(gweights', result.u(:,1))
corr(tweights', result.v(:,1))