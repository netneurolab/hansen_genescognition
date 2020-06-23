load('gene_expression.mat')
load('neurosynth.mat')
load('nodes.mat')
load('genes.mat')
load('terms.mat')
load('spins.mat')
load('coords.mat')

%% PLS analysis

X = zscore(expression125);
Y = zscore(cogact125); % or brainmap(1:111,:)

n = nodes.scale125.lefthem;
g = genes.scale125.stable;
t = terms.all; % or terms.brainmap.all

nnodes = length(n);
nterms = length(t);
ngenes = length(g);

% behav pls
option.method = 3;
option.num_boot = 10000;
option.num_perm = 0;
option.stacked_behavdata = Y(n,t); % or Y

exp{1} = X(:,g);

result = pls_analysis(exp, nnodes, 1, option);

%% spin test

spins = scale125spins;
nspins = 10000;
s_spins = zeros(nterms,nspins);
option.method = 3;
option.num_boot = 0;
option.num_perm = 0;
exp{1} = X(:,g);
for k = 1:nspins    
    option.stacked_behavdata = zscore(Y(n(spins(:,k)),:));
    
    datamatsvd=rri_xcor(option.stacked_behavdata,exp{1},0);
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

sprob = zeros(nterms,1);

for k = 1:nterms
    sprob(k) = nnz(find(s_spins(k,:)>=result.s(k)))/nspins;
end  

%% cross validation

[rtrain,rtest] = fcn_crossval_pls_brain_obvs(exp,zscore(Y(n,t)),100,0.75,1,coords125(116:226,:));

figure;
boxplot([rtrain rtest])
set(gca,'xticklabel',{'train','test'})
ylabel('correlation')
title('cross validation')

%% get terms

[B,I] = sort(result.boot_result.orig_corr(:,lv),'descend');
npos = length(find(result.boot_result.orig_corr(:,lv) > 0));
nneg = length(find(result.boot_result.orig_corr(:,lv) < 0));
tpos = I(1:floor(0.25*npos));
tneg = I(end-floor(0.25*nneg):end);
pos_terms = t(tpos);
neg_terms = t(tneg);

%% distribute scores

% intrinsic (resting-state) networks
load('rsn.mat')
load('rsn_names.mat')
rsn = rsn(109:end);

figure;
subplot(2,1,1)
boxplot(-result.usc(:,1),rsn)
xticklabel(rsn_names)
title('gene scores')
subplot(2,1,2)
boxplot(-result.vsc(:,1),rsn)
xticklabel(rsn_names)
title('term scores')

% von economo cytoarchitectonic classes

load('voneconomo.mat')
load('ve_names.mat')
voneconomo = voneconomo(109:end); % left hem only

figure;
subplot(2,1,1)
boxplot(-result.usc(:,1),voneconomo)
xticklabel(ve_names)
title('gene scores')
subplot(2,1,2)
boxplot(-result.vsc(:,1),voneconomo)
xticklabel(ve_names)
title('term scores')

% mesulam levels of laminar differentiation
% level names from https://github.com/MICA-MNI/micaopen/tree/master/MPC

mesulam = table2array(readtable('mesulam_mapping.csv'));
mesulam = mesulam(109:end); % left hem only
mesulam_names = {'paralimbic','heteromodal','unimodal','idiotypic'};

figure;
subplot(2,1,1)
boxplot(result.usc(:,1),mesulam)
xticklabel(mesulam_names)
title('gene scores')
subplot(2,1,2)
boxplot(result.vsc(:,1),mesulam)
xticklabel(mesulam_names)
title('term scores')
