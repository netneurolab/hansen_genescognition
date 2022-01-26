% This function calculates the ratio of genes preferentially expressed
% in different cell types, according to Seidlitz 2020 Nat Comms.
% I wrote this function after the paper was published (which is to say,
% I never actually used it in my project). It's up here because it might
% be helpful to others.

% input:
%   GeneListSub: gene subset to do ctd (string with gene names)
%   GeneListFull: all possible genes (string with gene names)
%   nperm: number of permutations for null model

% output:
%   ctd: cell type deconvolution
%   ctd_null: permuted null model

function [ctd, ctd_null, ctd_pvals, cellnames] = fcn_ctd(GeneListSub, GeneListFull, nperm)

% get cell types

celltypes = table2cell(readtable('celltypes_PSP.csv')); % load specific cell type expression
genenames = cellstr(celltypes(:,1));

for k = 1:length(genenames)                                                   % for each gene with specific cell type expression
    if ismember(genenames(k),GeneListFull)                                    % if gene overlaps with GeneListFull
        celltypes(k,3) = num2cell(find(strcmp(GeneListFull,genenames(k))));   % add index of gene
    else
        celltypes(k,3) = {0};                                                 % otherwise, add 0
    end
end
bad_idx = cell2mat(celltypes(:,3))==0; % remove genes not in GeneListFull
celltypes(bad_idx,:) = [];
genenames_infull = cellstr(celltypes(:,1));

[cellnames,~,i] = unique(cellstr(celltypes(:,2))); % index genes by which cell type they're expressed in

ntypes = length(cellnames);

% get cell type deconvolution
ctd = zeros(ntypes,1);
for k = 1:ntypes
    ctd(k) = length(intersect(GeneListSub,genenames_infull(i==k)))/length(GeneListSub);
end

% get null model
ctd_null = zeros(ntypes,nperm);
for k = 1:nperm
    y = datasample(1:length(GeneListFull),length(GeneListSub),'Replace',false);                         % get random gene set the size of the positive gene set
    for j = 1:ntypes                                                                                      % for each cell type
        ctd_null(j,k) = length(intersect(GeneListFull(y),genenames_infull(i==j)))/length(GeneListSub);  % find ratio of genes expressed in cell type to all genes
    end
end

% get pvals
ctd_pvals = zeros(ntypes,1);
for k = 1:ntypes
    ctd_pvals(k) = (1 + sum(abs(ctd_null(k,:) - mean(ctd_null(k,:))) >= abs(ctd(k) - mean(ctd_null(k,:))))) / (nperm + 1);
end

% plot
figure;
o = [1,3,6,7,2,4,5]; % order in which cell types appear
scatter(1:ntypes,ctd(o),30,'filled')
hold on
boxplot(squeeze(ctd_null(o,:))')
set(gca,'xticklabel',cellnames(o))
xtickangle(90)
title('cell type deconvolution')
ylabel('ratio')

end
