% this function cross-validates the correlation between gene and term
% scores by assigning observations (brain areas) to training/testing sets
% using a distance-based method. 

% Inputs:
% x - cell array containing node x genes expression matrix
% y - node x term functional activation matrix
% nsplit - number of train/test splits
% trainpct - 0<trainpct<1, percent of nodes in training set
% lv - latent variable from PLS to use
% coords - nodes x 3 (x,y,z) coordinates of brain regions

% Outputs:
% rtrain - nsplit x 1 vector of correlations from the training set
% rtest - nsplit x 1 vector of correlations from the testing set


function [rtrain,rtest] = fcn_crossval_pls_brain_obvs(x,y,nsplit,trainpct,lv,coords)

    % set up behav PLS
    option.method = 3;
    option.num_perm = 0;
    option.num_boot = 0;

    rtrain = zeros(nsplit,1);
    rtest = zeros(nsplit,1);
    
    nnodes = size(x{1},1);
    
    P = dist(coords'); % nnodes x nnodes matrix of Euclidean distances between every pair of nodes

    for k = 1:nsplit
        
        % pick a random node
        myNode = randperm(length(coords));
        myNode = myNode(1);
        
        % sort nodes by distance to myNode
        distances = P(myNode,:); % distance of each node to my randomly chosen source node
        [~,idx] = sort(distances);
        
        % make training/testing X and Y matrices based on the 75% of nodes
        % closest to the source node/25% of nodes furtherst from the source
        % node
        trainx{1} = x{1}(idx(1:floor(trainpct*nnodes)),:); % training expression matrix
        xtest = x{1}(idx(floor(trainpct*nnodes)+1:end),:); % testing expression matrix

        option.stacked_behavdata = y(idx(1:floor(trainpct*nnodes)),:); % training functional activation matrix
        ytest = y(idx(floor(trainpct*nnodes)+1:end),:); % testing functional activation matrix

        trainresult = pls_analysis(trainx,floor(trainpct*nnodes),1,option);

        % correlate scores on training set
        rtrain(k) = corr(trainresult.usc(:,lv),trainresult.vsc(:,lv));

        % project weights, correlated predicted scores in the test set
        rtest(k) = corr(xtest*trainresult.u(:,lv),ytest*trainresult.v(:,lv));

    end
end
    
    
    
    