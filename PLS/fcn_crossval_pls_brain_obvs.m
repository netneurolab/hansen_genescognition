function [rtrain,rtest] = fcn_crossval_pls_brain_obvs(x,y,nsplit,trainpct,lv,coords)
    % scpt_crossval_r.m (from Misic) as reference
    % cross-validate corr(usc,vsc) where observations are brain areas
    % 0<trainpct<1
    % coords - brain areas x 3
    
    % due to spatial autocorrelation of brain areas i.e. non independent
    % observations, selection of train and test sets is non-random

    option.method = 3;
    option.num_perm = 0;
    option.num_boot = 0;

    rtrain = zeros(nsplit,1);
    rtest = zeros(nsplit,1);
    
    nnodes = size(x{1},1);
    
    P = dist(coords');

    for k = 1:nsplit
        
        % pick a random node
        myNode = randperm(length(coords));
        myNode = myNode(1);
        
        % sort nodes by distance to myNode
        distances = P(myNode,:);
        [~,idx] = sort(distances);
        
        % pick the (1-trainpct) furthest to it

        trainx{1} = x{1}(idx(1:floor(trainpct*nnodes)),:);
        xtest = x{1}(idx(floor(trainpct*nnodes)+1:end),:);

        option.stacked_behavdata = y(idx(1:floor(trainpct*nnodes)),:);
        ytest = y(idx(floor(trainpct*nnodes)+1:end),:);

        trainresult = pls_analysis(trainx,floor(trainpct*nnodes),1,option);

        % correlate scores on training set
        rtrain(k) = corr(trainresult.usc(:,lv),trainresult.vsc(:,lv));

        % project weights, correlated predicted scores in the test set
        rtest(k) = corr(xtest*trainresult.u(:,lv),ytest*trainresult.v(:,lv));

    end
end
    
    
    
    