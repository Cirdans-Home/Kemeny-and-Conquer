%% BATCH TEST FOR TOEPLITZ CLUSTER

% Load src files
addpath('../hm-toolbox/')
addpath('../src/')

 matrix = {'Pisa.mat',...
    'big.mat',...
    'cage10.mat',...
    'cage11.mat',...
    'gre_1107.mat',...
    'nopoly.mat',...
    'pesa.mat',...
    'usroads-48.mat',...
    'USpowerGrid.mat',...
    'minnesota.mat'};

fprintf("Batch test of the different algorithms\n")
fprintf("Test run on the date %s\n",string(datetime));

for i = 1:length(matrix)
    % Load Matrix
    load(sprintf("../matrix/%s",matrix{i}));
    % Build Markov chain
    if exist('Problem','var')
        A = spones(Problem.A);
        matrixname = Problem.name;
	G = graph(A);
    	[bin,binsize] = conncomp(G);
    	idx = binsize(bin) == max(binsize);
    	SG = subgraph(G, idx);
    	A = SG.adjacency();
    else
        matrixname = matrix{i}(1:end-4);
    end
    n = size(A,1);
    p = dissect(A); % Reorder sparse matrix
    PA = A(p,p);
    PA = spdiags(PA*ones(n,1),0,n,n)\PA;
    [pi,~] = eigs(PA',1,'largestabs','MaxIterations',10000);
    pi = pi/sum(pi);

    fprintf("\nWorking on matrix %s of size n = %d\n",matrixname,n);
    %% Direct solution
    tic;
    try
        kval = kemenydirect(PA);
    catch
        kval = -1;
    end
    time0 = toc;
    fprintf("\tDirect Kemeny computation: k = %f time = %1.2f\n",kval,time0);
    %% SPARSE Solution
    tic;
    try
        k = recursivekemeny(PA,pi);
    catch
        k = -1;
    end
    time3 = toc;
    fprintf("\t(Sparse) Recursive Kemeny computation: k = %f time = %1.2f\n",k,time3);
    fprintf("\t(Sparse) Abs Error is %e Rel Error is %e\n",abs(kval-k),abs(kval-k)/kval);
    %% SPARSE-DIRECT Solution
    tic;
    try
        k = recursivekemenydirect(PA,pi);
    catch
        k = -1;
    end
    time4 = toc;
    fprintf("\t(Sparse-Direct) Recursive Kemeny computation: k = %f time = %1.2f\n",k,time4);
    fprintf("\t(Sparse-Direct) Abs Error is %e Rel Error is %e\n",abs(kval-k),abs(kval-k)/kval);
    %% HODLR Solution
    tic;
    PA = hodlr(PA);
    time2 = toc;
    tic;
    try
        k = recursivekemenyhodlr(PA,pi);
    catch
        k = -1;
    end
    time1 = toc;
    fprintf("\t(HODRL) Compression time %1.2f\n",time2);
    fprintf("\t(HODLR) Recursive Kemeny computation: k = %f time = %1.2f\n",k,time1);
    fprintf("\tAbs Error is %e Rel Error is %e\n",abs(kval-k),abs(kval-k)/kval);

    clear Problem PA
end

fprintf("Test completed on the date %s\n",string(datetime));
