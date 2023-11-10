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
    'usroads-48.mat'};

for i = 1:length(matrix)
    % Load Matrix
    load(sprintf("../matrix/%s",matrix{i}));
    % Build Markov chain
    if exist('Problem','var')
        A = spones(Problem.A);
        matrixname = Problem.name;
    else
        matrxiname = matrix{i}(1:end-4);
    end
    n = size(A,1);
    p = dissect(A); % Reorder sparse matrix
    PA = A(p,p);
    PA = spdiags(PA*ones(n,1),0,n,n)\PA;
    [pi,~] = eigs(PA',1,'largestabs','MaxIterations',10000);
    pi = pi/sum(pi);

    fprintf("\nWorking on matrix %s\n",matrixname);
    %% Direct solution
    tic;
    try
        kval = kemenydirect(PA);
    catch
        kval = -1;
    end
    time0 = toc;
    fprintf("Direct Kemeny computation: k = %f time = %1.2f\n\n",kval,time0);
    %% SPARSE Solution
    tic;
    try
        k = recursivekemeny(PA,pi);
    catch
        k = -1;
    end
    time3 = toc;
    fprintf("Recursive Kemeny computation: k = %f time = %1.2f\n\n",k,time3);
    fprintf("Abs Error is %e Rel Error is %e\n",abs(kval-k),abs(kval-k)/kval);
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
    fprintf("(HODRL) Compression time %1.2f\n",time2);
    fprintf("(HODLR) Recursive Kemeny computation: k = %f time = %1.2f\n\n",k,time1);
    fprintf("Abs Error is %e Rel Error is %e\n",abs(kval-k),abs(kval-k)/kval);

    clear Problem PA
end