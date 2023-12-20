%% Low-Precision randomized test
% Test with Hutch++

clear; clc; close all hidden;

addpath("../HutchPlusPlus/core/")
addpath("../src/")

matrix_list = {'USpowerGrid.mat','nopoly.mat','minnesota.mat'};


h1 = waitbar(0,"Test Progress");
h2 = waitbar(0,"Repetitions");
for i=1:length(matrix_list)
    load(sprintf("../matrix/%s",matrix_list{i}));
    
    % Generate test-problem
    A = spones(Problem.A);
    G = graph(A);
    [bin,binsize] = conncomp(G);
    idx = binsize(bin) == max(binsize);
    SG = subgraph(G, idx);
    A = SG.adjacency();
    n = size(A,1);
    p = dissect(A);
    PA = A(p,p);
    DA = spdiags(1./sqrt(PA*ones(n,1)),0,n,n);
    PA = DA*PA*DA;

    % Full-Kemeny
    tic;
    ktrue = kemenydirect(PA);
    time0 = toc;
    % Randomized kemeny
    delta = 1/4;
    epsilon = 1e-1;
    l = @(delta,epsilon) round(sqrt(log(1/delta))/epsilon + log(1/delta));
    attempts = 100;
    khutch = 0;
    tic;
    waitbar(0,h2);
    for att = 1:attempts
        khutch = khutch + kemenyfullestimate(PA,l(delta,epsilon),"ITERATIVE");
        waitbar(att/attempts,h2);
    end
    khutch = khutch/attempts;
    time1 = toc/attempts;
    
    fprintf("%s & %d & %1.2f & %1.2f & %1.2e\\\\\n", ...
        Problem.name,n,time0,time1,abs(ktrue-khutch)/abs(ktrue));
    waitbar(i/length(matrix_list),h1)
end

rmpath("../HutchPlusPlus/core/")
rmpath("../src/")
