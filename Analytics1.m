%% Import the data from the .txt file and convert it into a matrix
cities = readtable('10cityGCs.txt');
cities = table2array(cities);

%% Labelling of different cities
% The different cities are given the following labels:
% 1: Edinburgh  2: Glasgow      3: Cardiff  4: Bristol      5: Nottingham
% 6: Birmingham 7:Sheffield     8:Leeds     9:Manchester    10:London

cityNamesPt1 = ["Edinburgh","Glasgow","Cardiff","Bristol","Nottingham"];
cityNamesPt2 = ["Birmingham","Sheffield","Leeds","Manchester","London"];
cityNames = [cityNamesPt1,cityNamesPt2];

%% Change the number below to create the necessary plots for a particular
% city
cityNo = 1;                     % Edinburgh in this case
cityName = cityNames(cityNo);

%% Extracting edge set from the matrix for each of the cities:
% Edinburgh, Glasgow and Cardiff

city = cities(cities(:,1)==cityNo,2:3);
N = max(max(city));

%% Converting the edge sets into adjacency matrices for each city

A = sparse(city(:,1),city(:,2),1,N,N);

%% Creating initial vectors of all ones for each adjacency matrix,
% to be used as initial estimates of the eigenvectors of each matrix
% corresponding to their maximum eigenvalues

Ones = ones(N,1);
V = ones(N,1);

%% Running the power method

for i = 1:1e2
    V = A*V;
    V = V/norm(V,2);
end

%% Returning an estimate of the Perron-Frobenius eigenvalue (i.e. spectral
% radius) via the Rayleigh quotient

rho = (V'*A*V)/norm(V,2);
fprintf("Spectral radius = %f for %s.\n",rho,cityName)

%% Calculating the D-matrices for each city, where D = diag(A*w) with A 
% being the adjacency matrix for each city and w being a vector of ones

D = diag(A*ones(N,1));
I = eye(N);

D = sparse(D);
I = sparse(I);

%% Calculating Katz centralities for each matrix, and also the 
% Non-Backtracking Walk (NBTW) Centrality for a value of alpha (for each
% matrix) slightly less than the maximum eigenvalue
eps = 1e-2;
a = rho*(1 - eps);
Katz = (I - a*A)\Ones;
NBTW = (I - a*A + a^2 *(D - I))\Ones;

indices = 1:N;
[KatzSort, katzIndices] = sort(Katz);
[nbtwSort, nbtwIndices] = sort(NBTW);

%% Plotting one type of centrality against the other, for each city

figure(1);
scatter(katzIndices,nbtwIndices,'.r')

xlabel('Katz centrality (rank)','FontSize',12);
ylabel('NBTW centrality (rank)','FontSize',12);
title(cityName+" (r = "+num2str(corr(katzIndices,nbtwIndices))+")")
axis square
savefig(cityName+'1.fig')
saveas(gcf,cityName+'1.png')

%% Changing alpha to demonstrate that Q, the matrix corresponding to each 
% type of centrality blows up

figure(2)
i = 0:0.1:3;
eps = 10.^(-i);
a = (1-eps).*(1/rho); 
qkConds = zeros(1,length(i));
for j = 1:length(i)
    qkConds(j) = condest(inv(I - a(j)*A));
    fprintf("i = %f / %f\n",i(j),max(i))
end
plot(1-eps,qkConds,'b')

xlim([0.8,1])
xlabel("$\alpha/(1/\rho)$","Interpreter","LaTeX")
ylabel("$\textrm{cond}\left(\left(\mathbf{I}-\alpha\mathbf{A}\right)^{-1}\right)$","Interpreter","LaTeX")
title("Katz Centrality ("+cityName+")")
hold off
savefig(cityName+'2.fig')
saveas(gcf,cityName+'2.png')
%% Changing alpha to demonstrate that Q, the matrix corresponding to each 
% type of centrality blows up

figure(3)
eps = 10.^(-i);
a = (1-eps).*(1/rho); 
qnConds = zeros(1,length(i));
for j = 1:length(i)
    qnConds(j) = condest(inv(I - a(j)*A + a(j)^2 * (D-I)));
    fprintf("i = %f / %f\n",i(j),max(i))
end
plot(1-eps,qnConds,'r')

xlim([0,1])
xlabel("$\alpha/(1/\rho)$","Interpreter","LaTeX")
ylabel("$\textrm{cond}\left((\mathbf{I}-\alpha\mathbf{A}+\alpha^{2}(\mathbf{D}-\mathbf{I}))^{-1}\right)$","Interpreter","LaTeX")
title("Non-backtracking Walk Centrality ("+cityName+")")
hold off
savefig(cityName+'3.fig')
saveas(gcf,cityName+'3.png')

%% Calculating the Fiedler eigenvalue and eigenvector for discrete Laplacian (D - A)

figure(4)
[eigVecs,eigVals] = eig(full(D - A));
eigVals = diag(eigVals);
eigValFE = eigVals(2);
eigVecFE = eigVecs(:,2);

fprintf("Fiedler eigenvalue = %f for %s.\n",eigValFE,cityName)

Graph = graph(A);
GraphPlot = plot(Graph);

Colours = jet(N);

for i = 1:N
    Max = max(eigVecFE);
    Min = min(eigVecFE);
    Range = Max - Min;
        
    % if construct to avoid possibility of a 1 or (N+1) index - this is 
    % fine, as the number of nodes is very large
    if eigVecFE(i) == Max
       whichColour = N;            
    elseif eigVecFE(i) == Min
       whichColour = 1;
    else
       whichColour = ceil(N*(eigVecFE(i) - Min)/Range);
    end       

    highlight(GraphPlot,i,'NodeColor',Colours(whichColour,:));
end

bar = colorbar;
bar.Label.String = 'Element of Fiedler eigenvector';
colormap(jet)
caxis([Min,Max])
title(cityName)
hold off
savefig(cityName+'4.fig')
saveas(gcf,cityName+'4.png')