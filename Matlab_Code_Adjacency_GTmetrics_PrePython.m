%% MATLAB Part 1: Adjacency Matrix extraction from 3D data
% Required Input: 
% 1) Binary tif image
% 2) Label tif image
% 3) Analyze Region csv file

% Read 3D-Tiff + MorphoLibJ Analysis Regions3D (or Avizo, then change analysis region file) data file
clear;
clc;
fileID = fopen('TestData_KirchhoffInput.txt','w');

% Load binary 3D tif, label 3D tif and corresponding Analyze Regions csv
Verbunden = tiffreadVolume('TestData_binary.tif');
stack = tiffreadVolume('TestData_labels.tif');
AnalyzeRegionsData = readmatrix('Silber_labels-morpho.csv');
workspaceID = "TestData_workspace";

% Read center of gravity (cogs, column 28,29,30) + equivalent radius (column 31) from analyze region file
coordinatesXYZ = AnalyzeRegionsData(:,28:30);
diameter = 2*AnalyzeRegionsData(:,31);

% INITIALIZE + COMPUTE
% Initialize Particlenumbers based on Image stack
[maxParticle, maxInd] = max(stack, [], "all", "linear");
% Initialize Check Matrix
contactArea = zeros(maxParticle,maxParticle);
% Initialize Direction Moving 
dx = [2,-2,0,0,0,0];
dy = [0,0,2,-2,0,0];
dz = [0,0,0,0,2,-2];
% Initialize Middle Point
dxM = [1,-1,0,0,0,0];
dyM = [0,0,1,-1,0,0];
dzM = [0,0,0,0,1,-1];
% Initialize NumContacts
NumContacts = 0;
% Setup size of Volume
sizeX = length(stack(:,1,1));
sizeY = length(stack(1,:,1));
sizeZ = length(stack(1,1,:));

% Get value of Verbunden Image Stack (Because sometimes instead of 0 and 1 it is 0 and 255)
tempStack = Verbunden(1:50,1:50,1:50);
tempValue = max(tempStack,[],"all","linear");

for x = 3:(sizeX-2)
    for y = 3:(sizeY-2)
        for z = 3:(sizeZ-2)
            particleIndex1 = stack(x,y,z);
            if particleIndex1 ~= 0  % if not 0
                for k = 1:6
                    x2 = x + dx(k);
                    y2 = y + dy(k);
                    z2 = z + dz(k);
                    xM = x + dxM(k);
                    yM = y + dyM(k);
                    zM = z + dzM(k);
                    particleIndex2 = stack(x2,y2,z2);
                    particleIndexM = Verbunden(xM,yM,zM);
            
                    if particleIndex2 ~= 0 && particleIndex2~= particleIndex1 &&  particleIndexM == tempValue
                        contactArea(particleIndex1,particleIndex2) = contactArea(particleIndex1,particleIndex2)+1;
                    end
                    


                end
            end
            
        end
    end
end

% Cleaning of false-positives
particleIndex1 = 0;
particleIndex2 = 0;
for particleIndex1 = 1:maxParticle
    for particleIndex2 = (particleIndex1 + 1):maxParticle
        tempCA1 = contactArea(particleIndex1,particleIndex2);
        tempCA2 = contactArea(particleIndex2,particleIndex1);
        tempProduct = tempCA1*tempCA2;
        if tempProduct == 0
            contactArea(particleIndex1,particleIndex2) = 0;
            contactArea(particleIndex2,particleIndex1) = 0;
        end

    end

end

% Export to .txt for Kirchhoff-Code

for m = 1:length(coordinatesXYZ(:,1))
    fprintf(fileID,'%s\n','node');
    fprintf(fileID,'%14.12f\n',coordinatesXYZ(m,1));
    fprintf(fileID,'%14.12f\n',coordinatesXYZ(m,2));
    fprintf(fileID,'%14.12f\n',coordinatesXYZ(m,3));
    fprintf(fileID,'%14.12f\n',diameter(m));
end
Export partile nodes + contact voxel
particleIndex1 = 0;
particleIndex2 = 0;
for particleIndex1 = 1:maxParticle
    
    for particleIndex2 = (particleIndex1+1):maxParticle % NICHT particleIndex1 = 1:maxParticle --> Doppelter Export von
        tempCA1 = contactArea(particleIndex1,particleIndex2);
        if tempCA1 > 0
            % fprintf(fileID,'%s\n','edge');
            % fprintf(fileID,'%d\n',particleIndex1);
            % fprintf(fileID,'%d\n',particleIndex2);
            % fprintf(fileID,'%d\n',tempCA1);
            NumContacts = NumContacts + 1; % total Number of Contacts
            %RealContactArea(particleIndex1) = tempCA1; % in Voxel-Units
            
        end
        
    end

end

fclose(fileID); % close Kirchhoff Input file
save(workspaceID); % Save workspace

%% MATLAB Part 2: Build graph G(N,E) and calculate graph metrics

% Load from workspace if not open
 load(workspaceID,"contactArea");
G = graph(logical(contactArea));
GNumNodes = numnodes(G)
GNumEdges = numedges(G)
GDegree = G.degree; % Degree of the graph == coordination number k
CCO = centrality(G,'closeness'); % closeness centrality
CCB = centrality(G,'betweenness'); % betweenness centralityclead
CCE = centrality(G,'eigenvector'); % Eigenvector centrality
CCP = centrality(G,'pagerank'); % Pagerank centrality
Gdens = nnz(adjacency(G))./numel(adjacency(G)) % graph density
GDists = distances(G); % All distances between nodes
B = squareform(GDists); % takes upper triangle of a matrix without the diagonal
DistanceArray = B(~isinf(B));
L = mean(DistanceArray);
GNetDiam = max(DistanceArray);
%L = mean(GDists(~isinf(GDists))) wrong way to calculate L
RecDist = 1./DistanceArray; % Reciprocal Distances
%RecDist = RecDist(~isinf(RecDist)); % Remove infs
GEff = (sum(RecDist))/(GNumNodes*(GNumNodes-1))
%GNetDiam = max(GDists(~isinf(GDists))) % network diameter
e = size(G.Edges,1); % get number of vertices
k = mean(G.degree);
% Save workspace
filename_without_extension = extractBefore(workspaceID, ".");
workspaceID2 = filename_without_extension + "_GTA.mat";
save(workspaceID2);

%% MATLAB Part 3: Build workspace of only connected graph, i.e. Prepare subgraphs with largest connected component for Python Analysis

[bin,binsize] = conncomp(G,'Type','weak');
idx = binsize(bin) == max(binsize);
SG = subgraph(G, idx);
adjac = SG.adjacency;
save("TestData_GTA_CentralityAnalysis_Subgraph.mat")

