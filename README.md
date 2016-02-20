We also give your the PPI network from HPRD database in the file HPRDNetworkWithComplexes.txt. The file HPRDID.txt contain the mapping from gene symbols/names to RefSeq IDs that are used in HPRDNetworkWithComplexes.txt.

How to compile:

g++ createGraph.cpp -o createGraph
g++ extractModules_OptDis.cpp -fopenmp -O3 -g -o extractModules_OptDis
g++ extractModules_Dense.cpp -o extractModules_Dense
g++ createModules.cpp -o createModules

First, run createGraph to make a base folder that contain expression profile and network information that OptDis and Dense programs can extract the modules:

./createGraph HPRDNetworkWithComplexes.txt TEST.txt TEST.txt TEST HPRDID.txt 61 21

The 1st parameter is the network file, 2nd and 3rd parameter are the expression profiles (before I intended to use copy number here, but i leave one for copy number profiles). The 4th parameter is the base folder. The 5th parameter here is the network id file. The 6th parameter is the number of samples in the positive class. The 7th parameter is the number of samples in the negative class.

To extract the connected subnetwork marker as in the ISMB's paper, you run the following command:

./extractModules_OptDis TEST 61 4

The 1st parameter here is the base folder that was created by the above createGraph program. The 2nd parameter here is the number samples in the positive class. The 3rd parameter here is the maximum size of a subnetwork.

To extract the densely connected subnetwork marker as in the ECCB's paper, you run the following command:

./extractModules_Dense TEST 61 0.7

Again the 1st parameter here is the base folder that was created by the above createGraph program. The 2nd parameter here is the number samples in the positive class. The 3rd parameter here is the minimum density constraint. 

To create a subnetwork activity file:

for densely connected subnetworks:
./createModules TEST 61 TEST_Dense.csv Dense

So the 1st parameter here is the base folder that was created by the above createGraph program. The 2nd parameter here is the number samples in the positive class. The 3rd parameter here is output subnetwork activity file. The 4th parameter is the name of the method; Dense here for densely connected subnetworks and OptDis for connected subnetworks. The activity of a subnetwork in a sample is the average of the expressions of genes from that subnetwork.

for connected subnetworks:
./createModules TEST 61 TEST_OptDis.csv OptDis

