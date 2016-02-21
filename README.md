## OPTIMALLY DISCRIMINATIVE SUBNETWORK MARKERS

The expression profile of the sample data is in the file TEST.txt. It has 61 positive samples and 21 negative samples. Note that you need to put all the samples from a class together. Here I put all the samples in positive class first then come the samples from negative class.
We also give your the PPI network from HPRD database in the file HPRDNetworkWithComplexes.txt. The file HPRDID.txt contain the mapping from gene symbols/names to RefSeq IDs that are used in HPRDNetworkWithComplexes.txt.

## How To Compile

```
g++ createGraph.cpp -o createGraph
g++ extractModules_OptDis.cpp -fopenmp -O3 -g -o extractModules_OptDis
g++ createModules.cpp -o createModules
```

## Usage

First, run createGraph to make a base folder that contain expression profile and network information that OptDis and Dense programs can extract the modules:

```
./createGraph HPRDNetworkWithComplexes.txt TEST.txt TEST.txt TEST HPRDID.txt 61 21
```

The 1st parameter is the network file, 2nd and 3rd parameter are the expression profiles (before I intended to use copy number here, but i leave one for copy number profiles). The 4th parameter is the base folder. The 5th parameter here is the network id file. The 6th parameter is the number of samples in the positive class. The 7th parameter is the number of samples in the negative class.

To extract the connected subnetwork marker as in the ISMB's paper, you run the following command:

```
./extractModules_OptDis TEST 61 4
```

The 1st parameter here is the base folder that was created by the above createGraph program. The 2nd parameter here is the number samples in the positive class. The 3rd parameter here is the maximum size of a subnetwork.

To create a subnetwork activity file:

```
./createModules TEST 61 TEST_OptDis.csv OptDis
```

