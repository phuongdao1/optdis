## OPTIMALLY DISCRIMINATIVE SUBNETWORK MARKERS

OptDis (OPTIMALLY DISCRIMINATIVE SUBNETWORK MARKERS) is a tool for computing the provably optimal subnetworks for classification of samples from different classes. The discriminative score is calculated as the difference between the total distance between samples from different classes and the total distance between samples from the same class. Our algorithm is based on [color-coding paradigm] (#colorcoding), which allows for identifying the optimally discriminative subnetwork markers for any given error probability. The implementation here is based on [our paper](#citation) where we fix the error probability at 0.01.

## Sample Input Files

The expression profile of the sample data is in the file TEST.txt. It has 61 positive samples and 21 negative samples. Note that you need to put all the samples from a class together. all the samples in positive class should be ahead of the samples from negative class. We also give your the PPI network from [HPRD database](http://www.hprd.org) in the file HPRDNetworkWithComplexes.txt. The file HPRDID.txt contain the mapping from gene symbols/names to RefSeq IDs that are used in HPRDNetworkWithComplexes.txt.

## How To Compile

```
g++ createGraph.cpp -o createGraph
g++ extractModules_OptDis.cpp -O3 -g -o extractModules_OptDis
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

So the 1st parameter here is the base folder that was created by the above createGraph program. The 2nd parameter here is the number samples in the positive class. The 3rd parameter here is output subnetwork activity file. The 4th parameter is always OptDis for connected subnetworks. The activity of a subnetwork in a sample is the average of the expressions of genes from that subnetwork.


## CITATION

<a name="citation"></a>
1. Phuong Dao, Kendric Wang, Colin Collins, Martin Ester, Anna Lapuk and S. Cenk Sahinalp. "Optimally discriminative subnetwork markers predict response to chemotherapy". Bioinformatics. 2011 Jul 1; 27(13): i205–i213. [[LINK]](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3117373/)

<a name="colorcoding"></a>
2. Noga Alon, Raphael Yuster, Uri Zwick. "Color-coding". Journal of the ACM (JACM) 42 (4), 844-856 (1995).
[[LINK]](http://dl.acm.org/citation.cfm?id=210337)

## CONTACTS

Please report any problems directly to the github issue tracker. Also, you can also send your feedbacks to phuongdao1@gmail.com.

## LICENSE

OptDis is distributed under GNU GPL license version 3.
