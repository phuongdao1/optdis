#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <deque>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

using namespace std;

const int LINELIMIT=16384;
const int MAXMODULES=30000;
const int MAXDIMS=500;
const int MAXGENES=20000;
const int MAXSAMPLES=500;
char* DELIM="\t\n";
const double OVERLAPFACTOR=2.0;
const int MINSIZE=3;

typedef set<int> module;
vector< module> moduleArr;

int n;  // the number of genes
int m;  // the number of samples
int numN; //the number of normal samples
int maxSize=0;
int i,j,k,l;
map<string, int> g2id;
map<int, string> id2g;
set<int> visited;
char cmdline[256];
char filename[256];
char filename1[256];
char line[LINELIMIT];
char line1[LINELIMIT];
double gainArr[MAXMODULES];
int moduleOrder[MAXMODULES];
int numModules=0;
double denseThreshold;

char* token;
char* token1;
double **g;

double** eArr[MAXSAMPLES];

FILE *f;
FILE *f1;
FILE *f2;

double informationGain(double* v)
{
    int i,j,k;
    int numBins=(int)ceil(log(m*1.0));
    //int numBins=2;
    double l=1000000.0;
    double r=-1000000.0;
    for (i=0;i<m;i++){
        if (v[i]<l)
            l=v[i];
        if (v[i]>r)
            r=v[i];
    }

    double binLength=(r-l)/(numBins*1.0);
    double prob0=numN/(m*1.0);
    double prob1=(m-numN)/(m*1.0);
    double entropy=-prob1*log(prob1)-prob0*log(prob0);

    double maxR;
    double minThres=0.00001;
    int totalInBin;
    double pInBin;
    double prob1InBin;
    double prob0InBin;
    int num1InBin;
    int num0InBin;
    for (int i=1;i<=numBins;i++)
    {
        if (i==numBins)
            maxR=r+minThres;
        else
            maxR=l+i*binLength;

        totalInBin=0;
        for (j=0;j<m;j++)
        if ((v[j]>=l+(i-1)*binLength)&&(v[j]<maxR))
            totalInBin++;
        pInBin=totalInBin/(m*1.0);

        if (totalInBin>0)
        {
            num0InBin=0;
            for (j=0;j<numN;j++)
            if ((v[j]>=l+(i-1)*binLength)&&(v[j]<maxR))
                num0InBin++;
            prob0InBin=num0InBin/(totalInBin*1.0);

            num1InBin=0;
            for (j=numN;j<m;j++)
            if ((v[j]>=l+(i-1)*binLength)&&(v[j]<maxR))
                num1InBin++;
            prob1InBin=num1InBin/(totalInBin*1.0);

            if (prob1InBin>minThres)
                entropy=entropy+pInBin*prob1InBin*log(prob1InBin);
            if (prob0InBin>minThres)
                entropy=entropy+pInBin*prob0InBin*log(prob0InBin);
        }
    }
    return entropy;
}

void quickSort(int *a, int lo, int hi)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
    int i=lo, j=hi, h;
    double x=gainArr[a[(lo+hi)/2]];

    //  partition
    //cout<<"fail where"<<endl;
    while (i<=j)
    {
        while (gainArr[a[i]]>x) i++;
        while (gainArr[a[j]]<x) j--;
        //cout<<"fail here"<<endl;
        if (i<=j)
        {
            h=a[i]; a[i]=a[j]; a[j]=h;
            i++; j--;
        }
    }

    //  recursion
    if (lo<j) quickSort(a, lo, j);
    if (i<hi) quickSort(a, i, hi);
}

void greedyExpand()
{
    double arr[MAXSAMPLES];
    double arr1[MAXSAMPLES];
    double ig[MAXGENES];
    int i,j,k,l;
    double gain,maxGain,currentGain;
    int save;
    module current;
    double numEdges,newEdges;
    set<int> neighbors;
    set<int>::iterator iter;
    set<int>::iterator iter1;
    for (i=0;i<n;i++)
    {
        numEdges=0.0;
        current.clear();
        neighbors.clear();

        //cout<<"start "<<i<<endl;

        for (j=0;j<m;j++)
                arr[j]=0.0;
        //
        l=i;
        do{
            current.insert(l);

            for (j=0;j<m;j++)
                arr[j]+=eArr[j][l][0];



            for (j=0;j<n;j++)
            if (g[l][j]>0.0)
            if (current.find(j)==current.end())
                neighbors.insert(j);

            //cout<<"Size "<<neighbors.size()<<endl;

            // calculate current gain
            for (j=0;j<m;j++)
                arr1[j]=arr[j]/(sqrt(current.size()*1.0));
            currentGain=informationGain(arr1);
            //cout<<"current Gain "<<currentGain<<endl;

            l=-1;
            maxGain=-1000001.0;
            for (iter=neighbors.begin();iter!=neighbors.end();iter++){
                //calculate new gain
                newEdges=0.0;
                for (j=0;j<m;j++)
                    arr1[j]=(arr[j]+eArr[j][*iter][0])/(sqrt((current.size()+1)*1.0));
                gain=informationGain(arr1);
                //cout<<*iter<<" gains "<<gain<<endl;
                for (iter1=current.begin();iter1!=current.end();iter1++)
                    newEdges+=g[*iter][*iter1];


                if ((gain>maxGain)&&(gain>currentGain)&&((numEdges+newEdges/(current.size()*(current.size()-1)/2.0))>=denseThreshold)){
                    maxGain=gain;
                    l=*iter;
                }
            }
            //cout<<"choose "<<l<<endl;
        } while ((l>=0)&&(neighbors.size()>0));
        if (current.size()>=3)
        {
            moduleArr.push_back(current);
            cout<<"MODULE START"<<endl;
            for (iter=current.begin();iter!=current.end();iter++){
                cout<<*iter<<endl;
                visited.insert(*iter);
            }
            cout<<"MODULE END"<<endl;
            gainArr[numModules]=currentGain;
            numModules++;
        }
    }
    for (i=0;i<numModules;i++)
        moduleOrder[i]=i;
    quickSort(moduleOrder,0,numModules-1);
}

int main(int argc,char* argv[])
{
    srand ( time(NULL) );

    sscanf(argv[2],"%d",&numN);
    sscanf(argv[3],"%lf",&denseThreshold);
    double w;
    sprintf(filename,"%s/nodes.txt",argv[1]);
	n=0;
	f=fopen(filename,"r");
    while (fgets(line, LINELIMIT, f)>0){
        token = strtok(line, DELIM);
        sscanf(token,"%d",&i);
        if (i>n)
            n=i;
    	token = strtok(NULL, DELIM);
    	id2g[i]=token;
    	g2id[token]=i;
    }
	fclose(f);
	n++;

	//cout<<"before anything "<<n<<endl;

    g=(double**)malloc(sizeof(double*)*n);
    int i,j;
    for (i=0;i<n;i++)
    {
        g[i]=(double*)malloc(sizeof(double)*n);
        cout<<i<<endl;
    }
    for (i=0;i<n;i++)
    for (j=0;j<n;j++)
    	g[i][j]=0.0;

    //cout<<"after anything"<<endl;

    sprintf(filename,"%s/edges.txt",argv[1]);
	f=fopen(filename,"r");
    while (fgets(line, LINELIMIT, f)>0){
        //cout<<"iii"<<endl;
        token = strtok(line, DELIM);
        sscanf(token,"%d",&i);
    	token = strtok(NULL, DELIM);
    	sscanf(token,"%d",&j);
    	token = strtok(NULL, DELIM);
    	sscanf(token,"%lf",&w);
    	//totalEdges+=1.0;
    	//expTotalEdges+=w;
    	g[i][j]=w;
    	g[j][i]=w;
    }
	fclose(f);

	//cout<<"reading expression file"<<endl;

    sprintf(filename,"%s/expression.txt",argv[1]);
	f=fopen(filename,"r");
    fgets(line, LINELIMIT, f);
        //cout<<"iii"<<endl;
        token = strtok(line, DELIM);
        m=-1;
        do
        {
            m++;
            token = strtok(NULL, DELIM);
        } while (token!=NULL);
	fclose(f);
	cout<<"NUM DIMENSIONS "<<m<<endl;
	/*
	for (i=0;i<MAXGENES;i++)
	for (j=0;j<m;j++)
        active[i][j]=false;
    */
    for (i=0;i<MAXSAMPLES;i++)
    {
        eArr[i]=(double**)malloc(sizeof(double*)*n);

        for (j=0;j<n;j++)
        	eArr[i][j]=(double*)malloc(sizeof(double)*2);

        for (j=0;j<n;j++)
        for (k=0;k<2;k++)
        	eArr[i][j][k]=0.0;
    }

    //cout<<"before reading expression"<<endl;
    sprintf(filename,"%s/expression.txt",argv[1]);
	f=fopen(filename,"r");
	for (i=0;i<n;i++)
	{
        fgets(line, LINELIMIT, f);
        token = strtok(line, DELIM);
        for (j=0;j<m;j++)
        {
            token = strtok(NULL, DELIM);
            sscanf(token,"%lf",&eArr[j][i][0]);
        }
	}
	fclose(f);

    //cout<<"after reading expresion

    //cout<<"done reading"<<endl;

    greedyExpand();

    vector< module >::iterator iter;
    set<int>::iterator iter1;
    set<int>::iterator iter2;

    int top50;
    /*
    if (numModules<=50)
        top50=numModules;
    else
        top50=50;
    */

    bool chosen[MAXMODULES];

    for (i=0;i<numModules;i++){
	chosen[i]=false;
    }


    sprintf(filename,"%s/Dense/modules.txt",argv[1]);
    f=fopen(filename,"w");
    //fprintf(f,"Final number of patterns:%d\n\n",top50);


    top50=0;
    visited.clear();
    int newGenes=0;
    for (i=0;i<numModules;i++){
        newGenes=0;
        for (iter1=(moduleArr[moduleOrder[i]]).begin();iter1 != (moduleArr[moduleOrder[i]]).end();iter1++){
                if (visited.find(*iter1)==visited.end())
            newGenes++;
        }
        //cout<<"New Genes "<<newGenes<<" "<<endl;

        if (OVERLAPFACTOR*newGenes>=(moduleArr[moduleOrder[i]]).size())
        //if ((gainArr[moduleOrder[i]]<coExArr[moduleOrder[i]])||(gainArr[moduleOrder[i]]<=1.15*coExArr[moduleOrder[i]]))
        //if (gainArr[moduleOrder[i]]>1.5*coExArr[moduleOrder[i]])
        {
            top50++;
            chosen[i]=true;
            //cout<<i<<" "<<top50<<endl;
            //cout<<"New Module "<<newGenes<<" "<<i<<endl;
            for (iter1=(moduleArr[moduleOrder[i]]).begin();iter1 != (moduleArr[moduleOrder[i]]).end();iter1++){
                        visited.insert(*iter1);
                        //cout<<*iter1<<endl;
            }
        }

        if (top50==50)
            break;
    }


    //top50=numModules;

    fprintf(f,"Final number of patterns:%d\n\n",top50);
	
	int topSN = 1;
	for (i=0;i<numModules;i++) {
		if (chosen[i]) {
			fprintf(f,"Module %d\n", topSN);
			fprintf(f,"\t Size:%d\tDiscriminative Score:%lf\n",moduleArr[moduleOrder[i]].size(),gainArr[moduleOrder[i]]);
			fprintf(f,"\t");
			for (iter1=(moduleArr[moduleOrder[i]]).begin();iter1 != (moduleArr[moduleOrder[i]]).end();iter1++){
				fprintf(f,"%d ",*iter1);
			}
			fprintf(f,"\n");
			fprintf(f,"\n");
			topSN++;
		}

	}
	fclose(f);

    //cout<<"stupid at the end"<<endl;

    sprintf(cmdline,"rm %s/Dense/modules/*",argv[1]);
    cout<<cmdline<<endl;
    system(cmdline);

    //cout<<"stupid at the end"<<endl;
    cout<<numModules<<" top50 "<<top50<<endl;
    int curr=0;
    for (i=0;i<numModules;i++)
    if (chosen[i])
    {
        curr++;
        sprintf(filename,"%s/Dense/modules/%d",argv[1],curr);
        //cout<<filename<<endl;
        f=fopen(filename,"w");
        for (iter1=(moduleArr[moduleOrder[i]]).begin();iter1 != (moduleArr[moduleOrder[i]]).end();iter1++)
            //fprintf(f,"%d\n",*iter1);
	{
            //fprintf(f,"%s\n",id2g[(*iter1)].c_str());
		fprintf(f,"%s\n",id2g[(*iter1)].c_str());
                //if (upGenes.find(*iter1)!=upGenes.end())
                    //fprintf(f,"1\n");
                //else
                    //fprintf(f,"-1\n");
	}
            fclose(f);
    }


    for (i=0;i<n;i++)
        free(g[i]);
    free(g);
}
