#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

using namespace std;

const int LINELIMIT=1000000;
const int MAXMODULES=100000;
const int MAXDIMS=600;
const int MAXGENES=30000;
const int MAXSAMPLES=600;
const int MAXEDGES=100000;
const int MINSIZE=3;
const int MAXSIZES=5;
const double OVERLAPFACTOR=2.0;
/*const double MINEDGETHRESHOLD=0.00000001;*/
const double MINEDGETHRESHOLD=0.7;

char* DELIM="\t\n";
const double untouched=-1000000000.0;

typedef set<int> module;
vector< module> moduleArr;

int n;  // the number of genes
int m;  // the number of samples
int numN; //the number of normal samples
//int MAXSIZE=0;
int MAXSIZE=5;
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
double* table[MAXGENES];
int* already[MAXGENES][MAXSIZES];
int* whichVertex[MAXGENES];
int* whichComb[MAXGENES];
//module* moduleTable[MAXGENES];
int color[MAXGENES];
double maxSeperation[MAXGENES];
double seperation[MAXGENES];
double maxSeperationSize[MAXGENES][MAXSIZES];
module optimalModule[MAXGENES];
module optimalModuleSize[MAXGENES][MAXSIZES];
int maxComs;
int numModules=0;
int numEdges;
set<int> upGenes;
set<int> downGenes;

char* token;
char* token1;
double **g;

double** eArr[MAXSAMPLES];

FILE *f;
FILE *f1;
FILE *f2;

void calculateVertexWeights(int nn,int nu)
{
    int i,j,k;

    double inClass=0.0;
    double inClass1=0.0;
    double beClass=0.0;

    for (k=0;k<n;k++){
        inClass=0.0;
        inClass1=0.0;

        for (i=0;i<nn;i++)
		for (j=i+1;j<nn;j++)
		inClass1+=fabs(eArr[i][k][0]-eArr[j][k][0]);
        inClass1/=(nn*(nn-1)*1.0)/2.0;
        inClass+=inClass1;

	inClass1=0.0;
        for (i=nn;i<nn+nu;i++)
        for (j=i+1;j<nn+nu;j++)
            inClass1+=fabs(eArr[i][k][0]-eArr[j][k][0]);
        inClass1/=(nu*(nu-1)*1.0)/2.0;
        inClass+=inClass1;

        inClass/=2.0;

        beClass=0.0;
        for (i=0;i<nn;i++)
        for (j=nn;j<nn+nu;j++)
            beClass+=fabs(eArr[i][k][0]-eArr[j][k][0]);
        beClass/=(nn*nu*1.0);

        if (fabs(inClass)<0.000001)
            seperation[k]=-1.0;
        else
            //weight[k]=totalBeClass/totalInClass+ctotalBeClass/ctotalInClass;
            seperation[k]=beClass-inClass;


        //printf("vertex %d with weight %lf %lf %lf\n",k,weight[k],totalBeClass,totalInClass);
    }
}

double calculateModuleWeights(int nn,int nu,module* mo)
{
    int i,j,k;

    double inClass=0.0;
    double inClass1=0.0;
    double beClass=0.0;
    double ai,aj;

    //for (k=0;k<n;k++){
        inClass=0.0;
        inClass1=0.0;

        set<int>::iterator iter;
        for (i=0;i<nn;i++)
        for (j=i+1;j<nn;j++)
        {
            ai=0.0;
            aj=0.0;
            for (iter=mo->begin();iter!=mo->end();iter++){
                    ai+=eArr[i][*iter][0]/(mo->size()*1.0);
                    aj+=eArr[j][*iter][0]/(mo->size()*1.0);
            }

            inClass1+=fabs(ai-aj);
        }
        inClass1/=(nn*(nn-1)*1.0)/2.0;
        inClass+=inClass1;

        for (i=nn;i<nn+nu;i++)
        for (j=i+1;j<nn+nu;j++)
        {
            ai=0.0;
            aj=0.0;
            for (iter=mo->begin();iter!=mo->end();iter++){
                    ai+=eArr[i][*iter][0]/(mo->size()*1.0);
                    aj+=eArr[j][*iter][0]/(mo->size()*1.0);
            }

            inClass1+=fabs(ai-aj);
        }
            //inClass1+=fabs(eArr[i][k][0]-eArr[j][k][0]);
        inClass1/=(nu*(nu-1)*1.0)/2.0;
        inClass+=inClass1;

        inClass/=2.0;

        beClass=0.0;
        for (i=0;i<nn;i++)
        for (j=nn;j<nn+nu;j++)
        {
            ai=0.0;
            aj=0.0;
            for (iter=mo->begin();iter!=mo->end();iter++){
                    ai+=eArr[i][*iter][0]/(mo->size()*1.0);
                    aj+=eArr[j][*iter][0]/(mo->size()*1.0);
            }

            beClass+=fabs(ai-aj);
        }
            //beClass+=fabs(eArr[i][k][0]-eArr[j][k][0]);
        beClass/=(nn*nu*1.0);

        if (fabs(inClass)<0.000001)
            return -1.0;
        else
            //weight[k]=totalBeClass/totalInClass+ctotalBeClass/ctotalInClass;
            return beClass-inClass;
}


void quickSort(int *a, int lo, int hi)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
    int i=lo, j=hi, h;
    double x=gainArr[a[(lo+hi)/2]];

    while (i<=j)
    {
        while (gainArr[a[i]]>x) i++;
        while (gainArr[a[j]]<x) j--;

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

void initializeTable(int MAXSIZE)
{
    maxComs=1;
    int i;
    for (i=1;i<=MAXSIZE;i++)
        maxComs=maxComs*2;
    maxComs--;

    for (i=0;i<n;i++)
    {
        table[i]=(double*)malloc(sizeof(double)*(maxComs+1));
        whichVertex[i]=(int*)malloc(sizeof(int)*(maxComs+1));
        whichComb[i]=(int*)malloc(sizeof(int)*(maxComs+1));
        //moduleTable[i]=(module*)malloc(sizeof(module)*(maxComs+1));
    }

/*
    for (i=0;i<n;i++)
    for (j=0;j<MAXSIZE;j++)
    {
        already[i][j]=(int*)malloc(sizeof(int)*(maxComs+1));
        //whichVertex[i]=(int*)malloc(sizeof(int)*(maxComs+1));
        //whichComb[i]=(int*)malloc(sizeof(int)*(maxComs+1));
        //moduleTable[i]=(module*)malloc(sizeof(module)*(maxComs+1));
    }
*/

}

void cleanTable(int MAXSIZE)
{
    int i;
    #pragma omp parallel for
    for (i=0;i<n;i++){
        int j;
        for (j=0;j<=maxComs;j++)
            table[i][j]=untouched;
        //moduleTable[i][j].clear();
    }

    //#pragma omp barrier
}

void randomColor()
{
    int i;
    for (i=0;i<MAXGENES;i++)
        color[i]=rand() % MAXSIZE;
}

void putTogether(module *m,module *m1,module *m2){
     m->clear();
     set<int>::iterator iter;
     for (iter=m1->begin();iter!=m1->end();iter++)
        m->insert(*iter);
     for (iter=m2->begin();iter!=m2->end();iter++)
        m->insert(*iter);
}

void transfer(module *m1,module *m2){
     m2->clear();
     set<int>::iterator iter;
     for (iter=m1->begin();iter!=m1->end();iter++)
     {
        m2->insert(*iter);
     }
}

void parse(module *m, int vertex, int comb)
{

    if ((1<<color[vertex])==comb)
        m->insert(vertex);
    else{
        parse(m,whichVertex[vertex][comb],whichComb[vertex][comb]);
        parse(m,vertex,comb^whichComb[vertex][comb]);
    }
}



void extractModules()
{
    calculateVertexWeights(numN,m-numN);
    set<int>::iterator iter;
    initializeTable(MAXSIZE);
    double epsilon=0.01;
    int numColorings=(int)floor((log(1/epsilon)/log(exp(1.0)))*exp(MAXSIZE));
    int i,j,k,i1,i2;

    int ei[MAXEDGES];
    int ej[MAXEDGES];

    int subset[MAXSIZE+1][maxComs+1];
    int subsetWithColor[MAXSIZE+1][MAXSIZE+1][maxComs+1];

    int* containColors1[MAXSIZE][MAXSIZE][MAXSIZE][MAXSIZE];
    int* containColors2[MAXSIZE][MAXSIZE][MAXSIZE][MAXSIZE];

    for (i=0;i<MAXSIZE;i++)
    for (j=0;j<MAXSIZE;j++)
    for (i1=0;i1<MAXSIZE;i1++)
    for (i2=0;i2<MAXSIZE;i2++)
    {
        containColors1[i][j][i1][i2]=(int*)malloc(sizeof(int)*(2*maxComs+1));
        containColors2[i][j][i1][i2]=(int*)malloc(sizeof(int)*(2*maxComs+1));
    }

    for (i=0;i<MAXSIZE;i++)
    for (j=0;j<MAXSIZE;j++)
    for (i1=0;i1<MAXSIZE;i1++)
    for (i2=0;i2<MAXSIZE;i2++)
    for (k=0;k<2*maxComs+1;k++)
    {
                containColors1[i][j][i1][i2][k]=0;
                containColors2[i][j][i1][i2][k]=0;
    }



    for (i=0;i<=MAXSIZE;i++)
        subset[i][0]=0;

    for (i=0;i<=MAXSIZE;i++)
    for (j=0;j<=MAXSIZE;j++)
            subsetWithColor[i][j][0]=0;

    int test;
    int count1;
    int count2;

    for (int i=1;i<=maxComs;i++)
    {
         count1=0;
         for (j=0;j<MAXSIZE;j++)
         {
             test = i & (1 << j);
             if (test>0)
                count1++;
         }
         //cout<<count1<<" "<<i<<endl;
         subset[count1][0]++;
         subset[count1][subset[count1][0]]=i;
         for (j=0;j<MAXSIZE;j++)
         {
             test = i & (1 << j);
             if (test>0)
             {
                //cout<<i<<" XXX "<<j<<" "<<count1<<" "<<subsetWithColor[count1][j][0]<<endl;
                subsetWithColor[count1][j][0]++;
                subsetWithColor[count1][j][subsetWithColor[count1][j][0]]=i;
             }
         }
    }


    int at1[MAXSIZE+1];
    int at2[MAXSIZE+1];
    for (int i=1;i<maxComs;i++)
    for (int k=1;k<maxComs;k++)
    if ((i+k)==(i^k))
    {
         count1=0;
         for (j=0;j<MAXSIZE;j++)
         {
             test = i & (1 << j);
             if (test>0)
             {
                at1[count1]=j;
                count1++;
             }
         }
         count2=0;
         for (j=0;j<MAXSIZE;j++)
         {
             test = k & (1 << j);
             if (test>0)
             {
                at2[count2]=j;
                count2++;
             }
         }
         int i1,i2;
         //cout<<count1<<" "<<i<<endl;
         for (i1=0;i1<count1;i1++)
         for (i2=0;i2<count2;i2++){
             containColors1[at1[i1]][count1-1][at2[i2]][count2-1][0]++;
             containColors1[at1[i1]][count1-1][at2[i2]][count2-1][containColors1[at1[i1]][count1-1][at2[i2]][count2-1][0]]=i;
             containColors2[at1[i1]][count1-1][at2[i2]][count2-1][0]++;
             containColors2[at1[i1]][count1-1][at2[i2]][count2-1][containColors2[at1[i1]][count1-1][at2[i2]][count2-1][0]]=k;
         }
    }


    numEdges=0;
    for (i=0;i<n;i++)
    for (j=i+1;j<n;j++)
    if (g[i][j]>=MINEDGETHRESHOLD){
        ei[numEdges]=i;
        ej[numEdges]=j;
        numEdges++;
    }

    cout<<"NUM EDGES "<<numEdges<<endl;

    //int numColorings=ceil(log(1/epsilon)*exp(MAXSIZE)); //FIX THIS;
    int j1,j2,si,sj,sizei,sizej,total;

    for (i=0;i<n;i++)
        maxSeperation[i]=untouched;

    for (i=0;i<n;i++)
    for (j=0;j<MAXSIZE;j++)
        maxSeperationSize[i][j]=untouched;

    for (i=1;i<=numColorings;i++){
        randomColor();
        cleanTable(MAXSIZE);

        //initialize the table
        for (k=0;k<n;k++)
        {
            //cout<<k<<" "<<color[k]<<" "<<(1<<color[k])<<endl;
            //table[k][1>>color[k]]=seperation[k];
            table[k][1<<color[k]]=seperation[k];
            //cout<<seperation[k]<<" "<<k<<endl;

            /*
            already[k][0][0]=1;
            already[k][0][1]=1<<color[k];
            */

        }

        for (j=2;j<=MAXSIZE;j++){
            //#pragma omp parallel for
            for (k=0;k<numEdges;k++){
                //cout<<k<<" "<<numEdges<<"EDGE"<<endl;
                int sizei,sizej,si,sj,total,i1;

                /*
                int th_id = omp_get_thread_num();
                std::ostringstream ss;
                ss << "Hello World from thread " << th_id <<std::endl;
                std::cout << ss.str();
                */

                for (sizei=1;sizei<j;sizei++){
                    sizej=j-sizei;
                    for (i1=1;i1<=containColors1[color[ei[k]]][sizei-1][color[ej[k]]][sizej-1][0];i1++)
                    {
                        si=containColors1[color[ei[k]]][sizei-1][color[ej[k]]][sizej-1][i1];
                        sj=containColors2[color[ei[k]]][sizei-1][color[ej[k]]][sizej-1][i1];
                        //cout<<color[ei[k]]<<" "<<si<<" "<<sizei<<" "<<color[ej[k]]<<" "<<sj<<" "<<sizej<<endl;
                        if ((table[ei[k]][si]>untouched)&&(table[ej[k]][sj]>untouched)){
                                   total=si+sj;
                                   if (table[ei[k]][total]<table[ei[k]][si]+table[ej[k]][sj]){
                                                table[ei[k]][total]=table[ei[k]][si]+table[ej[k]][sj];
                                                whichVertex[ei[k]][total]=ej[k];
                                                whichComb[ei[k]][total]=sj;

                                                //cout<<ei[k]<<" "<<total<<" "<<table[ei[k]][total]<<endl;
                                    }

                                    if (table[ej[k]][total]<table[ei[k]][si]+table[ej[k]][sj]){
                                                table[ej[k]][total]=table[ei[k]][si]+table[ej[k]][sj];
                                                whichVertex[ej[k]][total]=ei[k];
                                                whichComb[ej[k]][total]=si;
                                    }
                        }
                    }

                }
            }
            //#pragma omp barrier
        }

        for (j=0;j<n;j++)
        {
            int size;
            for (size=MINSIZE;size<=MAXSIZE;size++)
            for (k=1;k<=subset[size][0];k++)
            if (table[j][subset[size][k]]>untouched)
            if (maxSeperationSize[j][size]<table[j][subset[size][k]]/(size*1.0))
            {
                maxSeperationSize[j][size]=table[j][subset[size][k]]/(size*1.0);
                optimalModuleSize[j][size].clear();
                parse(&optimalModuleSize[j][size],j,subset[size][k]);
            }
        }
    }

    cout<<"Almost done"<<endl;
    numModules=0;


    double saveMax,save;
    int saveSize,size;
    for (j=0;j<n;j++)
    {
        saveMax=untouched;
        for (size=MINSIZE;size<=MAXSIZE;size++)
        if (maxSeperationSize[j][size]>untouched)
        {
            save=calculateModuleWeights(numN,m-numN,&optimalModuleSize[j][size]);
            if (save>saveMax){
                saveMax=save;
                saveSize=size;
            }
        }
        if (saveMax>untouched){
                //cout<<"sizeMax "<<sizeMax<<endl;
                gainArr[numModules]=saveMax;
                moduleOrder[numModules]=numModules;
                moduleArr.push_back(optimalModuleSize[j][saveSize]);
                numModules++;
        }
    }
    quickSort(moduleOrder,0,numModules-1);
}

int main(int argc,char* argv[])
{
    srand ( time(NULL) );

    sprintf(cmdline,"mkdir %s/OptDis/",argv[1]);
    system(cmdline);

    sprintf(cmdline,"mkdir %s/OptDis/modules/",argv[1]);
    system(cmdline);

    sscanf(argv[3],"%d",&MAXSIZE);
    sscanf(argv[2],"%d",&numN);
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

    g=(double**)malloc(sizeof(double*)*n);
    int i,j;
    for (i=0;i<n;i++)
        g[i]=(double*)malloc(sizeof(double)*n);
    for (i=0;i<n;i++)
    for (j=0;j<n;j++)
    	g[i][j]=0.0;
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
	if (w>=MINEDGETHRESHOLD){
    	g[i][j]=w;
    	g[j][i]=w;
	}
    }
	fclose(f);

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

    //greedyExpand();

    double avg=0.0;
    for (i=0;i<n;i++){
        avg=0.0;
        for (j=numN;j<m;j++)
            avg+=eArr[j][i][0];
        avg/=((m-numN)*1.0);
        if (avg>0.0)
            upGenes.insert(i);
        else
            downGenes.insert(i);
    }
    extractModules();

    vector< module >::iterator iter;
    set<int>::iterator iter1;
    set<int>::iterator iter2;

    sprintf(filename,"%s/OptDis/modules.txt",argv[1]);
    f=fopen(filename,"w");
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

    top50=0;
    visited.clear();
    int newGenes=0;
    for (i=0;i<numModules;i++){
        newGenes=0;
        for (iter1=(moduleArr[moduleOrder[i]]).begin();iter1 != (moduleArr[moduleOrder[i]]).end();iter1++){
                if (visited.find(*iter1)==visited.end())
            newGenes++;
        }

        if (OVERLAPFACTOR*newGenes>=(moduleArr[moduleOrder[i]]).size())
        //if ((gainArr[moduleOrder[i]]<coExArr[moduleOrder[i]])||(gainArr[moduleOrder[i]]<=1.15*coExArr[moduleOrder[i]]))
        //if (gainArr[moduleOrder[i]]>1.5*coExArr[moduleOrder[i]])
        {
            top50++;
            chosen[i]=true;
            //cout<<i<<" "<<top50<<endl;
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
				fprintf(f,"%d(%lf) ",*iter1,seperation[*iter1]);
			}
			fprintf(f,"\n");
			fprintf(f,"\n");
			topSN++;
		}

	}
	fclose(f);

    sprintf(cmdline,"rm %s/OptDis/modules/*",argv[1]);
    cout<<cmdline<<endl;
    system(cmdline);

    cout<<numModules<<" top50 "<<top50<<endl;
    int curr=0;
    for (i=0;i<numModules;i++)
    if (chosen[i])
    {
        curr++;
        sprintf(filename,"%s/OptDis/modules/%d",argv[1],curr);
        //cout<<filename<<endl;
        f=fopen(filename,"w");
        for (iter1=(moduleArr[moduleOrder[i]]).begin();iter1 != (moduleArr[moduleOrder[i]]).end();iter1++)
            //fprintf(f,"%d\n",*iter1);
	{
            //fprintf(f,"%s\n",id2g[(*iter1)].c_str());
		fprintf(f,"%s\t",id2g[(*iter1)].c_str());
                if (upGenes.find(*iter1)!=upGenes.end())
                    fprintf(f,"1\n");
                else
                    fprintf(f,"-1\n");
	}
            fclose(f);
    }

    for (i=0;i<n;i++)
        free(g[i]);
    free(g);
}
