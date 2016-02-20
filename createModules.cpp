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
#define MATHLIB_STANDALONE
#include <time.h>

using namespace std;

const int LINELIMIT=1000000;
const int MAXMODULES=100000;
//const int MAXDIMS=300;
const int MAXGENES=30000;
const int MAXSAMPLES=600;
char* DELIM="\t\n";

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
double expressionArr[MAXMODULES];
int numModules=0;

char* token;
char* token1;
double **g;

double** eArr[MAXSAMPLES];

FILE *f;
FILE *f1;
FILE *f2;


int main(int argc,char* argv[])
{
    srand ( time(NULL) );

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
    	g[i][j]=w;
    	g[j][i]=w;
    }
	fclose(f);

    //cout<<"down till here111"<<endl;

    sprintf(filename,"%s/expression.txt",argv[1]);
	f=fopen(filename,"r");
    fgets(line, LINELIMIT, f);
        //cout<<"iii"<<endl;
        token = strtok(line, DELIM);
        //cout<<line<<endl;
        m=-1;
        do
        {
            m++;
            token = strtok(NULL, DELIM);
        } while (token!=NULL);
	fclose(f);
	//cout<<m<<endl;
	//cout<<"down till here222"<<endl;
	//cout<<"NUM DIMENSIONS "<<m<<endl;
	/*
	for (i=0;i<MAXGENES;i++)
	for (j=0;j<m;j++)
        active[i][j]=false;
    */
    for (i=0;i<MAXSAMPLES;i++)
    {
        //cout<<i<<" good"<<endl;
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
        //cout<<line<<endl;
        token = strtok(line, DELIM);
        for (j=0;j<m;j++)
        {
            token = strtok(NULL, DELIM);
            sscanf(token,"%lf",&eArr[j][i][0]);
        }
	}
	fclose(f);

    //cout<<"down till here"<<endl;

   module mod;
   do{


            numModules++;
            sprintf(filename1,"%s/%s/modules/%d",argv[1],argv[4],numModules);
            cout<<filename1<<endl;

            mod.clear();
            try
            {

                //inputFile.open(filename1);
                ifstream inputFile;
                inputFile.exceptions(ifstream::eofbit | ifstream::badbit | ifstream::failbit | ifstream::goodbit);
                inputFile.open(filename1);

                //if (inputFile.is_open())
                //cout<<filename1;
                inputFile.close();
            }
            catch (std::exception& e)
            {
                //cout<<"gets to here"<<endl;
                numModules--;
                break;
            }
            f1=fopen(filename1,"r");
            //cout<<"it comes inside here"<<endl;
        //cout<<"module";
            while (fgets(line1, LINELIMIT, f1)>0){
                token1 = strtok(line1, DELIM);
                mod.insert(g2id[token1]);
                //cout<<"M"<<(int)mod.size()<<endl;
                if (maxSize<((int)mod.size()))
                    maxSize=((int)mod.size());
                //cout<<" "<<token1;
            }
            fclose(f1);

            moduleArr.push_back(mod);
            //numModules++;
    } while (true);

    numModules=moduleArr.size();
    set<int>::iterator iter1;

    f=fopen(argv[3],"w");

    fprintf(f,"CANCER");
    for (i=1;i<=numModules;i++)
    fprintf(f,",m%i",i);
    fprintf(f,"\n");

    for (i=0;i<m;i++)
    {
            if (i<numN)
                fprintf(f,"0");
            else
                fprintf(f,"1");
        for (j=0;j<numModules;j++)
        {
            expressionArr[j]=0.0;
            for (iter1=moduleArr[j].begin();iter1 != moduleArr[j].end();iter1++)
                expressionArr[j]+=(eArr[i][*iter1][0]);
            expressionArr[j]/=moduleArr[j].size()*1.0;
            fprintf(f,",%lf",expressionArr[j]);
        }
        fprintf(f,"\n");
    }
    fclose(f);


    //cout<<"after reading expresion
}
