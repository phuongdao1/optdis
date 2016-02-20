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
double NULLMARKER=-999999.1;
const int MAXSAMPLES=600;
const int MAXGENES=30000;
const char* DELIM="\t\n";
const char* DELIM1="|";
const char* DELIM2=" \t\n ";
const char* DELIM3=",";

int n;  // the number of genes
int m;  // the number of dimensions
int CHECK=448;
int nn;
int nu;
int i,j,k,l;
map<string, int> g2id;
map<string, int> stringMap;
map<string, int> g2id2;
map<string, int> g2id1;
map<int, string> id2g;
char* line;
char cmdline[256];
char filename[256];
char command[256];
char line1[256];
bool* expressed[MAXSAMPLES];
double** eArr[MAXSAMPLES];
double** cArr[MAXSAMPLES];
//double** eArr1[MAXSAMPLES];
double meanArr[MAXGENES];
double stdArr[MAXGENES];
double cmeanArr[MAXGENES];
double cstdArr[MAXGENES];
char* tempString1;
char* tempString2;
char* token;
char* token1;
int* present;
double **ppi;
bool allPresent[MAXGENES];
//bool passFiltered[MAXGENES];
double thres=0.025;
//double thres1=30.0;

FILE *f;
FILE *f1;
bool new1,new2;
int id1,id2;

void allocateGraph(double*** a,int n)
{
    *a=(double**)malloc(sizeof(double*)*n);
    int i,j;
    for (i=0;i<n;i++)
        (*a)[i]=(double*)malloc(sizeof(double)*n);
    for (i=0;i<n;i++)
    for (j=0;j<n;j++)
    	(*a)[i][j]=0.0;
}

void deallocateGraph(double*** a,int n)
{
    //*a=(double**)malloc(sizeof(double*)*n);
    int i;
    for (i=0;i<n;i++)
        free((*a)[i]);
    free(*a);
}

void resetExpressionArray()
{
	int i,j;
	for (i=0;i<n;i++)
		present[i]=0;
}


void processAllSamples()
{
	int i,j,k,ck,l;
	for (i=0;i<n;i++)
		if (allPresent[i])
		{
		k=0;
		ck=0;
	
		for (j=0;j<nn+nu;j++)
		if (eArr[j][i][0]>(NULLMARKER+1.0))
		k++;

		if (k>0) {
			for (j=0;j<nn+nu;j++)
			if (eArr[j][i][0]>(NULLMARKER+1.0))
				meanArr[i]+=eArr[j][i][0];
		
			meanArr[i]=0.0;
			for (j=0;j<nn+nu;j++)
			if (eArr[j][i][0]>(NULLMARKER+1.0))
				meanArr[i]+=eArr[j][i][0];
			meanArr[i]/=(k)*1.0;
		
			double var=0.0;
			for (j=0;j<nn+nu;j++)
			if (eArr[j][i][0]>(NULLMARKER+1.0))
			var+=(eArr[j][i][0]-meanArr[i])*(eArr[j][i][0]-meanArr[i]);
			var/=(k)*1.0;
			stdArr[i]=sqrt(var);
		}
	
		for (j=0;j<nn+nu;j++)
		if (cArr[j][i][0]>(NULLMARKER+1.0))
		ck++;

		if (ck>0) {
			for (j=0;j<nn+nu;j++)
			if (cArr[j][i][0]>(NULLMARKER+1.0))
				cmeanArr[i]+=cArr[j][i][0];
		
			cmeanArr[i]=0.0;
			for (j=0;j<nn+nu;j++)
			if (cArr[j][i][0]>(NULLMARKER+1.0))
				cmeanArr[i]+=cArr[j][i][0];
			cmeanArr[i]/=(ck)*1.0;
		
			double var=0.0;
			for (j=0;j<nn+nu;j++)
			if (cArr[j][i][0]>(NULLMARKER+1.0))
				var+=(cArr[j][i][0]-cmeanArr[i])*(cArr[j][i][0]-cmeanArr[i]);
			var/=(ck)*1.0;
			cstdArr[i]=sqrt(var);
		}


	}

	for (i=0;i<n;i++) {
		if (allPresent[i])
			for (j=0;j<nn+nu;j++)
				if (eArr[j][i][0]>(NULLMARKER+1.0)) {
					if (fabs(stdArr[i])<0.0001)
						eArr[j][i][0]=(eArr[j][i][0]-meanArr[i]);
					else
						eArr[j][i][0]=(eArr[j][i][0]-meanArr[i])/stdArr[i];
				}
	}
	
	for (i=0;i<n;i++) {
		if (allPresent[i])
			for (j=0;j<nn+nu;j++)
				if (cArr[j][i][0]>(NULLMARKER+1.0)) {
					if (fabs(cstdArr[i])<0.0001)
						cArr[j][i][0]=(cArr[j][i][0]-cmeanArr[i]);
					else
						cArr[j][i][0]=(cArr[j][i][0]-cmeanArr[i])/cstdArr[i];
				}
	}
}

int compare (const void * a, const void * b)
{
  if ((*(double*)a - *(double*)b)>=0.0)
    return 1;
  return -1;
}



void processSample(int index)
{
	int i,j;


    double a1[MAXGENES+1];
    int n11=0;

    for (i=0;i<n;i++)
    if (allPresent[i]>0)
    if (eArr[index][i][0]>NULLMARKER+1.0)
    {
        n11++;
        a1[n11-1]=eArr[index][i][0];
    }


    qsort(a1,n11,sizeof(double),compare);
    int pivot1=(int)round(n11*(1-thres));
    int pivot2=(int)round(n11*thres);

    //cout<<"pivot"<<" "<<a1[pivot1]<<endl;
    double p1;
    int numDEs=0;
    	for (i=0;i<n;i++)
    	if (allPresent[i])
    	{

                if (eArr[index][i][0]>a1[pivot1])
                    expressed[index][i]=true;
                if ((eArr[index][i][0]<a1[pivot2])&&(fabs(eArr[index][i][0])>0.0001))
                    expressed[index][i]=true;


                if (expressed[index][i])
                    numDEs++;
    	}
    	cout<<index<<" Number of differentially expressed genes: "<<n11<<" "<<numDEs<<endl;
}

void processGene(int index)
{
	int i,j;

    if (allPresent[index]>0)
    {
        double a1[MAXSAMPLES+1];
        int n11=0;

        for (i=0;i<m;i++)
        if (eArr[i][index][0]>NULLMARKER+1.0)
        {
            n11++;
            a1[n11-1]=eArr[i][index][0];
        }


        qsort(a1,n11,sizeof(double),compare);
        int pivot1=(int)floor(n11*(1-thres));
        //int pivot2=(int)round(n11*thres);

        //cout<<"pivot"<<" "<<a1[pivot1]<<endl;
        //double p1;
        int numDEs=0;
        for (i=0;i<m;i++){

                    if (eArr[i][index][0]>=a1[pivot1])
                        expressed[i][index]=true;
                    //if ((eArr[index][i][0]<a1[pivot2])&&(fabs(eArr[index][i][0])>0.0001))
                        //expressed[index][i]=true;


                    //if (expressed[index][i])
                        //numDEs++;
        }
            //cout<<index<<" Number of differentially expressed genes: "<<n11<<" "<<numDEs<<endl;
    }
}

double max(double a,double b)
{
	if (a>b)
		return a;
	return b;
}

double informationGain(double* v,int id)
{
    int i,j,k,m;
    m=nn+nu;

    int numBins=4;
    double l=1000000.0;
    double r=-1000000.0;
    for (i=0;i<m;i++){
        if (v[i]<l)
            l=v[i];
        if (v[i]>r)
            r=v[i];
    }


    double binLength=(r-l)/(numBins*1.0);
    double prob0=nn/(m*1.0);
    double prob1=nu/(m*1.0);
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
            for (j=0;j<nn;j++)
            if ((v[j]>=l+(i-1)*binLength)&&(v[j]<maxR))
                num0InBin++;
            prob0InBin=num0InBin/(totalInBin*1.0);

            num1InBin=0;
            for (j=nn;j<m;j++)
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

double variance(double* v)
{
    int i,j,k,m;
    m=nn+nu;
    double mean=0.0;
    double variance=0.0;

    for (i=0;i<m;i++)
        mean+=v[i]/(m*1.0);

    for (i=0;i<m;i++)
        variance+=(v[i]-mean)*(v[i]-mean)/(m*1.0);

    return variance;
}

int main(int argc,char* argv[])
{
    sprintf(command,"mkdir %s",argv[4]);
    system(command);
    sprintf(command,"mkdir %s/Dense",argv[4]);
    system(command);
    sprintf(command,"mkdir %s/Dense/modules",argv[4]);
    system(command);
    sprintf(command,"mkdir %s/OptDis",argv[4]);
    system(command);
    sprintf(command,"mkdir %s/OptDis/modules",argv[4]);
    system(command);
    

    srand ( time(NULL) );
    //cout<<"norm "<<Rf_pnorm5(1.95,0.0,1.0,0,0)<<endl;
    //sscanf(argv[3],"%lf",&thres);
    n=0;
    sscanf(argv[6],"%d",&nn);
    sscanf(argv[7],"%d",&nu);
     m=nn+nu;
    line=(char*)malloc(sizeof(char)*LINELIMIT);
    tempString1=(char*)malloc(sizeof(char)*LINELIMIT);
    tempString2=(char*)malloc(sizeof(char)*LINELIMIT);
    present=(int*)malloc(sizeof(int)*MAXGENES);

    //cout<<"error down there"<<endl;

    // READING STRING MAP FIRST FOR IDS
    f=fopen(argv[5],"r");
    while (fgets(line, LINELIMIT, f)>0){
    	token = strtok(line, DELIM2);
        if (g2id[token]==0){
    		//new1=true;
    		n++;
    		g2id[token]=n;
    		id2g[n]=token;
    	}
    	//cout<<token<<"id"<<endl;
    	//cout<<token<<endl;
    	i=g2id[token];
    	//cout<<i<<endl;
    	token = strtok(NULL, DELIM2);
    	//cout<<token<<endl;
    	stringMap[token]=i;
    }
    fclose(f);

    //cout<<"error up there"<<endl;

    allocateGraph(&ppi,n);
    //allocateGraph(&diffGraph,n);

    cout<<"m = "<<m<<endl;
    for (i=0;i<MAXSAMPLES;i++)
    {
        eArr[i]=(double**)malloc(sizeof(double*)*n);

        for (j=0;j<n;j++)
        	eArr[i][j]=(double*)malloc(sizeof(double)*2);

        for (j=0;j<n;j++)
        	eArr[i][j][0]=NULLMARKER;
    }
    for (i=0;i<MAXSAMPLES;i++)
    {
        cArr[i]=(double**)malloc(sizeof(double*)*n);

        for (j=0;j<n;j++)
        	cArr[i][j]=(double*)malloc(sizeof(double)*2);

        for (j=0;j<n;j++)
        	cArr[i][j][0]=NULLMARKER;
    }
    /*
    for (i=0;i<MAXSAMPLES;i++)
    {
        deArr[i]=(bool**)malloc(sizeof(bool*)*n);

        for (j=0;j<n;j++)
        	deArr1[i][j]=(bool*)malloc(sizeof(bool)*m);

        for (j=0;j<n;j++)
        for (k=0;k<m;k++)
        	deArr1[i][j][k]=false;
    }
    */


    for (i=0;i<MAXSAMPLES;i++)
    {
        expressed[i]=(bool*)malloc(sizeof(bool)*n);

        for (j=0;j<n;j++)
        	expressed[i][j]=false;
    }

    /*
    for (i=0;i<MAXSAMPLES;i++)
    {
        graph[i]=(double**)malloc(sizeof(double*)*n);
        for (j=0;j<n;j++)
        	graph[i][j]=(double*)malloc(sizeof(double)*n);
        for (k=0;k<n;k++)
        	graph[i][j][k]=0.0;

    }
    */

    //cout<<"reach here"<<endl;
    f=fopen(argv[1],"r");
    //fgets(line, LINELIMIT, f);
    double weight;
    while (fgets(line, LINELIMIT, f)>0){
    	//h1=&line[0];
    	//t=&h1;
    	
// DEBUG: Outputs all PPI edges with corresponding weights    	
//     	cout<<line;
    	
    	token = strtok(line, DELIM2);
    	//token = strtok(NULL, DELIM2);
    	id1=stringMap[token];
    	//token = strtok(NULL, DELIM);
    	//token = strtok(NULL, DELIM);
    	//id1=g2id[token];
    	//token = strtok(NULL, DELIM2);
    	token = strtok(NULL, DELIM2);
    	id2=stringMap[token];
    	token = strtok(NULL, DELIM2);
    	sscanf(token,"%lf",&weight);
    	//cout<<id1<<" id "<<id2<<endl;
    	//cout<<id1<<" "<<id2<<" "<<weight<<endl;
    	if ((id1>0)&&(id2>0)){
    	//if ((ppi[id1-1][id2-1]<0.0001)||(ppi[id1-1][id2-1]>weight)){
    	    //cout<<argv[1]<<" "<<weight<<endl;
            ppi[id1-1][id2-1]=weight/1000.0;
            ppi[id2-1][id1-1]=weight/1000.0;

    	}
    }
    fclose(f);
    //allocateGraph(&posGraph,n);
    //allocateGraph(&negGraph,n);


    //////////////////// READ GENE EXPRESSION FILE
    f=fopen(argv[2],"r");
    int index;
    double v;
    fgets(line, LINELIMIT, f);
    fgets(line, LINELIMIT, f);
    double temp[MAXGENES];
    double bestIG[MAXGENES];
    int savedProbe[MAXGENES];

    for (i=0;i<MAXGENES;i++){
        bestIG[i]=-1000000.0;
    }

    int probeSetID;

    while (fgets(line, LINELIMIT, f)>0)
	if (strlen(line)>2*(nn+nu))
	{
            //cout<<line<<endl;
            //token = strtok(line, DELIM);
            //sscanf(token,"%d",&probeSetID);
	    probeSetID=0;
    		token = strtok(line, DELIM);
    		if (g2id.find(token)!=g2id.end()){

    			id1=g2id[token]-1;

    			present[id1]=1;
    			allPresent[id1]=true;


    			for (index=0;index<nn+nu;index++)
    			{
    			    token = strtok(NULL, DELIM);
                    sscanf(token,"%lf",&v);
                    temp[index]=v;
    			}

    			if (variance(temp)>bestIG[id1])
    			{
                        bestIG[id1]=variance(temp);
                        savedProbe[id1]=probeSetID;
                        for (index=0;index<nn+nu;index++)
                        //if (eArr[index][id1][0]<(NULLMARKER+1.0))
                            //eArr[index][id1][0]=v;
                        //else
                            eArr[index][id1][0]=temp[index];
    			}


    		}
    }
    fclose(f);

    //////////////////// READ COPY NUMBER VARIANT FROM FILE
    f=fopen(argv[3],"r");
    fgets(line, LINELIMIT, f);
    fgets(line, LINELIMIT, f);

    for (i=0;i<MAXGENES;i++){
        bestIG[i]=-1000000.0;
    }


    while (fgets(line, LINELIMIT, f)>0)
	if (strlen(line)>2*(nn+nu))
	{
            token = strtok(line, DELIM);
            //sscanf(token,"%d",&probeSetID);
	    probeSetID=0;
    		//token = strtok(NULL, DELIM);
    		if (g2id.find(token)!=g2id.end()){

    			id1=g2id[token]-1;

    			present[id1]=1;
    			allPresent[id1]=true;


    			for (index=0;index<nn+nu;index++)
    			{
    			    token = strtok(NULL, DELIM);
                    sscanf(token,"%lf",&v);
                    //cout<<token<<endl;
                    temp[index]=v;
    			}

    			if (variance(temp)>bestIG[id1])
    			{
                        bestIG[id1]=variance(temp);
                        savedProbe[id1]=probeSetID;
                        for (index=0;index<nn+nu;index++)
                        //if (eArr[index][id1][0]<(NULLMARKER+1.0))
                            //eArr[index][id1][0]=v;
                        //else
                            cArr[index][id1][0]=temp[index];
    			}


    		}
    }
    fclose(f);

    f=fopen("Probe2Gene.txt","w");
    for (i=0;i<n;i++)
    if (allPresent[i])
      fprintf(f,"%d\t%s\n",savedProbe[i],id2g[i+1].c_str());
    fclose(f);

    processAllSamples();

    cout<<CHECK<<" ";
    for (i=0;i<nn;i++)
    	cout<<" x"<<i+1<<"x "<<eArr[i][CHECK][0];
    cout<<endl;
    for (i=0;i<nu;i++)
    	cout<<" x"<<i+1<<"x "<<eArr[i+nn][CHECK][0];
    cout<<endl;

    for (i=0;i<n;i++){
    	processGene(i);
    }

    sprintf(filename,"%s/edges.txt",argv[4]);
    f=fopen(filename,"w");
    for (i=0;i<n;i++)
    if (present[i]>0)
    for (j=i+1;j<n;j++)
    if (present[j]>0)
    if (ppi[i][j]>0.00001)
    	fprintf(f,"%d\t%d\t%lf\n",i,j,ppi[i][j]);
    fclose(f);

    // REPLACE SPACES BY _s
    for (i=0;i<n;i++)
    if (id2g[i+1].find_first_of(" ")!=string::npos)
    {
    	id2g[i+1].erase(id2g[i+1].find_first_of(" "));
    	g2id[id2g[i+1]]=i+1;
    }

    int nms=0;
    int ms=0;

    sprintf(filename,"%s/nodes.txt",argv[4]);
    f=fopen(filename,"w");
    for (i=0;i<n;i++)
    	fprintf(f,"%d\t%s\n",i,id2g[i+1].c_str());
    fclose(f);

    sprintf(filename,"%s/population.txt",argv[4]);
    f=fopen(filename,"w");
    for (i=0;i<n;i++)
    	fprintf(f,"%s\n",id2g[i+1].c_str());
    fclose(f);

    sprintf(filename,"%s/dimensions.txt",argv[4]);
    f=fopen(filename,"w");
    for (i=0;i<nu;i++)
    fprintf(f,"%d\tPATIENT%d\n",i,i);
    fclose(f);

    sprintf(filename,"%s/attributes.txt",argv[4]);
    f=fopen(filename,"w");
    for (i=0;i<n;i++)
    {
    	fprintf(f,"%d",i);
    	    int count0=0;
    	    int count1=1;
    	    for (k=0;k<nn;k++)
    	    if (expressed[k][i])
                count1++;
            else
                count0++;

    	for (j=0;j<nu;j++)
    	{


            if ((expressed[j+nn][i])&&(count0>=2*count1))
                fprintf(f,"\t%d",1);
            else if ((!expressed[j+nn][i])&&(count1>=2*count0))
                fprintf(f,"\t%d",1);
            else
                fprintf(f,"\t%d",0);

    	}
        fprintf(f,"\n");
    }
    fclose(f);

    sprintf(filename,"%s/attributes1.txt",argv[4]);
    f=fopen(filename,"w");
    for (i=0;i<n;i++)
    {
    	fprintf(f,"%d",i);
    	for (j=0;j<nn+nu;j++)
    	if (expressed[j][i])
                fprintf(f,"\t%d",1);
        else
                fprintf(f,"\t%d",0);
        fprintf(f,"\n");
    }
    fclose(f);

    sprintf(filename,"%s/expression.txt",argv[4]);
    f=fopen(filename,"w");
    for (i=0;i<n;i++)
    {
    	fprintf(f,"%d",i);
    	for (j=0;j<nn+nu;j++)
    	if (eArr[j][i][0]>(NULLMARKER+1.0))
    		fprintf(f,"\t%.4lf",eArr[j][i][0]);
        else
            fprintf(f,"\t0.0");
    	fprintf(f,"\n");
    }
    fclose(f);

    sprintf(filename,"%s/copynumber.txt",argv[4]);
    f=fopen(filename,"w");
    for (i=0;i<n;i++)
    {
    	fprintf(f,"%d",i);
    	for (j=0;j<nn+nu;j++)
    	if (cArr[j][i][0]>(NULLMARKER+1.0))
    		fprintf(f,"\t%.4lf",cArr[j][i][0]);
        else
            fprintf(f,"\t0.0");
    	fprintf(f,"\n");
    }
    fclose(f);


    sprintf(filename,"%s/expression1.txt",argv[4]);
    f=fopen(filename,"w");
    for (i=0;i<n;i++)
    {
    	fprintf(f,"%s",id2g[i+1].c_str());
    	for (j=0;j<nn+nu;j++)
    	if (eArr[j][i][0]>(NULLMARKER+1.0))
    		fprintf(f,"\t%.4lf",eArr[j][i][0]);
        else
            fprintf(f,"\t0.0");
    	fprintf(f,"\n");
    }
    fclose(f);

    cout<<"n = "<<n<<endl;
    cout<<"nms = "<<nms<<endl;
    cout<<"ms = "<<ms<<endl;

    deallocateGraph(&ppi,n);

    free(line);
    free(tempString1);
    free(tempString2);
    free(present);
    return 0;
}
