// test.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include"stdio.h"
#include"string.h"
#include"math.h"
#include"stdlib.h"
#include "matrix.h"

/************** Matrix calculate Need ***************/
#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
#define STD std
#else
#define STD
#endif

#ifndef _NO_TEMPLATE
typedef matrix<double> Matrix;
#else
typedef matrix Matrix;
#endif

#ifndef _NO_EXCEPTION
#  define TRYBEGIN()	try {
#  define CATCHERROR()	} catch (const STD::exception& e) { \
						cerr << "Error: " << e.what() << endl; }
#else
#  define TRYBEGIN()
#  define CATCHERROR()
#endif


/* Global variables */
#define NBUS 5			
#define NLINE 7		

int nL,nBus,nSH;	
float X[NBUS];		//???
int L;				//???
double def[2*NBUS];	//???

struct Line     //线路结构体
 {
 int Num,NumI,NumJ;   //线路号  左节点名   右节点名
 float R,X,B,K;       //电阻    电抗       电纳       变比（K等于1为普通支路， 不等于1为变压器支路的变比）   
 };
struct Bus           //节点结构体
 {
 int Num ;          
 float Volt,Phase,GenP,GenQ,LoadP,LoadQ;   
 int Type;
 };
struct Delta{           //不平衡量
 int n;         
 float P,Q,U2;   
 };
struct Uinit{			//电压初始值
	int n;
	float e,f;
};
struct Jcb{			//小Jaccobi矩阵
	float H[NBUS][NBUS],N[NBUS][NBUS],J[NBUS][NBUS],L[NBUS][NBUS],R[NBUS][NBUS],S[NBUS][NBUS];
};


int _tmain(int argc, _TCHAR* argv[])
{
	FILE *fp;
	FILE *fpout;
	struct Bus sBus[NBUS];
	struct Line sLine[NLINE];
	struct Delta sDelta[NBUS]={0};
	struct Uinit sU[NBUS];
	struct Jcb sJcb={0};

	int time=0, i, j, li, lj, k;		//time为迭代次数，for循环使用ij，Y矩阵的计算中下标用li lj; k一个临时变量
	float d1, g, b, r, x, m;
	float YG[NLINE][NLINE],YB[NLINE][NLINE];   //YB[NBUS]溢出


/*********** SRART READING FILE **************/
	// in.txt文件是否存在、
	if((fp=fopen("in.txt","r"))==NULL){
		printf("Can not open the file named 'in.txt' \n");
		exit(0);
	};
	fscanf(fp,"%d,%d,%d",&nBus,&nL,&nSH);	//5,7,0
	//read bus data to sBus
	for(i=0;i<nBus;i++){
		sBus[i].Num=sBus[i].Type=0;sBus[i].Volt=1.0;
		sBus[i].Phase=sBus[i].GenP=sBus[i].GenQ=sBus[i].LoadP=sBus[i].LoadQ=0.0;
		fscanf(fp,"%d,%f,%f,%f,%f,%f,%f,%d",&sBus[i].Num,&sBus[i].Volt,&sBus[i].Phase,&sBus[i].GenP,&sBus[i].GenQ,&sBus[i].LoadP,&sBus[i].LoadQ,&sBus[i].Type);
	};
	//read line data to sLine
	for(i=0;i<nL;i++){
		sLine[i].Num=sLine[i].NumI=sLine[i].NumJ=0;
		sLine[i].R=sLine[i].X=sLine[i].B=0.0;sLine[i].K=1.0;
		fscanf(fp,"%2d %3d %3d %f %f %f %f",&sLine[i].Num,&sLine[i].NumI,&sLine[i].NumJ,&sLine[i].R,&sLine[i].X,&sLine[i].B,&sLine[i].K);
	}
	//判断文件是否关闭成功
	if(fclose(fp)==0)	printf("Reading in.txt...OK.\n");  
	else	perror("fclose");	

/******** Make Y Matrix ***********/
	/* 初始化Y数组*/		//这里ij嵌套循环与NumI NumJ可能重叠?
	for(i=0;i<nL;i++)
		for(j=0;j<nL;j++){
			YG[i][j]=0.0;
			YB[i][j]=0.0;
		};
	/* 运算定义 */
	for(int l=0; l<nL; l++){
		li=sLine[l].NumI-1;
		lj=sLine[l].NumJ-1;
		r=sLine[l].R;
		x=sLine[l].X;
		d1=r*r+x*x;
		g=r/d1;
		b=-x/d1;	
		m=sLine[l].K;
	/* 普通支路 */
		if(fabs(sLine[l].K-1.0)<0.000001){     
			YG[li][li]= YG[li][li] + g;
			YG[lj][lj]= YG[lj][lj] + g;
			YB[li][li]= YB[li][li] + b + sLine[l].B;
			YB[lj][lj]= YB[lj][lj] + b + sLine[l].B;
			YG[li][lj]= YG[li][lj] - g;
			YG[lj][li]= YG[lj][li] - g;
			YB[li][lj]= YB[li][lj] - b;
			YB[lj][li]= YB[lj][li] - b;
		}
	/* 变压器支路 */
		else{
			YG[li][li]= YG[li][li] + g/m+g*(m-1)/m;
			YG[lj][lj]= YG[lj][lj] + g/m+g*(1-m)/(m*m);
			YB[li][li]= YB[li][li] + b/m+b*(m-1)/m	+ sLine[l].B;
			YB[lj][lj]= YB[lj][lj] + b/m+b*(1-m)/(m*m) + sLine[l].B;
			YG[li][lj]= YG[li][lj] - g/m;
			YG[lj][li]= YG[lj][li] - g/m;
			YB[li][lj]= YB[li][lj] - b/m;
			YB[lj][li]= YB[lj][li] - b/m;
		}
	}

	/* 输出 Y 矩阵 */	
	if(NULL==(fp=fopen("result.txt","w"))){printf("Can not open 'result.txt'\n");exit(0);}
	else printf("Writing Y matrix to result.txt...");
	fprintf(fp,"------Y Matrix------\n");
	for(i=0;i<nBus;i++)
		for(j=0;j<nBus;j++)
			if(fabs(YB[i][j]-0.0)>0.000001)
				fprintf(fp,"Y(%3d,%-3d)=(%10.5f,%10.5f)\n",i+1,j+1,YG[i][j],YB[i][j]);

	/*判断文件是否关闭成功*/
	if(fclose(fp)==0)	printf("OK.\n");//system("notepad result.txt"); 
	else	perror("fclose");


/**************** 节点电压初值 *******************/
	for(i=0; i<nBus; i++){
		sU[i].n=sBus[i].Num;
		sU[i].e=sBus[i].Volt*cos(sBus[i].Phase);
		sU[i].f=sBus[i].Volt*sin(sBus[i].Phase);
	}

/************* 计算不平衡量deltaP Q U**************/
	//TODO: 写成函数形式
	for(i=0; i<nBus; i++){
		switch(sBus[i].Type){
			case 0:	//PQ
				sDelta[i].n=sBus[i].Num;
				sDelta[i].P=sBus[i].GenP-sBus[i].LoadP;
				sDelta[i].Q=sBus[i].GenQ-sBus[i].LoadQ;
				for(j=0; j<nBus; j++){
					sDelta[i].P-=(sU[i].e*(YG[i][j]*sU[j].e-YB[i][j]*sU[j].f)-sU[i].f*(YG[i][j]*sU[j].f+YB[i][j]*sU[j].e));
					sDelta[i].Q-=(sU[i].f*(YG[i][j]*sU[j].e-YB[i][j]*sU[j].f)-sU[i].e*(YG[i][j]*sU[j].f+YB[i][j]*sU[j].e));
				}break;
			case 1:	//PV
				sDelta[i].n=sBus[i].Num;
				sDelta[i].P=sBus[i].GenP-sBus[i].LoadP;
				sDelta[i].U2=sBus[i].Volt*sBus[i].Volt-(sU[i].e*sU[i].e+sU[i].f*sU[i].f);
				for(j=0; j<nBus; j++)
					sDelta[i].P-=(sU[i].e*(YG[i][j]*sU[j].e-YB[i][j]*sU[j].f)-sU[i].f*(YG[i][j]*sU[j].f+YB[i][j]*sU[j].e));
				break;
			case 2:	break;	//平衡节点给定 不参加迭代
		}
	}

/************* 形成雅各比矩阵 **************/
	for(i=0; i<nBus; i++){
		if(sBus[i].Type==2) continue;	//没有平衡节点
		for(j=0; j<nBus; j++){			
			if(sBus[j].Type==2) continue;  //没有平衡节点
			sJcb.H[i][j]= YG[i][j]*sU[i].f - YB[i][j]*sU[i].e;
			sJcb.N[i][j]= YG[i][j]*sU[i].e + YB[i][j]*sU[i].f;
			if(sBus[i].Type==0){//PQ
				sJcb.J[i][j]=-sJcb.N[i][j];sJcb.L[i][j]= sJcb.H[i][j];
			}

			if(i==j){
				for(k=0; k<nBus; k++){	//求和时是n个节点
					sJcb.H[i][j]+=YG[i][k]*sU[k].f + YB[i][k]*sU[k].e;
					sJcb.N[i][j]+=YG[i][k]*sU[k].e - YB[i][k]*sU[k].f;
					if(sBus[i].Type==0){//PQ的话算JL
						sJcb.J[i][j]+=YG[i][k]*sU[k].e - YB[i][k]*sU[k].f;
						sJcb.L[i][j]-=YG[i][k]*sU[k].f + YB[i][k]*sU[k].e;
					}
				}
				if(sBus[i].Type==1)	//PV的话算RS
					sJcb.R[i][j]=2*sU[i].f; sJcb.S[i][j]= 2*sU[i].e;
			}
		}
	}
	//生成并输出Jaccobi矩阵
	fp=fopen("result.txt","a");
		if(NULL==fp){printf("Can not open 'result.txt'\n");exit(0);}
	else printf("Writing Jaccobi matrix to result.txt...");
	fprintf(fp,"\n\n---------------Jaccobi Matrix (%d)------------------\n",time);

	for(i=0;i<nBus;i++){			if(sBus[i].Type==2) continue;	//没有平衡节点
		for(j=0, k=0;j<nBus;j++){	if(sBus[j].Type==2) continue;	//没有平衡节点
			if(k==0)	fprintf(fp, "%10.3f%10.3f",sJcb.H[i][j],sJcb.N[i][j]);
			else		fprintf(fp, "%10.3f%10.3f",sJcb.J[i][j],sJcb.L[i][j]);
			if(j==nBus-1 && k==0)	fprintf(fp, "\n"),j=0,k=1;
		}
		fprintf(fp,"\n");
	}
		printf("OK.\n");



/************* 解修正方程 ***************/
		//Matrix mJaccobi(2*nBus,2*nBus);
		//Matrix mDelta(2*nBus);
		//Matrix mU(2*nBus);
		struct Jcb *p=&sJcb;
		for(i=0;i<nBus;i++){
		for(j=0, k=0;j<nBus;j++){
			printf("%10.3f",(*p).H[i][j]);

		}
					printf("\n");
	}

//修正节点电压


//计算线路功率和平衡节点  PV节点功率



	system("pause");
	return 0;
}

