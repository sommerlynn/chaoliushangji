// aaa.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
struct Line     //��·�ṹ��
 {
 int Num,NumI,NumJ;   //��·��  ��ڵ���   �ҽڵ���
 float R,X,B,K;       //����    �翹       ����       ��ȣ�K����1Ϊ��֧ͨ·�� ������1Ϊ��ѹ��֧·�ı�ȣ�   
 };
struct Bus           //�ڵ�ṹ��
 {
 int Num ;          
 float Volt,Phase,GenP,GenQ,LoadP,LoadQ;   
 int Type;
 };
struct ssh           //�ڵ�ṹ��
 {
 int NumI,Num ;          
 float G,B;   
 };
 

#include"stdio.h"
#include"string.h"
#include"math.h"
#include"stdlib.h"

#define NBUS 5			//����Ͷ����˲���ѧ��
#define NLINE 7			

/* Global variables */
int nL,nBus,nSH;
float X[NBUS];
int L;
double def[2*NBUS];

int _tmain(int argc, _TCHAR* argv[])
{

FILE *fp;
FILE *fpout;

int i,j,k,l,h,n,v;
int i1,i2,i3,kp,kq;

float d1,d2,d3,d4,d5,d6,r,x,g,b,tt,LL,e,ps,qs,shsh,m;
struct Line sLine[NLINE];
struct Bus sBus[NBUS];
struct ssh sSH[NBUS];
float YG[NBUS][NBUS],YB[NBUS][NBUS];


i1=i2=i3=0;
d1=d2=d3=d4=d5=d6=ps=qs=0.0;

/* in.txt�ļ��Ƿ���� */
for(i=0;i<NBUS;i++)					//���������ɾ����
if((fp=fopen("in.txt","r"))==NULL)		//�����ڣ��ж���������
  { printf("Can not open the file named 'in.txt' \n");
    exit(0);
  };
//�������ȡ��һ�� �ڵ����� ��·���� SH����(��û��SH)
	fscanf(fp,"%d,%d,%d",&nBus,&nL,&nSH);

/*read bus data to sBus*/
for(i=0;i<nBus;i++){
	/* ��ʼ��sBus�ṹ�� */
 sBus[i].Num=sBus[i].Type=0;sBus[i].Volt=1.0;
 sBus[i].Phase=sBus[i].GenP=sBus[i].GenQ=sBus[i].LoadP=sBus[i].LoadQ=0.0;
	/* ��in.txt�ж�ȡNum Volt Phase GP GQ LP LQ Type */
 fscanf(fp,"%d,%f,%f,%f,%f,%f,%f,%d",sBus[i].Num,sBus[i].Volt,sBus[i].Phase,sBus[i].GenP,sBus[i].GenQ,sBus[i].LoadP,sBus[i].LoadQ,sBus[i].Type);
 //fscanf(fp,"%d,%f,%f,%f,%f,%f,%f,%d",&i1,&d1,&d2,&d3,&d4,&d5,&d6,&i2);
 //sBus[i].Num=i1; sBus[i].Volt=d1; sBus[i].Phase=d2; sBus[i].GenP=d3; sBus[i].GenQ=d4;sBus[i].LoadP=d5,sBus[i].LoadQ=d6;sBus[i].Type=i2;
 };

/*read line data to sLine*/
for(i=0;i<nL;i++){
 sLine[i].Num=sLine[i].NumI=sLine[i].NumJ=0;
 sLine[i].R=sLine[i].X=sLine[i].B=0.0;sLine[i].K=1.0;
 fscanf(fp,"%2d %3d %3d %f %f %f %f",sLine[i].Num,sLine[i].NumI,sLine[i].NumJ,sLine[i].R,sLine[i].X,sLine[i].B,sLine[i].K);
 //fscanf(fp,"%2d %3d %3d %f %f %f %f",&i1,&i2,&i3,&d1,&d2,&d3,&d4);
 //sLine[i].Num=i1;sLine[i].NumI=i2;sLine[i].NumJ=i3;sLine[i].R=d1;sLine[i].X=d2;sLine[i].B=d3;sLine[i].K=d4;
 }

/* read sh data ��û�����ݣ�*/
for(i=0;i<nSH;i++){
 sSH[i].Num=sSH[i].NumI=0;
 sSH[i].G=sSH[i].B=0.0;
 fscanf(fp,"%2d %3d %f",&i1,&i2,&d1);
 sSH[i].Num=i1;sSH[i].NumI=i2;sSH[i].B=d1;
 }
if(fp!=NULL) fclose(fp); 
/* end reading file */


/*Make Y Matrix*/
for(i=1;i<nBus+1;i++)
	for(j=1;j<nBus+1;j++){
		YG[i][j]=0.0;
		YB[i][j]=0.0;
	};

for(l=0; l<nL; l++){
	i=sLine[l].NumI;
	j=sLine[l].NumJ;
	r=sLine[l].R;
	x=sLine[l].X;
	d1=r*r+x*x;
	g=r/d1;			//
	b=-x/d1;		//
	m=sLine[l].K;	//���
 
	if(fabs(sLine[l].K-1.0)<0.000001){     //��֧ͨ·
		YG[i][i]= YG[i][i] +g;
		YG[j][j]= YG[j][j] +g;
		YB[i][i]= YB[i][i] +b +sLine[l].B;
		YB[j][j]= YB[j][j] +b+ sLine[l].B;
		YG[i][j]= YG[i][j] -g;
		YG[j][i]= YG[j][i] -g;
		YB[i][j]= YB[i][j] -b;
		YB[j][i]= YB[j][i] -b;
	}

   else                            //��ѹ��֧·
   {
		YG[i][i]= YG[i][i] + g/m+g*(m-1)/m;
		YG[j][j]= YG[j][j] + g/m+g*(1-m)/(m*m);
		YB[i][i]= YB[i][i] + b/m+b*(m-1)/m		+ sLine[l].B;
		YB[j][j]= YB[j][j] + b/m+b*(1-m)/(m*m)	+ sLine[l].B;
		YG[i][j]= YG[i][j] - g/m;
		YG[j][i]= YG[j][i] - g/m;
		YB[i][j]= YB[i][j] - b/m;
		YB[j][i]= YB[j][i] - b/m;



   }
 }

/* Check the Y matrix */
if((fp=fopen("GGBB.txt","w"))==NULL){
  printf("Can not open the file named 'GGBB.txt' \n");exit(0);}
fprintf(fp,"---Y Matrix---\n");
for(i=1;i<nBus+1;i++)for(j=1;j<nBus+1;j++)if(fabs(YB[i][j]-0.0)>0.000001)
fprintf(fp,"Y(%3d,%-3d)=(%10.5f,%10.5f)\n",i,j,YG[i][j],YB[i][j]);

//for(i=1;i<nBus+1;i++)
//	for(j=1;j<nBus+1;j++)
//		if(fabs(YB[i][j]-0.0)>0.000001){
//			fprintf(fp,

if(fp!=NULL) fclose(fp);


//�ڵ��ѹ����ֵ
float Volt[NBUS][2];	//Volt[i][0]=angle or ei; Volt[i][1]=amp or fi;

double Vei[NBUS],Vfi[NBUS];

//��Volt��phase��ֵ���ڵ��ѹ����Volt;
for (i=0; i<NBUS; i++)
	if (sBus[i].Type == 2){			//��ƽ��ڵ�
		Volt[i][0] = sBus[i].Phase;
		Volt[i][1] = sBus[i].Volt;
	}
	else{							//����ƽ��ڵ�
		Volt[i][0] = 0.0;
		Volt[i][1] = 1.0;
	};


//���㹦�ʲ�ƽ����



//�γ��ſ˱Ⱦ���



//����������



//�����ڵ��ѹ


//������·���ʺ�ƽ��ڵ�  PV�ڵ㹦��


	return 0;
}

