// test.cpp : �������̨Ӧ�ó������ڵ㡣
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
struct Delta{           //��ƽ����
 int n;         
 float P,Q,U2;   
 };
struct Uinit{			//��ѹ��ʼֵ
	int n;
	float e,f;
};
struct Jcb{			//СJaccobi����
	float H,N,J,L,R,S;
};


int _tmain(int argc, _TCHAR* argv[])
{
	FILE *fp;
	FILE *fpout;
	ofstream outfile;
	struct Bus sBus[NBUS];
	struct Line sLine[NLINE];
	struct Delta sDelta[NBUS]={0};
	struct Uinit sU[NBUS];
	struct Jcb sJcb[NBUS][NBUS]={0};

	int time=0, i, j, li, lj, k;		//timeΪ����������forѭ��ʹ��ij��Y����ļ������±���li lj; kһ����ʱ����
	float d1, g, b, r, x, m;
	float YG[NLINE][NLINE],YB[NLINE][NLINE];   //YB[NBUS]���


/*********** SRART READING FILE **************/
	// in.txt�ļ��Ƿ���ڡ�
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
	//�ж��ļ��Ƿ�رճɹ�
	if(fclose(fp)==0)	printf("Reading in.txt...OK.\n");  
	else	perror("fclose");	

/******** Make Y Matrix ***********/
	/* ��ʼ��Y����*/		//����ijǶ��ѭ����NumI NumJ�����ص�?
	for(i=0;i<nL;i++)
		for(j=0;j<nL;j++){
			YG[i][j]=0.0;
			YB[i][j]=0.0;
		};
	/* ���㶨�� */
	for(int l=0; l<nL; l++){
		li=sLine[l].NumI-1;
		lj=sLine[l].NumJ-1;
		r=sLine[l].R;
		x=sLine[l].X;
		d1=r*r+x*x;
		g=r/d1;
		b=-x/d1;	
		m=sLine[l].K;
	/* ��֧ͨ· */
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
	/* ��ѹ��֧· */
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

	/* ��� Y ���� */	
	if(NULL==(fp=fopen("result.txt","w"))){printf("Can not open 'result.txt'\n");exit(0);}
	else printf("Writing Y matrix to result.txt...");
	fprintf(fp,"------Y Matrix------\n");
	for(i=0;i<nBus;i++)
		for(j=0;j<nBus;j++)
			if(fabs(YB[i][j]-0.0)>0.000001)
				fprintf(fp,"Y(%3d,%-3d)=(%10.5f,%10.5f)\n",i+1,j+1,YG[i][j],YB[i][j]);

	/*�ж��ļ��Ƿ�رճɹ�*/
	if(fclose(fp)==0)	printf("OK.\n");//system("notepad result.txt"); 
	else	perror("fclose");


/**************** �ڵ��ѹ��ֵ *******************/
	for(i=0; i<nBus; i++){
		sU[i].n=sBus[i].Num;
		sU[i].e=sBus[i].Volt*cos(sBus[i].Phase);
		sU[i].f=sBus[i].Volt*sin(sBus[i].Phase);
	}

/************* ���㲻ƽ����deltaP Q U**************/
	//TODO: д�ɺ�����ʽ
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
			case 2:	break;	//ƽ��ڵ���� ���μӵ���
		}
	}

/************* �γ��Ÿ��Ⱦ��� **************/
	for(i=0; i<nBus; i++){
		if(sBus[i].Type==2) continue;	//û��ƽ��ڵ�
		for(j=0; j<nBus; j++){			
			if(sBus[j].Type==2) continue;  //û��ƽ��ڵ�
			sJcb[i][j].H= YG[i][j]*sU[i].f - YB[i][j]*sU[i].e;
			sJcb[i][j].N= YG[i][j]*sU[i].e + YB[i][j]*sU[i].f;
			if(sBus[i].Type==0){//PQ
				sJcb[i][j].J=-sJcb[i][j].N;sJcb[i][j].L= sJcb[i][j].H;
			}

			if(i==j){
				for(k=0; k<nBus; k++){	//���ʱ��n���ڵ�
					sJcb[i][j].H+=YG[i][k]*sU[k].f + YB[i][k]*sU[k].e;
					sJcb[i][j].N+=YG[i][k]*sU[k].e - YB[i][k]*sU[k].f;
					if(sBus[i].Type==0){//PQ�Ļ���JL
						sJcb[i][j].J+=YG[i][k]*sU[k].e - YB[i][k]*sU[k].f;
						sJcb[i][j].L-=YG[i][k]*sU[k].f + YB[i][k]*sU[k].e;
					}
				}
				if(sBus[i].Type==1)	//PV�Ļ���RS
					sJcb[i][j].R=2*sU[i].f; sJcb[i][j].S= 2*sU[i].e;
			}
		}
	}
	//���ɲ����Jaccobi����
	
		if(NULL==(fp=fopen("result.txt","a"))) {printf("Can not open 'result.txt'\n");exit(0);}
	else printf("Writing Jaccobi matrix to result.txt...");
	fprintf(fp,"\n\n---------------Jaccobi Matrix (%d)------------------\n",time);

	for(i=0;i<nBus;i++){			if(sBus[i].Type==2) continue;	//û��ƽ��ڵ�
		for(j=0, k=0;j<nBus;j++){	if(sBus[j].Type==2) continue;	//û��ƽ��ڵ�
			if(k==0)	fprintf(fp, "%10.3f%10.3f",sJcb[i][j].H,sJcb[i][j].N);
			else		fprintf(fp, "%10.3f%10.3f",sJcb[i][j].J,sJcb[i][j].L);
			if(j==nBus-1 && k==0)	fprintf(fp, "\n"),j=0,k=1;
		}
		fprintf(fp,"\n");
	}
	/*�ж��ļ��Ƿ�رճɹ�*/
	if(fclose(fp)==0)	printf("OK.\n");//system("notepad result.txt"); 
	else	perror("fclose");



/************* ���������� ***************/
		Matrix mJcb(2*nBus-2,2*nBus-2);
		Matrix mDelta(2*nBus-2,1);
		Matrix mU(2*nBus-2,1);
		
		/* ��ֵ��Jaccobi���� */
		struct Jcb *pJ=sJcb[0];
		for(i=0;i<nBus-1;i++){
			if(sBus[i].Type==2) pJ+=nBus;
			for(j=0;j<nBus-1;j++){
				if(sBus[j].Type==2) pJ++;
				mJcb(2*i  ,2*j  )=(*pJ).H;	
				mJcb(2*i  ,2*j+1)=(*pJ).N;
				mJcb(2*i+1,2*j  )=(*pJ).J;
				mJcb(2*i+1,2*j+1)=(*pJ).L;
				pJ++;
			}
		}

		/* ��ֵ����ƽ�������� */
		struct Delta *pD=sDelta;
		for(i=0;i<nBus-1;i++){
			if(sBus[i].Type==2)	pD++;
			mDelta(2*i  ,0)=(*pD).P;
			mDelta(2*i+1,0)=(*pD).Q;
			pD++;
		}
		//cout << "------mDelta-------\n" << mDelta << endl;
		/* �õ���ѹ�������� */
		mU= (!mJcb) * mDelta;
		outfile.open("result.txt",ios::app);
		outfile << "\n------------------- ��" << time << "�ε�ѹ������ -------------------" <<endl;
		outfile << std::right<< setiosflags(ios::fixed) << setprecision(3) << mU << endl;
		outfile.close();

//�����ڵ��ѹ
		struct Uinit *pU=sU;
		for(i=0;i<nBus-1;i++,pU++){
			if(sBus[i].Type==2) pU++;
			(*pU).f+=mU(2*i  ,0);
			(*pU).e+=mU(2*i+1,0);
		}
	/* ��� ������ĵ�ѹ�� */	
	if(NULL==(fp=fopen("result.txt","a"))){printf("Can not open 'result.txt'\n");exit(0);}
	else printf("Writing Correcting U to result.txt...");
	fprintf(fp,"----------��%d�ε�����ĵ�ѹ��-----------\n", time);
	for(i=0;i<nBus;i++)
		fprintf(fp,"�ڵ�%d��e=%10.5f,f=%10.5f\n",i+1, sU[i].e, sU[i].f);
	if(fclose(fp)==0)	printf("OK.\n");//system("notepad result.txt"); 
	else	perror("fclose");
//������·���ʺ�ƽ��ڵ�  PV�ڵ㹦��



	system("pause");
	return 0;
}

