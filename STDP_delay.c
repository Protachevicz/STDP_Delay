//This code reproduces the results presents in the article : 
//Plastic neural network with transmission delays promotes 
//equivalence between function and structure. 
//Chaos, Solitons $\&$ Fractals 171,  2023, 113480.

#include <unistd.h> 
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
//******  Randon numbers Parameters ********************//
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI acos(-1.0)
//******  Networks Parameters ********************//
#define rede 100      //Number of neurons in the subnetworks
#define subredes 4   //Number of Subnetworks 
#define numerocnl  1.0*rede*(rede-1)  //Non-local connection number for each subnetwork 
#define conexMAX rede*subredes   // probable maximum number of each neuron connections
#define p_exc 1.0    // % excitatory neurons (80%)
#define links 0.05*p_exc*rede*rede    // Maximal number of connection between the subnetworks 
//******   Parameters ********************//
#define NN rede*subredes   // number of nodes in the network
#define N 4      // numero of equations by neuron
#define n 1.0001*pow(10,7)  //total number of steps   

#define transient 1000 //in ms
#define NMD 400000 // maximal number of spikes
#define TRANSIENTE_KMEDIO 100 
//*********delay**********//
#define delay_in 0.01   //inicial internal delay - inicial, miliseconds // 0.01 ms - 10.0 ms,  
#define delay_fin 0.01  // final delay interno 
#define delay_inter_in 0.01  // delay externo inicial // 0.01 ms - 20.0 ms
#define delay_inter_fin 0.01   //final external delay
#define passo_delay 0.04   
#define passo_delay_inter 0.3  
#define repeticoes 1    // number of initial condictions 
#define delayMax 3200    // maximal delay values
#define delayDP 0.00    //standard deviation of delays
//******acoplamentos******//
#define razao 1       // relacion between the intensity of inhibitory and excitatory coupling
#define conex_inicial_exc 0.001 // initial intensity of connections
#define conex_inicial_inib razao*conex_inicial_exc
#define desvio_conex_inicial 0.000
#define limite_sup_exc 0.01  // limite superior dos pesos das conexões exc  // maximo 0.002 para 4 redes de 100
#define limite_sup_inib razao*limite_sup_exc  // limite superior das conexões inib
/**************parametros das chaves***************/
#define Ngif 51      // number of layers gif pra montar animação
#define passoM 0.5*pow(10,6)  // time step to save matrixes, aprox = n/Ngif
#define passoAD 1*pow(10,4)  // time step to save mean coupling
//******  STDP Parameters ********************//
#define ligar 0    //  n = turn off stdp,  0 = turn on stdp
#define A1 1.0
#define A2 0.5
#define tau1 1.8
#define tau2 6.0 
#define delta_exc 0.00001
//****** iSTDP Parameters ********************//
#define g 1.0
#define	alpha 1.0
#define	beta 10.0
#define	alpha_plus 0.94
#define	alpha_minus 1.1
#define delta_inib 0.00001
//*********************************************
float ran1(long *idum);
float gasdev(long *idum);
#define NR_END 1
#define FREE_ARG char*
void nrerror(char error_text[]);
int *vector(long nl,long nh);
int **imatrix(long nrl, long nrh, long ncl, long nch);
void free_vector(int *v, long nl, long nh);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

double *dvector(long nl,long nh);
void free_dvector(double *v, long nl, long nh);

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void derivs(double y[],double df[], double *Gexc, double *Gini,double *I);

FILE *o;
FILE *p;
FILE *q;
FILE *q2;
FILE *q3;
FILE *q4;
FILE *fases;
FILE *coup;
FILE *raster1,*raster2,*raster3,*raster4;

int main()
{
  int *conextotal,**listaADJ, **ADJ,chave2,contdeltat,contdeltatE,contdeltatI, Vr[NN+2],**delay,t2,auxx, contaux[subredes+1][subredes+2];
  double **listaconex,*df,*y,*I,deltatmedio,deltatmedioE,deltatmedioI, chave3,**SynMenexc,**SynMenini, chave4, chave5,aux1, acop_medioex[subredes+1][subredes+2];
  
  int i,j,ii,jj,aux,aux2,auxx2,t,cont,cont2,chave,contM,g_cont,rep,m;
  double x[N*NN+2],xant[N*NN+2],tempo,h,a[N*NN+2],b[N*NN+2],c[N*NN+2],pico[N*NN+2],Deltat,W1,W2,temp,nkk[NN+2],tudant[NN+2],tud[NN+2],W1_inib,W2_inib,gnorm,sum;
	
  float ***Matriz,***MatrizI, Imedio;

  double aux_I,kinicial,kfinal,R[subredes+1],real,compl,Rmedio[subredes+1],**nk,*phi, tpeak[NN+2][NN+2],tauex, delay_medio, delay_medio_inter,aux3, gex_medio, ginib_medio,kmedio,*Gexc,*Gini,GexcM,GiniM;
  int contR[subredes+1],*k,*kmax, nfibras, auxT, imprimir, contT[NN+2];  
  double GEX, GEX2, GINIB, GINIB2, K, K2, Rrep, Rrep2, Tmedio[NN+2],T;
  long idum; 
  
  conextotal=vector(1,NN+1);
  listaADJ=imatrix(1,NN+2,1,conexMAX+2);
  ADJ=imatrix(1,NN+2,1,conexMAX+2);
  listaconex=dmatrix(1,conexMAX+1,1,conexMAX+1);
  I=dvector(1,NN+1);
  y=dvector(1,N*NN+1);
  df=dvector(1,N*NN+1);
  Matriz=f3tensor(1,NN+2,1,NN+2,1,Ngif+2);
  MatrizI=f3tensor(1,NN+2,1,NN+2,1,Ngif+2);
  Gexc=dvector(1,NN+1);
  Gini=dvector(1,NN+1);
  k=vector(1,NN+1);
  kmax=vector(1,NN+1);
  phi=dvector(1,NN+1);
  nk=dmatrix(1,NMD+2,1,NN+2);
  delay=imatrix(1,NN+2,1,NN+2);
  SynMenexc=dmatrix(1,NN+2,1,delayMax+2);
  SynMenini=dmatrix(1,NN+2,1,delayMax+2);

  o=fopen("tempo_gex_gin_k_deltat_gex=max0.01_gi=ge_4NN100_pinter0.05_din0.0_dex0.0.dat","wt"); 
  p=fopen("Matriz_layers_gex=max0.01_gi=ge_4NN100_pinter0.05_din0.0_dex0.0.dat","wt"); 
    q=fopen("xindex_tempo_gex=max0.01_gi=ge_4NN100_pinter0.05_din0.0_dex0.0.dat","wt"); 
  q2=fopen("delayinter_delayinter_mMax_RmMax_Rms_Tmedio_gex=max0.01_gi=ge_4NN100_pinter0.05_din0.0_dex0.0.dat","wt"); 
  q3=fopen("POKxtempo_gex=max0.01_gi=ge_4NN100_pinter0.05_din0.0_dex0.0.dat","wt"); 
   q4=fopen("tempo_Gexc_Ginb_Iexc_Iinib_IexcM_IinibM_ItotalM_gex=max0.01_gi=ge_4NN100_pinter0.05_din0.0_dex0.0.dat","wt");       
    fases=fopen("tempo_neuron_phase_gex=max0.01_gi=ge_4NN100_pinter0.05_din0.0_dex0.0.dat","wt");       
   coup=fopen("post_pre_acop_acopdiff_gex=max0.01_gi=ge_4NN100_pinter0.05_din0.0_dex0.0.dat","wt");  
   raster1=fopen("Raster1.dat","wt");
   raster2=fopen("Raster2.dat","wt");
   raster3=fopen("Raster3.dat","wt");
   raster4=fopen("Raster4.dat","wt");

  if (o==NULL)
    {
      puts ("erro no arquivo");
      getchar();
      exit(1);
    }

  if (p==NULL)
    {
      puts ("erro no arquivo");
      getchar();
      exit(1);
    }

  if (q==NULL)
    {
      puts ("erro no arquivo");
      getchar();
      exit(1);
    }

idum=-123456789;

for(delay_medio=delay_in;delay_medio<=delay_fin;delay_medio=delay_medio+passo_delay)
for(delay_medio_inter=delay_inter_in;delay_medio_inter<=delay_inter_fin;delay_medio_inter=delay_medio_inter+passo_delay_inter)
{
    GEX= 0.0;
    GEX2=0.0;
    GINIB =0.0;
    GINIB2 = 0.0;
    K=0.0;
    K2= 0.0;
    Rrep =0.0;
    Rrep2= 0.0;
    auxT=0.0; 	  
    		
for(rep=1;rep<=repeticoes;rep=rep+1)     
{
  h=0.01; // integration step
  chave=1.0;
  chave2=1.0;
  gnorm=pow(beta,beta)*exp(-beta);  // iSTDP
  contM=1.0;
  tauex=2.728;    
    
  imprimir = 0;
      
        if(delay_medio==delay_fin)
		if(delay_medio_inter==delay_inter_fin)
		if(rep==repeticoes) 
			imprimir = 1;  
  
  for(i=1;i<=rede*subredes;i++)
    {
      conextotal[i]=0.0;         
      tud[i]=transient;
      tudant[i]=transient;
      Vr[i]=20.0;
      Tmedio[i]=0.0;
      contT[i]=0;
                  
      for(j=1;j<=rede*subredes;j++)
	     { 
		  ADJ[i][j]=1.0; 	              
	        if(i==j)
	         ADJ[i][j]=0.0; 	      
	     }
    }

for(aux=0;aux<=rede*subredes-rede;aux=aux+rede)   
  {	
    for(aux2=aux+1;aux2<=aux+rede;aux2=aux2+1)
      if(aux2>=p_exc*rede+aux+1) 
	Vr[aux2]=-75.0; 	
    
    sum=0.0;
    if(numerocnl!=1.0*rede*(rede-1))
      while(sum<rede*(rede-1)-numerocnl)  
	{
	  i=(int)rede*ran1(& idum)+aux+1;    // póstsynaptic
	  j=(int)rede*ran1(& idum)+aux+1;    // présynaptic   
	  
	  if(ADJ[i][j]>0)   
	    {
	      ADJ[i][j] = 0;   	 
	      sum=sum+1;    
	    }  
	}
  } 
/////////////////////////////////////////////
 if(subredes>1)       
   for(i=1;i<=subredes;i++) 
     for(j=1;j<=subredes;j++)
       if(i!=j)
	 for(aux=rede*(i-1)+1;aux<=rede*(i-1)+rede;aux=aux+1)
	   for(aux2=rede*(j-1)+1;aux2<=rede*(j-1)+rede;aux2=aux2+1)
	     {
	       ADJ[aux][aux2]=0.0;
	     }
 
 if(subredes>1)
   for(i=1;i<=subredes;i++) 
     for(j=1;j<=subredes;j++)
       if(i!=j)
	 {
	   nfibras= links;                   
	   
	   if(nfibras!=0)
	     { 
	       sum=0.0;
	       
	       while(sum<nfibras)
		 {
		   aux=(int)rede*ran1(& idum)+rede*(i-1)+1;           
		   aux2=(int)rede*p_exc*ran1(& idum)+rede*(j-1)+1;   
		   
		   if(ADJ[aux][aux2]<1.0)
		     {
		       ADJ[aux][aux2]=1.0;     // ADJ[póst][pré]
		       sum=sum+1;
		     }                          //printf("%d %d %d %d\n",i,j,sum,0);					    	
		 }
	     }
	 }
 
 for(i=1;i<=subredes;i++) //  ADJ[post][pré]
   for(j=1;j<=subredes;j++)
     for(aux=rede*(i-1)+1;aux<=rede*(i-1)+rede;aux=aux+1)   //póst
       for(aux2=rede*(j-1)+1;aux2<=rede*(j-1)+rede;aux2=aux2+1)  //pré
	 if(ADJ[aux][aux2]==1) 
	   {			
	     if(i==j)
	       {	
		 conextotal[aux2]=conextotal[aux2]+1;  
		 listaADJ[aux2][conextotal[aux2]]=aux;       
	         
		 listaconex[aux][aux2]=desvio_conex_inicial*gasdev(& idum)+conex_inicial_exc;  
	         
		 auxx=((delayDP)*gasdev(& idum) + delay_medio)*(1/h);	 
		 while(auxx<0 || h*auxx>2.0*delay_medio)             
		   auxx=((delayDP)*gasdev(& idum) + delay_medio)*(1/h);
		 
		 delay[aux][aux2]=auxx;	  	      
	       } 
	     else
	       {	 
		 conextotal[aux2]=conextotal[aux2]+1;  
		 listaADJ[aux2][conextotal[aux2]]=aux;       
	         
		 listaconex[aux][aux2]=desvio_conex_inicial*gasdev(& idum)+conex_inicial_exc;  
	         
		 auxx=((delayDP)*gasdev(& idum) + delay_medio_inter)*(1/h);	 
		 while(auxx<0 || h*auxx>2.0*delay_medio_inter)   
		   auxx=((delayDP)*gasdev(& idum) + delay_medio_inter)*(1/h);
		 
		 delay[aux][aux2]=auxx;	  	      
	       } 			
	   } 
 //****************  Initial Conditions **********************************************
 for(i=1;i<=NN;i++)     
   {   
     x[1+(i-1)*N]=-60.0+80.0*ran1(& idum);   
     x[2+(i-1)*N]=ran1(& idum);
     x[3+(i-1)*N]=ran1(& idum);
     x[4+(i-1)*N]=ran1(& idum);   
     
     xant[i]=x[1+(i-1)*N];
     I[i]=10.0+1.0*ran1(& idum); 
         
     for(ii=1;ii<=80;ii++)
     for(jj=1;jj<=80;jj++)
     {
		 if(I[ii]>I[jj])
		 {
		  aux_I= I[ii];
		  I[ii]=I[jj];
		  I[jj]=aux_I;
		 }	 
	 }
	 
     for(ii=101;ii<=180;ii++)
     for(jj=101;jj<=180;jj++)
     {
		 if(I[ii]>I[jj])
		 {
		  aux_I= I[ii];
		  I[ii]=I[jj];
		  I[jj]=aux_I;
		 }	 
	 } 
	 
for(ii=201;ii<=280;ii++)
     for(jj=201;jj<=280;jj++)
     {
		 if(I[ii]>I[jj])
		 {
		  aux_I= I[ii];
		  I[ii]=I[jj];
		  I[jj]=aux_I;
		 }	 
	 } 
	 
     for(ii=301;ii<=380;ii++)
     for(jj=301;jj<=380;jj++)
     {
		 if(I[ii]>I[jj])
		 {
		  aux_I= I[ii];
		  I[ii]=I[jj];
		  I[jj]=aux_I;
		 }	 
	 } 
	 
    for(ii=81;ii<=100;ii++)
     for(jj=81;jj<=100;jj++)
     {
		 if(I[ii]>I[jj])
		 {
		  aux_I= I[ii];
		  I[ii]=I[jj];
		  I[jj]=aux_I;
		 }	 
	 }
	  
     for(ii=181;ii<=200;ii++)
     for(jj=181;jj<=200;jj++)
     {
		 if(I[ii]>I[jj])
		 {
		  aux_I= I[ii];
		  I[ii]=I[jj];
		  I[jj]=aux_I;
		 }	 
	 }
     
          for(ii=281;ii<=300;ii++)
     for(jj=281;jj<=300;jj++)
     {
		 if(I[ii]>I[jj])
		 {
		  aux_I= I[ii];
		  I[ii]=I[jj];
		  I[jj]=aux_I;
		 }	 
	 }
	 
	 
     for(ii=381;ii<=400;ii++)
     for(jj=381;jj<=400;jj++)
     {
		 if(I[ii]>I[jj])
		 {
		  aux_I= I[ii];
		  I[ii]=I[jj];
		  I[jj]=aux_I;
		 }	 
	 }  
          
     pico[i]=x[1+(i-1)*N];   
     
     k[i]=0.0;            
     kmax[i]=0.0;	  
     nk[1][i]=0.0;        
     Gexc[i]=0.0;
     Gini[i]=0.0;
     for(j=1;j<=delayMax;j++) 
       {
	 SynMenexc[i][j]=0.0;
	 SynMenini[i][j]=0.0;
       }
   }	    

 for(aux=0;aux<=rede*subredes-rede;aux=aux+rede)   
   {
     cont=0;  
     for(i=1+aux;i<=p_exc*rede+aux;i=i+1) 
       {
	 cont=cont+1;
	 nkk[cont]=I[i];
       }
     
     for(i=1;i<=p_exc*rede;i=i+1) 
       for(j=1;j<=p_exc*rede-1;j=j+1) 
	 if(nkk[j+1]<=nkk[j])   

	   {
             temp=nkk[j];
             nkk[j]=nkk[j+1];
             nkk[j+1]=temp;
	   }
     
     cont=0.0;
     for(i=1+aux;i<=p_exc*rede+aux;i=i+1)  
       {
	 cont=cont+1;
	 I[i]=nkk[cont];           
       }   
          if(p_exc<1)
       {
	 cont=0;  
	 for(i=p_exc*rede+1+aux;i<=rede+aux;i=i+1) 
	   {
	     cont=cont+1;
	     nkk[cont]=I[i];
	   }
	 
	 for(i=1;i<=rede-p_exc*rede;i=i+1) 
	   for(j=1;j<=rede-p_exc*rede;j=j+1) 
	     if(nkk[j+1]<=nkk[j])   
	    	       {
		 temp=nkk[j];
		 nkk[j]=nkk[j+1];
		 nkk[j+1]=temp;
	       }
	 
	 cont=0.0;
	 for(i=p_exc*rede+1+aux;i<=rede+aux;i=i+1)  
	   {
	     cont=cont+1;
	     I[i]=nkk[cont];           
	
	   }     
       }   
     
   }
 tempo=0.0;
 deltatmedio=0.0;
 contdeltat=0;
 deltatmedioE=0.0;
 contdeltatE=0;
 deltatmedioI=0.0;
 contdeltatI=0;
 gex_medio=0.0;
 g_cont=0.0;
 ginib_medio=0.0;
 kmedio=0.0;
 t2=0;
  
 for(i=1;i<=rede*subredes;i++)
   for(j=1;j<=conextotal[i];j++)
     tpeak[i][j]=0.0;
  
 for(t=1;t<=n;t++)  
   {                   
     tempo=tempo+h;
     t2=t2+1; 
     
     aux=1.0;        
     for(i=1;i<=NN;i++)  
       { 							
	 Gexc[i]=Gexc[i]*exp(-h/tauex)+SynMenexc[i][t2];  
	 Gini[i]=Gini[i]*exp(-h/tauex)+SynMenini[i][t2];  		
	 
	 SynMenexc[i][t2]=0.0;
	 SynMenini[i][t2]=0.0;        
	 
	 y[1+(i-1)*N]=x[1+(i-1)*N];	
	 y[2+(i-1)*N]=x[2+(i-1)*N];		
	 y[3+(i-1)*N]=x[3+(i-1)*N];		
	 y[4+(i-1)*N]=x[4+(i-1)*N];	
	 
	 if(x[1+(i-1)*N]<-30.0)
	   pico[i]=0.0;
	 
	 if(x[1+(i-1)*N]>pico[i]) 
	   {  
	     //---------------------Saving the spike times--------------------------------//
	     tudant[i]=tud[i];  
	     
	     tud[i]=tempo-h-(xant[i]*h)/(x[1+(i-1)*N]-xant[i]); 	  
	 
	  if(tempo<10000 || tempo> 90000)
	  {
	  if(i>0 && i<101) fprintf(raster1, "%f %d\n", tempo, i);
	   if(i>100 && i<201) fprintf(raster2, "%f %d\n", tempo, i);
	    if(i>200 && i<301) fprintf(raster3, "%f %d\n", tempo, i);
	    if(i>300 && i<401) fprintf(raster4, "%f %d\n", tempo, i);
	    }
	   
	     if(tempo>n*h-1000)
			{
			Tmedio[i]=Tmedio[i]+(tud[i]-tudant[i]);
			contT[i]=contT[i]+1;
		    }
	     	     
	     if(imprimir == 1)
		   if(tempo> h*n-2*transient)         // raster plot
		     fprintf(q,"%f %d\n",tud[i],i);	
	     
	     if(tempo>0.01)
	       { 
		 k[i]=k[i]+1;     
		 nk[k[i]][i]=tempo;  
		 
		 if(tempo<100.0)
		   if(k[i]==1.0)
		     kinicial=tempo;
	       }					  
	     kmax[i]=kmax[i]+1;
	     
	     pico[i]=1000;   
	     
	     if(Vr[i]==20.0)   
	       {
		 for(j=1;j<=conextotal[i];j++) 	
		   if(t2+delay[listaADJ[i][j]][i]<=delayMax)
		     SynMenexc[listaADJ[i][j]][t2+delay[listaADJ[i][j]][i]]=SynMenexc[listaADJ[i][j]][t2+delay[listaADJ[i][j]][i]]+listaconex[listaADJ[i][j]][i];  	                  
		   else
		     SynMenexc[listaADJ[i][j]][(t2+delay[listaADJ[i][j]][i])-delayMax]=SynMenexc[listaADJ[i][j]][(t2+delay[listaADJ[i][j]][i])-delayMax]+listaconex[listaADJ[i][j]][i];  	                  
	       }
	     else 
	       {						  
		 for(j=1;j<=conextotal[i];j++) 
		   if(t2+delay[listaADJ[i][j]][i]<=delayMax)
		     SynMenini[listaADJ[i][j]][t2+delay[listaADJ[i][j]][i]]=SynMenini[listaADJ[i][j]][t2+delay[listaADJ[i][j]][i]]+listaconex[listaADJ[i][j]][i];  	                  
		   else
		     SynMenini[listaADJ[i][j]][(t2+delay[listaADJ[i][j]][i])-delayMax]=SynMenini[listaADJ[i][j]][(t2+delay[listaADJ[i][j]][i])-delayMax]+listaconex[listaADJ[i][j]][i];  
	       }					  
	   }	
	 
	 if(kmax[i]>=2)
	   aux=0.0;  
	 
	 xant[i]=x[1+(i-1)*N];	
	 
       } //ending of neuron i loop
            
     if(imprimir == 1)
	   if(tempo> h*n-2*transient)  
	     {	
	       GexcM=0.0;
	       GiniM=0.0;	
	       Imedio=0.0;		
	       for(i=1;i<=NN;i++)           
		 {
		   GexcM=GexcM+Gexc[i]*(20.0-x[1+(i-1)*N]);
		   GiniM=GiniM+Gini[i]*(-75.0-x[1+(i-1)*N]);
		   Imedio=GiniM+GexcM;
		 }
	       
	       fprintf(q4,"%f %f %f %f %f %f %f %f \n",tempo, Gexc[10],Gini[10],Gexc[10]*(20.0-x[10]),Gini[10]*(-75.0-x[10]),GexcM/NN,GiniM/NN, Imedio/NN);	
	     }
     // ------------ Runge-Kutta --------------- 
     derivs(y,df,Gexc,Gini,I);
     for(i=1;i<=N*NN;i++)
       {
	 a[i]=h*df[i];
	 y[i]=x[i]+a[i]/2.0;
       }
     
     derivs(y,df,Gexc,Gini,I);
     for(i=1;i<=N*NN;i++)
       {
	 b[i]=h*df[i];
	 y[i]=x[i]+b[i]/2.0;
       }
     
     derivs(y,df,Gexc,Gini,I);
     for(i=1;i<=N*NN;i++)
       {
	 c[i]=h*df[i];
	 y[i]=x[i]+c[i]; 
       }
     
     derivs(y,df,Gexc,Gini,I);
     for(i=1;i<=N*NN;i++)
       x[i]=x[i]+(a[i]+h*df[i])/6.0+(b[i]+c[i])/3.0; 	  
     
     //********************* print coupling, mean degree and deltat times****************************************
     if(t==chave2)
       {
	 if(tempo>0.01)  // transient+passoAD*h)
	   {  
	     cont=0;
	     sum=0.0;
	     for(i=1;i<=NN;i++)
	       for(j=1;j<=NN;j++)
		 {
		   if(listaconex[i][j]>0)
		     if(Vr[j]==20.0)
		       {
			 sum=sum+listaconex[i][j];
			 cont=cont+1;	
		       }
		 }
	     
	     if(imprimir == 1)
		   fprintf(o,"%f %f ",tempo,sum/cont); 
	     
	     if(tempo>=TRANSIENTE_KMEDIO) 
	       {       
		 gex_medio = gex_medio+sum/cont;
		 g_cont = g_cont + 1.0;
	       }
	     
	     cont2=0;
	     sum=0.0;
	     for(i=1;i<=NN;i++)
	       for(j=1;j<=NN;j++)
		 {
		   if(listaconex[i][j]>0)
		     if(Vr[j]==-75.0)
		       {	
			 sum=sum+listaconex[i][j];
			 cont2=cont2+1;		      
		       }
		 }  
	     
	     aux3=(cont+cont2)/NN;
	     
	     if(imprimir == 1)
		   fprintf(o,"%f  %f ",sum/cont2,aux3); 
	     
	     
	     if(tempo>=TRANSIENTE_KMEDIO) 
	       {
		 ginib_medio = ginib_medio + sum/cont2;
		 kmedio = kmedio + aux3;
	       }
	        
	     if(imprimir == 1)
		   fprintf(o,"%f %f %f\n",deltatmedio/contdeltat,deltatmedioE/contdeltatE,deltatmedioI/contdeltatI); 
	     
	     deltatmedio=0.0;
	     contdeltat=0;	
	     deltatmedioE=0.0;
	     contdeltatE=0;	
	     deltatmedioI=0.0;
	     contdeltatI=0;	
	   }
	 chave2=chave2+passoAD;	       
       }  	
     //----------------------Plasticity-----------------------------------------
     if(tempo>ligar)   
     if(aux<1) 
       for(i=1;i<=NN;i++) 
	 if(tempo>=tud[i] && tempo<(tud[i]+h))
	   {				
	     //---------------------- eSTDP and iSTDP-----------------------------------------
	     for(j=1;j<=NN;j++) 		
	       if(ADJ[i][j]==1.0)  
		 {
		   //---------------------potentiation from neuron j to i -------------------//
		   ///////////eSTDP			
		   if(Vr[j]==20.0)   ///////////////// se o pré é excitatório		
		     if(kmax[j]>0) 
		       {
			 Deltat=tud[i]-tud[j];
			 
			 if(Deltat!=0)
			   { 
			     if(Deltat<0) 
			       Deltat=tud[i]-tudant[j];
			       
			     W1=A1*exp(-Deltat/tau1); // Potentiation
			     W2=-A2*exp(-Deltat/tau2); // Depression 
			     
			     listaconex[i][j]=listaconex[i][j]+delta_exc*W1; 
			     
			     if(Deltat<30.0)     
			       {
				 deltatmedio=deltatmedio+Deltat;
				 contdeltat=contdeltat+1;
				 
				 if(Deltat<10.0)    
				   {
				     deltatmedioE=deltatmedioE+Deltat;
				     contdeltatE=contdeltatE+1;
				   } 
			       }
			     
			     if(listaconex[i][j]>limite_sup_exc) 
			       listaconex[i][j]=limite_sup_exc;  
			     
			     if(listaconex[i][j]<0.0)  
			       listaconex[i][j]=0.0;
			   }  
		       }
		   //////////iSTDP
		   if(Vr[j]==-75.0) 
		     if(kmax[j]>0) 
		       {
			 Deltat=tud[i]-tud[j];
			 
			 if(Deltat!=0)
			   { 
			     if(Deltat<0) 
			       Deltat=tud[i]-tudant[j];
			     
			     W1_inib= (g/gnorm)*pow(alpha_plus,beta)*fabs(Deltat)*pow(Deltat,(beta-1.0))*exp(-alpha_plus*(fabs(Deltat)));  // Potentiation
			     
			     W2_inib=(g/gnorm)*pow(alpha_minus,beta)*fabs(Deltat)*pow((-Deltat),(beta-1.0))*exp(-alpha_minus*(fabs(Deltat)));  // Depression
			     
			     listaconex[i][j]=listaconex[i][j]+delta_inib*W1_inib;  
			     
			     if(Deltat<30.0)   
			       {
				 deltatmedio=deltatmedio+Deltat;
				 contdeltat=contdeltat+1;
				 
				 if(Deltat>2.5)   
				   {
				     deltatmedioI=deltatmedioI+Deltat;
				     contdeltatI=contdeltatI+1;
				   }
			       }
			     if(listaconex[i][j]>limite_sup_inib) 
			       listaconex[i][j]=limite_sup_inib;
			     
			     if(listaconex[i][j]<0.0)
			       listaconex[i][j]=0.0; 
			   }
		       } 
		   //------------------------------depression from neuron i to j---------------------//
		   /////////eSTDP
		   if(Deltat!=0 && Vr[i]==20.0)//caso i excitatorio
		     if(ADJ[j][i]==1)    // verifica se a conexão existe
		       {
			 if(Deltat<30.0)          
			   {
			     deltatmedio=deltatmedio+Deltat;
			     contdeltat=contdeltat+1;
			     
			     if(Deltat<10.0)
			       {
				 deltatmedioE=deltatmedioE+Deltat;
				 contdeltatE=contdeltatE+1;
			       }
			   }   
			 W2=-A2*exp(-Deltat/tau2); 
			 
			 listaconex[j][i]=listaconex[j][i]+delta_exc*W2;
			 
			 if(listaconex[j][i]>limite_sup_exc)  
			   listaconex[j][i]=limite_sup_exc;
			 
			 if(listaconex[j][i]<0.0)    
			   listaconex[j][i]=0.0;
		       }	
		   
		   //////////iSTDP
		   if(Deltat!=0 && Vr[i]==-75.0) //caso i inibitorio
		     if(ADJ[j][i]==1)  
		       {
			
			 if(Deltat<30.0)   
			   {
			     deltatmedio=deltatmedio+Deltat;
			     contdeltat=contdeltat+1;
			     
			     if(Deltat>2.5)       
			       {
				 deltatmedioI=deltatmedioI+Deltat;
				 contdeltatI=contdeltatI+1;
			       }
			   }
			 
			 W2_inib=(g/gnorm)*pow(alpha_minus,beta)*fabs(Deltat)*pow((-Deltat),(beta-1.0))*exp(-alpha_minus*(fabs(Deltat))); 
			
			 listaconex[j][i]=listaconex[j][i]+delta_inib*W2_inib;
			 		 
			 if(listaconex[j][i]>limite_sup_inib)
			   listaconex[j][i]=limite_sup_inib;
			 
			 if(listaconex[j][i]<0.0)
			   listaconex[j][i]=0.0;
		       }  
		 } 
	   } // ending of stdp loop
     //******************save the layers of coupling matrix ********************************
     if(imprimir == 1)
	   if(t==chave)
	     {
	       sum=0.0;
	       for(i=1;i<=NN;i++)
		 for(j=1;j<=NN;j++)
		   {
		     Matriz[i][j][contM]=listaconex[i][j];
		     if(Vr[j]==-75.0)
		       Matriz[i][j][contM]=-listaconex[i][j];		  
		   }
	       chave=chave+passoM;
	       contM=contM+1;	    
	     }  
     
     if(t2==delayMax) // restart the delay matrix cycle
       t2=0;       
   }  // ending of timing loop
 
 /////////////////////////print the mean weight of the connections between the subnetworks
 if(subredes>1)     
   if(imprimir == 1)
	 for(i=1;i<=subredes;i++) 
	   for(j=1;j<=subredes;j++)
	     {
	       acop_medioex[i][j]=0.0;
	       contaux[i][j] = 0; 
	       
	       
		   for(aux=rede*(i-1)+1;aux<=rede*(i-1)+rede;aux=aux+1)
		     for(aux2=rede*(j-1)+1;aux2<=rede*(j-1)+rede;aux2=aux2+1)
		       if(ADJ[aux][aux2]==1)
		   	    {	      	
			     acop_medioex[i][j] = acop_medioex[i][j]+listaconex[aux][aux2];
			     contaux[i][j] = contaux[i][j]+1;
			    }
		 
		 acop_medioex[i][j] = acop_medioex[i][j]/contaux[i][j];	 
	     }
 for(i=1;i<=subredes;i++) 
	   for(j=1;j<=subredes;j++)
			{
			 if(i==j)	
			     fprintf(coup,"%d %d %f %f\n",i,j, 0.0, 0.0); 
			 else			 
		        if(acop_medioex[i][j]>=acop_medioex[j][i])
				   fprintf(coup,"%d %d %f %f\n",i,j,acop_medioex[i][j],acop_medioex[i][j]-acop_medioex[j][i]);
			    else
				   fprintf(coup,"%d %d %f %f\n",i,j,acop_medioex[i][j],0.0);
			}
 //******************print matrix ********************************
 if(imprimir == 1)
       for(i=1;i<=NN;i++)
	 for(j=1;j<=NN;j++)
	   {
	     fprintf(p,"%d %d ",i,j);   // i=pós ,  j=pré   
	     
	     if(ADJ[i][j]==1)  
	       {		
		 for(contM=1;contM<=Ngif;contM++)
		   fprintf(p," %.3f",Matriz[i][j][contM]);  
		 
	         fprintf(p,"\n");
	       }	  
	     else
	       {
		 for(contM=1;contM<=Ngif;contM++)
		   fprintf(p," %.3f",0.000); 
		 fprintf(p,"\n"); 
	       }
	   }	  

 // Parameter of Measure
 chave3=0.01;
 tempo=0.0;
 kfinal=n*h;
 
 for(m=1;m<=subredes;m=m+1) 
   {
     Rmedio[m]=0.0;
     contR[m]=0.0;
   }
 
 for(i=1;i<=NN;i=i+1)
   {
     phi[i]=0.0;
     kmax[i]=k[i];
     k[i]=1.0;	
     
     if(kmax[i]>1)    
       if(kfinal>nk[kmax[i]][i] && kfinal>kinicial)
	 kfinal=nk[kmax[i]][i];
   }
 
 chave4=kinicial;
 chave5=kfinal-45.0;  
 for(tempo=0.01;tempo<=n*h;tempo=tempo+1*h)  //
   {
     for(i=1;i<=NN;i=i+1)
       if(kmax[i]>1)
	 if(tempo>nk[k[i]][i] && tempo<nk[kmax[i]][i])
	   {
	     if(tempo<=nk[k[i]+1][i])	  
	       phi[i]=2*PI*(tempo-nk[k[i]][i])/(nk[k[i]+1][i]-nk[k[i]][i]); //+2*PI*k[i]
	     
	     if(tempo>=nk[k[i]+1][i])
	       k[i]=k[i]+1;		
	   }
     ////////////////////////////////////////////////////  print phase on file
     if(imprimir == 1)
	   {
	     if(tempo>=kinicial && tempo<=kinicial+40.0)   
	       if(tempo>=chave4)
		 {	 
		   for(i=1;i<=NN;i++)	    
		     fprintf(fases,"%f\t%d\t%f\n", tempo, i, phi[i]);
		   
		   chave4 = chave4+0.1;
		 }
	     
	     if(tempo>=kfinal-45.0 && tempo<=kfinal-5.0)
	       if(tempo>=chave5)
		 {
		   for(i=1;i<=NN;i++)  
		     fprintf(fases,"%f\t%d\t%f\n", tempo, i, phi[i]);
	           
		   chave5=chave5+0.1;
		 }
	   }	   
/////////////////////////////////////////////////
     if(tempo>=kfinal-1000 && tempo<=kfinal)
       {
	 
	 for(m=1;m<=subredes;m=m+1)  //for(m=1;m<=subredes;m=m+1) 
	   {
	     real=0.0;
	     compl=0.0;
	     cont=0;
	     
	     for(i=1;i<=NN;i++)
	       if(kmax[i]>1)
		 {
		   real=real+cos(m*phi[i]);
		   compl=compl+sin(m*phi[i]); 
		   
		   cont=cont+1;
		 }
	     
	     real=real/cont;
	     compl=compl/cont;
	     
	     R[m]=sqrt(real*real+compl*compl);
	     
	     
	     if(tempo>=kfinal-1000)
	       {
		 Rmedio[m]=Rmedio[m]+R[m];
		 contR[m]=contR[m]+1.0;
	       }
	   }
	   
	    if(imprimir == 1)
	    {
	       fprintf(q3,"%f ",tempo);	
	       
	       for(m=1;m<=subredes;m=m+1)   //for(m=1;m<=subredes;m=m+1)
		        fprintf(q3,"%f ",R[m]);
	       
	       fprintf(q3,"\n");
	       
	       chave3=chave3+10.0;	      
	     }	        
       } 
   }
 } 
 
  T=0.0;
  for(i=1;i<=NN;i++)	
	T=T+Tmedio[i]/contT[i]; 
 
/////////////  routine to printf delays and higher moment order parameter
 fprintf(q2,"%f %f ",delay_medio, delay_medio_inter);	
 
 aux1=Rmedio[1]/contR[1];
 aux2=1;
 aux3=0.0; 	 
 for(m=2;m<=subredes;m=m+1)
   {
     aux3=Rmedio[m]/contR[m];
     
     if(aux3>aux1)
       {
	 aux1=aux3;
	 aux2=m;
       } 			 
   }
 
 fprintf(q2,"%d %f ",aux2, aux1);
 
 for(m=1;m<=subredes;m=m+1)
   fprintf(q2,"%f ", Rmedio[m]/contR[m]);
 
 fprintf(q2,"%f \n", T/NN);   // printed file: delay_intra, delay_inter, maior m, R_medio(maior m), R_medio(m=1), R_medio(m=2) ... periodo medio da rede 
 //////////////////////
 
 GEX= GEX + gex_medio/g_cont;
 GEX2= GEX2+(gex_medio/g_cont)*(gex_medio/g_cont);
 GINIB = GINIB + ginib_medio/g_cont;
 GINIB2 = GINIB2 + (ginib_medio/g_cont)*(ginib_medio/g_cont);
 K=K+kmedio/g_cont;
 K2=K2+(kmedio/g_cont)*(kmedio/g_cont);
 
 } //ending delay loops
 
 GEX= GEX/repeticoes;
 GEX2=sqrt(GEX2/repeticoes - GEX*GEX);
 GINIB = GINIB/repeticoes;
 GINIB2 = sqrt(GINIB2/repeticoes - GINIB*GINIB);
 K=K/repeticoes;
 K2= sqrt(K2/repeticoes - K*K);
 free_dmatrix(listaconex,1,NN+1,1,NN+1);
 free_imatrix(listaADJ,1,NN+2,1,conexMAX+2);
 free_imatrix(ADJ,1,NN+2,1,conexMAX+2);
 free_vector(conextotal,1,NN+1); 
 free_dvector(y,1,N*NN+1);
 free_dvector(df,1,N*NN+1); 
 free_dvector(I,1,NN+1);
 free_f3tensor(Matriz,1,NN+2,1,NN+2,1,Ngif+2);
 free_f3tensor(MatrizI,1,NN+2,1,NN+2,1,Ngif+2);
 
 free_vector(k,1,NN+1); 
 free_vector(kmax,1,NN+1);
 free_dvector(phi,1,NN+1); 
 free_dmatrix(nk,1,NMD+2,1,NN+2);
 free_dvector(Gexc,1,NN+1); 
 free_dvector(Gini,1,NN+1); 
 free_dmatrix(SynMenexc,1,NN+2,1,delayMax+2);
 free_dmatrix(SynMenini,1,NN+2,1,delayMax+2);
 free_imatrix(delay,1,NN+2,1,NN+2);
 fclose(o);
 fclose(p);
 fclose(q);
 fclose(q2);
 fclose(q3);
 fclose(q4);
 fclose(fases);
 fclose(coup);
 
  return 0;
}

void derivs(double y[],double df[], double *Gexc, double *Gini,double *I)
{
  int i;
  double gna,gL,gk,vna,vk,vL,CM,am,bm,ah,bh,an,bn,Vr,Vr_inib;
  
  gna=120.0;
  gk=36.0;
  gL=0.3;
  vna=50.0;
  vk=-77.0;
  vL=-54.4; 
  CM=1.0;
  Vr=20.0;    
  Vr_inib=-75.0; 
  
  for(i=1;i<=NN;i=i+1)
    {
      am=0.1*(y[1+(i-1)*N]+40.0)/(1.0-exp((-y[1+(i-1)*N]-40.0)/10.0));
      bm=4.0*exp((-y[1+(i-1)*N]-65.0)/18.0);
      ah=0.07*exp((-y[1+(i-1)*N]-65.0)/20.0);
      bh=1.0/(1.0+exp((-y[1+(i-1)*N]-35.0)/10.0));
      an=0.01*(y[1+(i-1)*N]+55.0)/(1.0-exp((-y[1+(i-1)*N]-55.0)/10.0));
      bn=0.125*exp((-y[1+(i-1)*N]-65.0)/80.0);
      
      df[1+(i-1)*N]=(1.0/CM)*(I[i]-gna*(pow(y[2+(i-1)*N],3)*y[3+(i-1)*N]*(y[1+(i-1)*N]-vna))
			      -gk*(pow(y[4+(i-1)*N],4)*(y[1+(i-1)*N]-vk))-gL*(y[1+(i-1)*N]-vL)) + Gexc[i]*(Vr-y[1+(i-1)*N]) + Gini[i]*(Vr_inib-y[1+(i-1)*N]);
      df[2+(i-1)*N]=am*(1.0-y[2+(i-1)*N])-bm*y[2+(i-1)*N];
      df[3+(i-1)*N]=ah*(1.0-y[3+(i-1)*N])-bh*y[3+(i-1)*N];
      df[4+(i-1)*N]=an*(1.0-y[4+(i-1)*N])-bn*y[4+(i-1)*N];      
    } 
}

float ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  
  if(*idum<=0 || !iy)
    {
     if(-(*idum)<1) *idum=1;
     else *idum = -(*idum);
    for(j=NTAB+7;j>=0;j--)
      {
       k=(*idum)/IQ;
       *idum=IA*(*idum-k*IQ)-IR*k;
       if(*idum<0) *idum +=IM;
       if(j<NTAB) iv[j]=*idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQ;
   *idum=IA*(*idum-k*IQ)-IR*k;
   if(*idum<0) *idum += IM;
   j=iy/NDIV;
   iy=iv[j];
   iv[j]=*idum;
   if((temp=AM*iy)>RNMX) return RNMX;
   else return temp;
}

double *dvector(long nl,long nh)
{
   double *v;
   
   v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
   if (!v) nrerror("allocation failure in dvector()");
   return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh)
{
   free((FREE_ARG) (v+nl-NR_END));
}

void nrerror(char error_text[])
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

int *vector(long nl,long nh)
{
   int *v;
   
   v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
   if (!v) nrerror("allocation failure in dvector()");
   return v-nl+NR_END;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
{
   long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   int **m;

   m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
   if (!m) nrerror("allocation failure 1 in imatrix()");
   m += NR_END;
   m -= nrl;

   m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
   if (!m[nrl]) nrerror("allocation failure 2 in imatrix()");
   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   return m;
}

void free_vector(int *v, long nl, long nh)
{
   free((FREE_ARG) (v+nl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
   free((FREE_ARG) (m[nrl]+ncl-NR_END));
   free((FREE_ARG) (m+nrl-NR_END));
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
   long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   double **m;

   m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
   if (!m) nrerror("allocation failure 1 in matrix()");
   m += NR_END;
   m -= nrl;

   m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
   if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
   free((FREE_ARG) (m[nrl]+ncl-NR_END));
   free((FREE_ARG) (m+nrl-NR_END));
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;

	/* allocate pointers to pointers to rows */
	t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

float gasdev(long *idum)
{
  float ran1(long *idum);
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;  
  
  if(*idum<0) iset=0;
  if(iset==0) {
    do {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq>=1.0 || rsq==0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
} 
