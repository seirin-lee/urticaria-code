#include <stdio.h>
#include <stdlib.h>
#include <math.h>


///2arrarys Definition ----------------------------------
#ifdef __cplusplus
template<typename T>T**AllocMatrix(int u,int v)
{
    int i; T**a,*b;
    try { a=(T**)new char[(sizeof*a+sizeof*b*v)*u]; }
    catch (...) { a=0; }
    if (a) b=(T*)(a+u); else return 0;
    for (i=1;i<u;i++,b+=v) a[i]=b;
    return a;
}
#define ALLOC_MATRIX(T,U,V) AllocMatrix<T>(U,V)
#define FREE(X) delete[]X

#else

void*AllocMatrix(int s,int u,int v)
{
    int i,t=s*v; char**a,*b;

    a=(char**)malloc((sizeof*a+t)*u);
    if (a) b=(char*)(a+u); else return 0;
    for (i=1;i<u;i++,b+=t) a[i]=b;
    return a;
}
#define ALLOC_MATRIX(T,U,V) (T**)AllocMatrix(sizeof(T),U,V)
#define FREE(X) free(X)
#endif
//------------------------------------------------------------------


#define N 256  //Spatial grid number
#define M 4000 //Time
#define m 100 //temporal grid number
#define Nn 256 //plot grid number
#define Mn 20 //plot time number


//Diffusion coefficients
#define D1x 0.0000047
#define D1y 0.0000047
#define D4x 0.0000047
#define D4y 0.0000047


//Total histamin of Mast cell(U0) and Basophil(U0B)
#define U0 (600.0)
#define U0B (600.0)


/* Kinetic Paramaeters */
//Baosphil-related equation
#define delB (0.1)
#define gaB (5.2)
#define alpB (0.335)
#define alpB0 (0.00625)
#define muB (1.0)

//Tissue factor-related equation
#define delT (0.01)
#define gaT (7.0)
#define gaT0 (2.2)
#define alpT (1.0)
#define alpT0 (8.925)
#define muT (1.0)

//Mast cell-related equation
#define delM (0.01)
#define gaM (5.0)
#define alpM (0.865)
#define alpM0 (0.0035)
#define muM (1.0)

//Coagulation factor-related equation
#define gaC (20.0)
#define muC (1.0)
#define beta (20.0)


//Parameter for initial conditions
#define p_rangeB (0.01)
#define pp (0.9995)
#define p_rangeH (0.01)
#define u10 (delM/muM)
#define u20 (delT/muT)
#define u200 (0.67)
#define u30 (delB/muB)
#define u300 (u30*(1+p_rangeB*pp))
#define u40 (0.0)
//Wheal state function(Skin=S_w) parameter
#define sp (10.0)
#define sp2 (1.0)




#define pi 3.14159265359

void initial (double **u1, double **u2, double **u3, double **u4);
void output(double **u1, double **u2, double **u3, double **u4, double **u1_sum, int k);
void RungeKutta(double **u1, double **u2, double **u3, double **u4, double **u1_sum, double **u3_sum);

/*ADI法*/
void ReactionTerm(double **u1, double **u2, double **u3, double **u4, double **u1_sum, double **u3_sum, double **F1, double **F4);
void RightSide1(double **u1, double **F1, double alp1y, double **u4, double **F4, double alp4y);
void RightSide2(double **u1, double **F1, double alp1x, double **u4, double **F4, double alp4x);
void LU1(double **u1, double alp1x, double **u4, double alp4x);
void LU2(double **u1, double alp1y, double **u4, double alp4y);

double Rac1(double u1, double u2, double u3, double u4, double KaiM);
double Rac2(double u1, double u2, double u3, double u4, double KaiM, double KaiB);
double Rac3(double u1, double u2, double u3, double u4, double KaiM, double KaiB);
double Rac4(double u1, double u2, double u3, double u4, double KaiM, double KaiB);


double dx = 1.0/N;     
double dy = 1.0/N;
double dt = 1.0/m;      


int main()
{  FILE *fp;
	double **u1, **u2, **u3, **u4, **u1_sum, **u3_sum, **F1, **F4;
	double skin, KaiM, KaiB, Kai;
	int k,  i,  j;
	double alp1x = (D1x/2.0) * dt/(dx * dx);
	double alp1y = (D1y/2.0) * dt/(dy * dy);
    double alp4x = (D4x/2.0) * dt/(dx * dx);
    double alp4y = (D4y/2.0) * dt/(dy * dy);
    double kk1, kk2, kk3, kk4, ww1, ww2, ww3, ww4;
    

	
	
	// Initializing u1 and u2 arrays
    u1 = ALLOC_MATRIX(double, N,  N);
    u2 = ALLOC_MATRIX(double, N,  N);
    u3 = ALLOC_MATRIX(double, N,  N);
    u4 = ALLOC_MATRIX(double, N,  N);
    
    u1_sum = ALLOC_MATRIX(double, N,  N);
    u3_sum = ALLOC_MATRIX(double, N,  N);
    F1= ALLOC_MATRIX(double, N,  N);
    F4= ALLOC_MATRIX(double, N,  N);



	

    
	for (i=1 ; i<=N-1 ; ++i ) {
		for (j=1 ; j<=N-1 ; j++) {
			u1_sum[i][j]=0.0;
            u3_sum[i][j]=0.0;
								}
    }
    
    
    initial(u1, u2, u3, u4);
    k=0;
    output(u1, u2, u3, u4, u1_sum, k);

    
    
    for (k=1 ; k<=M ; k++){

        ReactionTerm(u1, u2, u3, u4, u1_sum,u3_sum,F1, F4);
        RungeKutta(u1, u2, u3, u4, u1_sum, u3_sum);
        RightSide1(u1, F1, alp1y, u4, F4, alp4y);
        LU1(u1, alp1x, u4, alp4x);
        RightSide2(u1, F1, alp1x, u4, F4, alp4x);
        LU2(u1, alp1y, u4, alp4y);

	
    if (k%(M/Mn)==0){
            output(u1, u2, u3, u4, u1_sum, k);
		}
		
		
		
       
    for (i=1 ; i<=N-1 ; i++ ) {
        for (j=1 ; j<=N-1 ; j++){
            u1_sum[i][j]=u1[i][j]+u1_sum[i][j];
            u3_sum[i][j]=u3[i][j]+u3_sum[i][j];
        }
    }
		
 
        
	}//k-end


	
	
// Cleaning memory
    FREE(u1);
    FREE(u2);
    FREE(u3);
    FREE(F1);
    FREE(u4);
    FREE(F4);
    FREE(u1_sum);
	FREE(u3_sum);

	
	return 0;
}

//===============END=============================








//============ FUNCTIONS===========================

void initial (double **u1, double **u2, double **u3, double **u4 ){
    int i, j ;
    double x, y;
    
    
    
for (i=1 ; i<= N-1 ; i++) {
		for (j=1; j<=N-1 ; j++){
            double eexf=(double)rand()/RAND_MAX;
            x=i*dx; y=j*dy;
            u1[i][j] = u10*(1.0+p_rangeH*eexf);
            u2[i][j] = u20*(1.0+0.01*eexf);
            u3[i][j] = u30*(1.0+p_rangeB*eexf);
            if(u3[i][j]<u300){
                u3[i][j]=0.0;
            }
            u4[i][j]=u40*(1.0+0.0*eexf);
            
        }
    }
}			


void output(double **u1, double **u2, double **u3, double **u4, double **u1_sum, int k){
	int i, j ;
    double skin;
    
	for (i=1 ; i<= N-1 ; i++) {
            for (j=1; j<=N-1 ; j++){
                skin=1.0/(1.0+exp(-sp2*(u4[i][j]-sp)));
    printf(" %f %f %f %f %f %f %f %f\n",  dt*k,  i*dx, j*dy, u1[i][j], u2[i][j], u3[i][j], u4[i][j], skin);
    
            }
         printf("\n");
        }
     printf("\n");

}


void RungeKutta(double **u1, double **u2, double **u3, double **u4, double **u1_sum, double **u3_sum){
    int i, j;
    double KaiM, KaiB;
    double kk1, kk2, kk3, kk4, ww1, ww2, ww3, ww4;

    for (i=1 ; i<=N-1 ; i++ ) {
                for (j=1 ; j<=N-1 ; j++) {
                  
                    if(u1_sum[i][j]<U0){
                        KaiM=1.0;
                    }
                    else{
                        KaiM=0.0;
                    }
                    
                    if(u3_sum[i][j]<U0B){
                        KaiB=1.0;
                    }
                    else{
                        KaiB=0.0;
                    }
                    
                    
                    kk1=Rac2(u1[i][j], u2[i][j], u3[i][j], u4[i][j], KaiM, KaiB);
                    kk2=Rac2(u1[i][j], u2[i][j]+kk1*dt*0.5, u3[i][j]+ww1*dt*0.5, u4[i][j], KaiM, KaiB);
                    kk3=Rac2(u1[i][j], u2[i][j]+kk2*dt*0.5, u3[i][j]+ww2*dt*0.5, u4[i][j], KaiM, KaiB);
                    kk4=Rac2(u1[i][j], u2[i][j]+kk3*dt, u3[i][j]+ww3*dt, u4[i][j], KaiM, KaiB);
                    
                    
                    if(u3[i][j]<u300){
                        u3[i][j]=0.0;
                    }
                    else{
                    ww1=Rac3(u1[i][j], u2[i][j], u3[i][j], u4[i][j], KaiM, KaiB);
                    ww2=Rac3(u1[i][j], u2[i][j]+kk1*dt*0.5, u3[i][j]+ww1*dt*0.5, u4[i][j], KaiM, KaiB);
                    ww3=Rac3(u1[i][j], u2[i][j]+kk2*dt*0.5, u3[i][j]+ww2*dt*0.5, u4[i][j], KaiM, KaiB);
                    ww4=Rac3(u1[i][j], u2[i][j]+kk3*dt, u3[i][j]+ww3*dt, u4[i][j], KaiM, KaiB);
                    
                    u3[i][j]=u3[i][j]+dt*(ww1+2.0*ww2+2.0*ww3+ww4)/6.0;
                    }
                                        
                    u2[i][j]=u2[i][j]+dt*(kk1+2.0*kk2+2.0*kk3+kk4)/6.0;
                }
        
    }
}

void ReactionTerm(double **u1, double **u2, double **u3, double **u4, double **u1_sum, double **u3_sum, double **F1, double **F4){
	int i, j;
    double KaiM, KaiB;
    
        
    
	for (i=1 ; i<=N-1 ; i++ ) {	
			    for (j=1 ; j<=N-1 ; j++) {
                  if(u1_sum[i][j]<U0){
                        KaiM=1.0;
                    }
                    else{
                        KaiM=0.0;
                    }
                    if(u3_sum[i][j]<U0B){
                            KaiB=1.0;
                    }
                    else{
                            KaiB=0.0;
                    }
                   
                    
                    F1[i][j]= (dt/2.0)*Rac1(u1[i][j], u2[i][j], u3[i][j], u4[i][j], KaiM);
                    F4[i][j]= (dt/2.0)*Rac4(u1[i][j], u2[i][j], u3[i][j], u4[i][j], KaiM, KaiB);
                }
        }
}



//Functions for ADI method

void RightSide1(double **u1, double **F1,  double alp1y, double **u4, double **F4,  double alp4y){
	int i, j ;
	double w1[N], w4[N];
	for (i=1 ; i<=N-1 ; i++) {	
		    for (j=1 ; j<=N-1 ; j++ ) {	
                if(j==1) {
				w1[j]=alp1y*u1[i][j]+(1.0-2*alp1y)*u1[i][j]+alp1y*u1[i][j+1]+F1[i][j];
                w4[j]=alp4y*u4[i][j]+(1.0-2*alp4y)*u4[i][j]+alp4y*u4[i][j+1]+F4[i][j];
                    }
                else if (j==N-1) {
                    w1[j]=alp1y*u1[i][j-1]+(1.0-2*alp1y)*u1[i][j]+alp1y* u1[i][N-1]+F1[i][j];
                    w4[j]=alp4y*u4[i][j-1]+(1.0-2*alp4y)*u4[i][j]+alp4y* u4[i][N-1]+F4[i][j];
                }
                else {
                    w1[j]=alp1y*u1[i][j-1]+(1.0-2*alp1y)*u1[i][j]+alp1y * u1[i][j+1]+F1[i][j];
                    w4[j]=alp4y*u4[i][j-1]+(1.0-2*alp4y)*u4[i][j]+alp1y * u4[i][j+1]+F4[i][j];
                }
                            }
        for(j=1; j<N;j++){
            u1[i][j]=w1[j]; u4[i][j]=w4[j];
        }
    }
}

//Functions for ADI method
void RightSide2(double **u1, double **F1,  double alp1x, double **u4, double **F4,  double alp4x){
	int i, j;
	double w1[N], w4[N];
     		for (j=1 ; j<=N-1 ; j++) {	
                for (i=1 ; i<=N-1 ; i++ ) {
                    if(i==1) {
                        w1[i]=alp1x*u1[1][j]+(1.0-2*alp1x)*u1[i][j]+alp1x*u1[i+1][j]+F1[i][j];
                        w4[i]=alp4x*u4[1][j]+(1.0-2*alp4x)*u4[i][j]+alp4x*u4[i+1][j]+F4[i][j];
                            }
                    else if (i==N-1) {
                        w1[i]=alp1x*u1[i-1][j]+(1.0-2*alp1x)*u1[i][j]+alp1x * u1[N-1][j]+F1[i][j];
                        w4[i]=alp4x*u4[i-1][j]+(1.0-2*alp4x)*u4[i][j]+alp4x * u4[N-1][j]+F4[i][j];
                            }
                    else {
                        w1[i]=alp1x*u1[i-1][j]+(1.0-2*alp1x)*u1[i][j]+alp1x * u1[i+1][j]+F1[i][j];
                        w4[i]=alp4x*u4[i-1][j]+(1.0-2*alp4x)*u4[i][j]+alp4x * u4[i+1][j]+F4[i][j];
                                }
                }
                for(i=1; i<N;i++){
                        u1[i][j]=w1[i];  u4[i][j]=w4[i];
                    
                }
        }
}

//Functions for ADI method
void LU1(double **u1, double alp1x, double **u4, double alp4x){
	int i, j;
	double bet1,  gam1[N], bet4,  gam4[N];
	for (j=1; j<=N-1; j++){

        bet1=1.0 +alp1x;
        bet4=1.0 +alp4x;
		u1[1][j] = u1[1][j] /bet1;
        u4[1][j] = u4[1][j] /bet4;
		

		for  (i=2; i<=N-1 ; i++){
			gam1[i] = -alp1x /bet1 ;
			bet1 = 1.0 +2.0*alp1x + alp1x * gam1[i];
            gam4[i] = -alp4x /bet4 ;
            bet4 = 1.0 +2.0*alp4x + alp4x * gam4[i];
			
			if (i==N-1) {bet1 = bet1-alp1x;  bet4 = bet4-alp4x;}
			u1[i][j]=(u1[i][j] + alp1x*u1[i-1][j])/bet1;
            u4[i][j]=(u4[i][j] + alp4x*u4[i-1][j])/bet4;
		    }

		for (i=(N -2) ; i>=1 ; i--){	
			u1[i][j] =  u1[i][j]-gam1[i+1]*u1[i+1][j];
            u4[i][j] =  u4[i][j]-gam4[i+1]*u4[i+1][j];
		    
        }
			
    }
}

//Functions for ADI method
void LU2(double **u1,double alp1y, double **u4, double alp4y){
	int i, j;	
	double bet1, gam1[N], bet4, gam4[N];
for (i=1; i<=N-1; i++){
 			bet1=1.0 +alp1y;
            u1[i][1] = u1[i][1] /bet1;
            bet4=1.0 +alp4y;
            u4[i][1] = u4[i][1] /bet4;

		for  (j=2; j<=N-1 ; j++){
			gam1[j] = -alp1y /bet1 ;
			gam4[j] = -alp4y /bet4 ;
			bet1 = 1.0 +2.0*alp1y + alp1y * gam1[j];
			bet4 = 1.0 +2.0*alp4y + alp4y * gam4[j];
			if (j==N-1) {bet1 = bet1-alp1y; bet4 = bet4-alp4y; }
			u1[i][j]=(u1[i][j] + alp1y*u1[i][j-1])/bet1;
            u4[i][j]=(u4[i][j] + alp4y*u4[i][j-1])/bet4;
			
		}//j-loop end
		for (j=(N -2) ; j>=1 ; j--){				
			u1[i][j] =  u1[i][j]-gam1[j+1]*u1[i][j+1];
            u4[i][j] =  u4[i][j]-gam4[j+1]*u4[i][j+1];
					}
	}
}





//u1=Mast cell, u2=TF, u3=Basophil, u4=Coagulation factors

double Rac1( double u1, double u2, double u3, double u4, double KaiM )
{
   
    
	return delM+gaM*u4*KaiM*(1-alpM*pow(u1,2)/(alpM0+pow(u1,2)))-muM*u1;
  


}

double Rac2( double u1, double u2, double u3, double u4, double KaiM, double KaiB )
{   

    return delT+gaT*((u3+u1)/(gaT0+u3+u1))*(1-alpT*pow(u1+u3,2)/(alpT0+pow(u1+u3,2)))-muT*u2;

}

double Rac3( double u1, double u2, double u3, double u4, double KaiM, double KaiB )
{
   
    
		return delB+gaB*(u2)*KaiB*(1-alpB*pow(u3,2)/(alpB0+pow(u3,2)))-muB*u3;
  
}

double Rac4( double u1, double u2, double u3, double u4, double KaiM, double KaiB )
{
    
    
   
    return gaC/(1+exp(-beta*(u2-u200)))-muC*u4;

}
