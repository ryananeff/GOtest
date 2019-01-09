#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <math.h>
int sumi(int *x, int n)
{
	int result;
	result=0;
	for (int i=0; i < n; i++){
		result += x[i];
	}
	return(result);
}
void fGSEAfunc2(int *tag, double *strength, int *n, double *alpha, double *ab){
	int Nh=sumi(tag,*n);
	double Nm;
	double ys[Nh],f2[Nh];
	double a,b,a1,b1,t1;
	if(Nh<2) return;
	Nm=(*n-Nh+0.0);
	int j=0,i0=-1,i1=0;
	double ysSum=0.0;
	for(int i=0; i < *n; i++){
		if(tag[i]==0) continue;
		i1 = i-i0-1;
		f2[j] = i1/Nm;
		i0=i;
		if(fabs(*alpha) < 0.01){ //0.01 is an arbitary threshold
			ys[j]=1.0;
		}else{
			t1=fabs(strength[i]);
			//Rprintf("# = %d, v1 = %f, v2 = %f\n", i , strength[i],t1 );
			if(fabs(*alpha-1.0) > 0.01) t1=pow(t1, *alpha);
			ys[j]=t1;
			ysSum+=t1;
		}
		j++;
	}
	ys[0]=ys[0]/ysSum;
	a=ys[0]-f2[0];
	a1=a;
	b=a-ys[0];
	b1=b;
	for(int i=1; i < Nh; i++){
		ys[i]=ys[i]/ysSum;
		a1+=ys[i]-f2[i];
		//Rprintf("%f ",ys[i]);
		//if( i % 5 ==0) Rprintf("\n");
		if(a1>a) a=a1;
		b1=a1-ys[i];
		if(b1<b) b=b1;
	}
	//Rprintf("\nv1 = %f, v2 = %f\n", a,b);
	*ab= (a + b) > 0 ? a : b;
	return;
}
