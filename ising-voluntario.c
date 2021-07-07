//VOLUNTARIO MODELO DE ISING


# include <stdio.h>
# include <math.h>
#include "gsl_rng.h" 


# define PI 3.1415926535
# define N 32
gsl_rng *tau; 

//Definimos las funciones que vamos a utilizar. 

double incremento_energia(int n, int m, int s[N][N]);
void conf_inicial_1(int s[N][N]);
void conf_inicial_2(int s[N][N]);
double magnet(int s[N][N]);
double energia(int s[N][N]);
double correlacion(int s[N][N],int i);

int main()
{

	int i,j,k,n,m,paso,cont;  //Variables auxiliares
    double x,y,z;  //Variables auxiliares
    extern gsl_rng *tau; //Puntero al estado del número aleatorio
    int semilla=91284892; //Semilla del generador de números aleatorios
	int	s[N][N]; //s: matriz donde guardamos los espines.
	double delta_E,p;	// delta_E=diferencia de energía.
	double T; //temperatura del sistema
	double eps; 
	double mag; //magnetización
	double mag2; //magnetización al cuadrado
	double e_media; //energía media
	double e_media2; //energía al cuadrado
	double calor; //calor específico
	double cor[N]; //correlación
	double sigma_mag, sigma_e;  //varianza magnetización, varianza energía
	FILE *f1,*f2,*f3,*f4,*f5;

	//Abrimos ficheros para guardar los datos 

	f1=fopen("ising-32.txt","w"); 
	f2=fopen("magnet-32.txt","w"); 
	f3=fopen("energia-32.txt","w"); 
	f4=fopen("calor-32.txt","w"); 
	f5=fopen("correlacion-32.txt","w"); 

	fprintf(f2,"N\tT\tMagnetización\tError\n");
	fprintf(f3,"N\tT\tEnergía\tError\n");
	fprintf(f4,"N\tT\tCalor específico\n");
	fprintf(f5,"N\tT\tDistancia\tCorrelación\n");


	tau=gsl_rng_alloc(gsl_rng_taus); //Inicializamos el puntero
    gsl_rng_set(tau,semilla); //Inicializamos la semilla

	paso=pow(N,2);   //definimos el paso Montecarlo.

	for(T=1.5;T<=3.5;T=T+0.2)
	{
	
		//Definimos la configuración inicial de espines.
	
		conf_inicial_1(s);  //espines ordenados
	 
		//conf_inicial_2(s);    //espines desordenados de manera aleatoria

		//Igualamos a 0 las magnitudes que vamos a calcular. 

		mag=0;
		mag2=0;
		e_media=0;
		e_media2=0;

		for(j=0;j<N;j++)
		{
			cor[j]=0;
		}
	
		//Hacemos evolucionar el sistema 10^6 pasos MonteCarlo.
	
		for(k=0;k<pow(10,6)*paso;k++)
	
		{
			n=gsl_rng_uniform_int(tau,N);
			m=gsl_rng_uniform_int(tau,N); 
	

			delta_E=incremento_energia(n,m,s);
		
	
			if(exp(-delta_E/T)<=1)
			{
				p=exp(-delta_E/T);
			}
			else 
			{
				p=1;
			}
		
			eps=gsl_rng_uniform(tau);
	
	
			if(eps<p)
			{
				s[n][m]=-1*s[n][m];
			
			}

			//Calculamos la magnetización, energía media y correlación cada 100 pMC.
		
			if(k%(100*paso)==0)
			{
				mag=mag+magnet(s);
				mag2=mag2+pow(magnet(s),2);
				e_media=e_media+energia(s);
				e_media2=e_media2+pow(energia(s),2);
	
				for(i=0;i<N;i++)
				{
					cor[i]=cor[i]+correlacion(s,i);		
					cor2[i]=cor2[i]+pow(correlacion(s,i),2);			
				}

			
			}
		
		}
		
		//Calculamos los promedios de la magnetización, energía, calor específico y correlación. También calculamos sus respectivas varianzas para calcular los errores.

		mag=mag/pow(10,4);
		mag2=mag2/pow(10,4);

		e_media=e_media/(2*N*pow(10,4));
		e_media2=e_media2/pow(10,4);

		sigma_mag=mag2-pow(mag,2);
		sigma_e=e_media2/(2*N)-pow(e_media,2)*2*N;

		calor=(e_media2-pow(2*N*e_media,2))/(pow(N,2)*T);

		for(i=0;i<N;i++)
		{
			cor[i]=cor[i]/(pow(10,4)*pow(N,2));

		}

		//Imprimimos en los ficheros los datos obtenidos para estas magnitudes.


		fprintf(f2,"%i\t%lf\t%lf\t%lf\n",N,T,mag,sigma_mag);
		fprintf(f3,"%i\t%lf\t%lf\t%lf\n",N,T,e_media,sigma_e);
		fprintf(f4,"%lf\t%lf\n",T,calor);	
		
		for(i=0;i<N;i++)
		{
			fprintf(f5,"%i\t%lf\t%i\t%lf\n",N,T,i,cor[i]);

		}
	}

	fclose(f1);
	fclose(f2);
	fclose(f3);
	fclose(f4);
	fclose(f5);

	return 0;

}

//La función "conf_inicial_1" calcula la configuración inicial de la red con todos los espines orientados en la misma dirección.

void conf_inicial_1(int s[N][N])
{
	int i,j;
	
	for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				s[i][j]=1;		
			}
		}

	return;
}

//La función "conf_inicial_2" calcula la configuración inicial de la red de manera aleatoria 

void conf_inicial_2(int s[N][N])
{
	int i,j;
	extern gsl_rng *tau; 
	int semilla=892918;
		
	
	tau=gsl_rng_alloc(gsl_rng_taus); //Inicializamos el puntero
    gsl_rng_set(tau,semilla); //Inicializamos la semilla
	
	for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				s[i][j]=2*gsl_rng_uniform_int(tau,2)-1;		
			}
		}
	return;

}
//La siguiente función calcula el incremento de energía entre el estado inicial y el estado final teniendo en cuenta las condiciones de contorno.

double incremento_energia(int n, int m, int s[N][N])

{
	double delta_E;


	if(n==0 && m!=0 && m!=N-1)
	{
		delta_E=2*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][m+1]+s[n][m-1]);	
	}

	if(m==0 && n!=0 && n!=N-1)
	
	{
		delta_E=2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][m+1]+s[n][N-1]);	
	}
	
	if(n==N-1 && m!=0 && m!=N-1)
	{
		delta_E=2*s[n][m]*(s[0][m]+s[n-1][m]+s[n][m+1]+s[n][m-1]);	
	}

	if(m==N-1 && n!=0 && n!=N-1)
	{
		delta_E=2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][0]+s[n][m-1]);
	}

	if(n==0 && m==0)
	{
		delta_E=2*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][m+1]+s[n][N-1]);
	}

	if(n==0 && m==N-1)
	{
		delta_E=2*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][0]+s[n][m-1]);
	}

	if(n==N-1 && m==0)
	{
		delta_E=2*s[n][m]*(s[0][m]+s[n-1][m]+s[n][m+1]+s[n][N-1]);
	}

	if(n==N-1 && m==N-1)
	{
		delta_E=2*s[n][m]*(s[0][m]+s[n-1][m]+s[n][0]+s[n][m-1]);
	}
	
	if(n!=0 && n!=N-1 && m!=0 && m!=N-1)
	{
		delta_E=2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][m+1]+s[n][m-1]);
	}
	
	return delta_E;
}

//La función "magnet" calcula la magnetización del sistema. 

double magnet(int s[N][N])
{
	int i,j;	
	double m;
	double sum=0;
	
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			sum=sum+s[i][j];

		}
	}

	m=fabs(sum)/(pow(N,2));	
	
	return m;
}

//La función "energia" calcula la energía del sistema. 

double energia(int s[N][N])
{
	int i,j,k;
	double E=0;

	for(i=1;i<N-1;i++)
	{
		for(j=1;j<N-1;j++)
		{
			E=E+s[i][j]*(s[i][j+1]+s[i][j-1]+s[i+1][j]+s[i-1][j]);

		}
	}

	for(j=1;j<N-1;j++)
	{

		E=E+s[0][j]*(s[0][j+1]+s[0][j-1]+s[1][j]+s[N-1][j])+s[N-1][j]*(s[N-1][j+1]+s[N-1][j-1]+s[0][j]+s[N-2][j]);;

	}

	for(i=1;i<N-1;i++)
	{
		E=E+s[i][0]*(s[i][1]+s[i][N-1]+s[i+1][0]+s[i-1][0])+s[i][N-1]*(s[i][0]+s[i][N-2]+s[i+1][N-1]+s[i-1][N-1]);;
	
	}

	E=E+s[0][0]*(s[0][1]+s[0][N-1]+s[1][0]+s[N-1][0]);
	E=E+s[0][N-1]*(s[0][0]+s[0][N-2]+s[1][N-1]+s[N-1][N-1]);
	E=E+s[N-1][0]*(s[N-1][1]+s[N-1][N-1]+s[0][0]+s[N-2][0]);
	E=E+s[N-1][N-1]*(s[N-1][0]+s[N-1][N-2]+s[0][N-1]+s[N-2][N-1]);
	
	return -0.5*E;
}

//La función "correlacion" calcula la correlación entre los espines del sistema. 

double correlacion(int s[N][N],int i)
{
	int n,m;
	double sum=0;
	
	for(m=0;m<N;m++)
	{
		for(n=0;n<N-i;n++)
		{

			sum=sum+s[n][m]*s[n+i][m];

		}
	
	}
	
	return sum;



}


