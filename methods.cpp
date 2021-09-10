#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "methods.h"


void Sonda::escribirSondeo(std::ofstream& f)
{
	using namespace std;
	f << fixed << setprecision(3);
	f << "Tiempo\t U1\t U2\t U3\t rho\t v\t T\t P\t a\t H\t E\t M\t Para posición: x= \n" << punto_sondeo[0].x << "\n";
	for (int i{ 0 }; i<llamadas; i++) 
	{
		f << left << tiempos[i] << setw(2) << "\t ";
		f << left << punto_sondeo[i].U1 << setw(2) << "\t ";
		f << left << punto_sondeo[i].U2 << setw(2) << "\t ";
		f << left << punto_sondeo[i].U3 << setw(2) << "\t ";
		f << left << punto_sondeo[i].rho << setw(2) << "\t ";
		f << left << punto_sondeo[i].v << setw(2) << "\t ";
		f << left << punto_sondeo[i].T << setw(2) << "\t ";
		f << left << punto_sondeo[i].P << setw(2) << "\t ";
		f << left << punto_sondeo[i].a << setw(2) << "\t ";
		f << left << punto_sondeo[i].H << setw(2) << "\t ";
		f << left << punto_sondeo[i].E << setw(2) << "\t ";
		f << left << punto_sondeo[i].M << setw(2) << "\t ";
		f << "\n";
	};
}

inline void Frontera::actualizaFrontera(Nodo& n0, const Nodo& n1, const Nodo& n2, const double& del_t_x, const METHOD method) // Int lineal
{	
	static Gradiente dU{};
	static const double H_in{ n_ghost.H };
	// El flujo másico se fija con las condiciones en el punto sónico. Los factores 0.634 y 0.833 se fijan de la solución exacta
	static const double flumas{ 0.634*n_ghost.rho*s.Amin*sqrt(Gamma*R*0.833*n_ghost.T) };
	// n[i].U1 -= dU[i].p_U1*del_t_x;
	
	if (n0.entrada == true)    // SI ES ENTRADA
	{
		if (libre == false) // Si forzamos condiciones a la entrada
		{
			if ( n0.M >= 1.0 )	// Se fijan las 3 variables a la entrada
			{
				n0.U1 = n_ghost.U1;
				n0.U2 = flumas;
				n0.v = n0.U2/n0.U1;
				n0.U3 = n0.A*(n0.P/(Gamma-1) + 0.5*n0.rho*n0.v*n0.v);
			}
			else				// Se fijan 2 variables a la entrada, la otra fluctua (U2, U3 fijas)
			{
				//if (method != METHOD::ROE)	n0.U1 = 2*n1.U1 - n2.U1;
				if (method == METHOD::MACCORMACK)	n0.U1 = 2*n1.U1 - n2.U1;
				else if (method == METHOD::HLL)
				{
					SolverHLLS(n_ghost, n0, n1, dU);
					n0.U1 -= dU.p_U1*del_t_x;
					n_ghost = n0;
				}
				else if (method == METHOD::HLLC)
				{
					SolverHLLCS(n_ghost, n0, n1, dU);
					n0.U1 -= dU.p_U1*del_t_x;
					n_ghost = n0;
				}
				n0.U2 = flumas;

				n0.v = n0.U2 / n0.U1;
				n0.rho = n0.U1 / n0.A;
				n0.P = n0.rho*((Gamma-1)/Gamma)*(H_in - 0.5*n0.v*n0.v);//H_in
				n0.U3 = n0.A*(n0.P/(Gamma-1) + 0.5*n0.rho*n0.v*n0.v);
			}
		}
		else
		{
			if (method == METHOD::MACCORMACK)
			{
				n0.U1 = 2*n1.U1 - n2.U1;
				n0.U2 = 2*n1.U2 - n2.U2;
				n0.U3 = 2*n1.U3 - n2.U3;
			}
			else if (method == METHOD::HLL)
			{
				n_ghost = n0;
				SolverHLLS(n_ghost, n0, n1, dU);
				n0.U1 -= dU.p_U1*del_t_x;
				n0.U2 -= dU.p_U2*del_t_x;
				n0.U3 -= dU.p_U3*del_t_x;
			}
			else if (method == METHOD::HLLC)
			{
				n_ghost = n0;
				SolverHLLCS(n_ghost, n0, n1, dU);
				n0.U1 -= dU.p_U1*del_t_x;
				n0.U2 -= dU.p_U2*del_t_x;
				n0.U3 -= dU.p_U3*del_t_x;
			}
		}
	}
	else    // SI ES SALIDA
	{
		if (libre == false)	// Fijo una y dos fluctuan (salida forzada)
		{ 
			// if (method != METHOD::ROE)
			// {
			// 	n0.U1 = 2*n2.U1 - n1.U1;
			// 	n0.U2 = 2*n2.U2 - n1.U2;
			// }
			if (method == METHOD::MACCORMACK)
			{
				n0.U1 = 2*n2.U1 - n1.U1;
				n0.U2 = 2*n2.U2 - n1.U2;
			}
			else if (method == METHOD::HLL)
			{
				SolverHLLS(n2, n0, n_ghost, dU);
				n0.U1 -= dU.p_U1*del_t_x;
				n0.U2 -= dU.p_U2*del_t_x;
				n_ghost = n0;
			}
			else if (method == METHOD::HLLC)
			{
				SolverHLLCS(n2, n0, n_ghost, dU);
				n0.U1 -= dU.p_U1*del_t_x;
				n0.U2 -= dU.p_U2*del_t_x;
				n_ghost = n0;
			}
			n0.v = n0.U2 / n0.U1;
			n0.rho = n0.U1 / n0.A;
			n0.U3 = n0.U1 * ( s.P[s.N_puntos_definidos-1]/(n0.rho*(Gamma-1)) + 0.5*n0.v*n0.v);
		}
		else 				// Si es salida libre las 3 fluctuan
		{
			if (method != METHOD::ROE)
			{
				n0.U1 = 2*n2.U1 - n1.U1;
				n0.U2 = 2*n2.U2 - n1.U2;
				n0.U3 = 2*n2.U3 - n1.U3;
			}
		}
	}
}

void Malla::inicializar()
{
	int c{ 0 };
	double aux{};

	for(int i=0; i<s.N_nodos; i++)
	{
		n[i].A = areaTubo(n[i].x);
		
		//LINEAL A TRAMOS PARA INICIALIZAR
		aux = (s.nodos_definidos[c+1]-s.nodos_definidos[c])*delta_x;
		n[i].T = s.T[c] - (s.T[c] - s.T[c+1])/aux * ( n[i].x - n[s.nodos_definidos[c]].x );			
		n[i].v = s.v[c] - (s.v[c] - s.v[c+1])/aux * ( n[i].x - n[s.nodos_definidos[c]].x );
		n[i].P = s.P[c] - (s.P[c] - s.P[c+1])/aux * ( n[i].x - n[s.nodos_definidos[c]].x );
		//FERMI PARA INICIALIZAR
		//n[i].T = s.T[c+1] + (s.T[c]-s.T[c+1])/( 1 + 0.001*exp(8*n[i].x));
		//n[i].v = s.v[c+1] + (s.v[c]-s.v[c+1])/( 1 + 0.001*exp(8*n[i].x));
		//n[i].P = s.P[c+1] + (s.P[c]-s.P[c+1])/( 1 + 0.001*exp(8*n[i].x));

		n[i].rho = n[i].P/(n[i].T*R);
		n[i].H = (R/(Gamma-1))*n[i].T+n[i].P/n[i].rho+0.5*n[i].v*n[i].v;
		n[i].E = n[i].P/(Gamma-1) + 0.5*n[i].rho*n[i].v*n[i].v;
		n[i].a = sqrt(Gamma*n[i].P/n[i].rho);
		n[i].M = n[i].v/n[i].a;
					
		n[i].U1 = n[i].rho*n[i].A;
		n[i].U2 = n[i].rho*n[i].A*n[i].v;
		n[i].U3 = n[i].E*n[i].A;

		if ( i == s.nodos_definidos[c+1] ) 
		{
			c++;
		}
	}
		
	delta_t = buscaDelta_t();
}

inline void promediosRoe(const Nodo& n_L, const Nodo& n_R, double& v_prom, double& H_prom, double& a_prom)
{
	static double sqrtI, sqrtD;
	
	sqrtI = sqrt(n_L.rho*n_L.A);		sqrtD = sqrt(n_R.rho*n_R.A);

	v_prom = (n_L.v*sqrtI + n_R.v*sqrtD)/
				(sqrtI + sqrtD);

	sqrtI = sqrt(n_L.rho);		sqrtD = sqrt(n_R.rho);
	
	H_prom = (n_L.H*sqrtI + n_R.H*sqrtD)/
				(sqrtI + sqrtD);
	a_prom = sqrt((Gamma-1)*(H_prom-0.5*v_prom*v_prom));
}

inline void flujo(const Nodo& n, std::array<double, 3>& F)
{
	static double aux, aux1;
	aux = n.U2 / n.U1;
	aux1 = n.U2*aux;
	F[0] = n.U2;
	F[1] = aux1 + (Gamma-1) * (n.U3 - 0.5 * aux1);
	F[2] = Gamma * n.U3*aux - 0.5*(Gamma-1) * aux1*aux;
	//F = { n.U2,
	//	  aux1 + (Gamma-1) * (n.U3 - 0.5 * aux1),
	//	  Gamma * n.U3*aux - 0.5*(Gamma-1) * aux1*aux };
}

inline void flujoHLL(const Nodo& n_L, const Nodo& n_R, const double& v_prom, const double& a_prom, std::array<double, 3>& F)
{
	static std::array<double, 3> F_R;
	static double S1, S2;			// Signal speeds
	static double lambda1, lambda3;			// Autovalores Roe
	static double aux;

	// Autovalores Roe
	lambda1 = v_prom - a_prom; 
	lambda3 = v_prom + a_prom;

	// Saco velocidades máxima y mínima de señal
	aux = n_L.v - n_L.a;
	S1 = aux <= lambda1 ? aux :  lambda1;
	aux = n_R.v + n_R.a;
	S2 = aux >= lambda3 ? aux :  lambda3;
	
	if (S1 >= 0) 
	{
		flujo(n_L, F);
	}
	else if (S2 <= 0)
	{
		flujo(n_R, F);
	}
	else
	{
		flujo(n_R, F_R);
		flujo(n_L, F);

		aux = 1./(S2 - S1);
		F[0] = ( S2*F[0] - S1*F_R[0] + S1*S2*(n_R.U1 - n_L.U1) )*aux;
		F[1] = ( S2*F[1] - S1*F_R[1] + S1*S2*(n_R.U2 - n_L.U2) )*aux;
		F[2] = ( S2*F[2] - S1*F_R[2] + S1*S2*(n_R.U3 - n_L.U3) )*aux;
	}
}

inline void flujoHLLS(const Nodo& n_L, const Nodo& n_R, double& v_prom, double& H_prom, double& a_prom, std::array<double, 3>& F, const bool& minus)
{
	static std::array<double, 3> F_R{0, 0, 0};
	static std::array<double, 3> S1S2_H{0, 0, 0};
	static double S1, S2;					// Signal speeds
	static double deltaA, P_prom;
	static double aux, cte;

	//static int used{ 0 };
	//static double S1_old, S2_old;
	static double S1l, S1r, S2l, S2r;
	//static double N_usos{ static_cast<int>(!s.salida_libre)*5./(areaTubo(1+s.delta_x) - s.Amin) };

	if ( n_L.rho*n_R.rho < TOL )
	{
		v_prom = H_prom = a_prom = 0.0;
	}

	// Para modifiscación del término fuente
	deltaA = n_R.A - n_L.A; // Variación de área en i+1/2
	P_prom =  0.5*(n_R.P + n_L.P);
	cte = P_prom*deltaA;

	// Autovalores Roe como estimadores de velocidades máxima y mínima de señal
	S1 = v_prom - a_prom;
	S2 = v_prom + a_prom;

	S1l = n_L.v - n_L.a;
	S1r = n_R.v - n_R.a;
	S2l = n_L.v + n_L.a;
	S2r = n_R.v + n_R.a;

	if ( cte == 0 ) // Si no hay término fuente, utilizo siempre los estimadores de HLL
	{
		S1 = S1l < S1 ? S1l : S1;
		S1 = S1r < S1 ? S1r : S1;

		S2 = S2r > S2 ? S2r : S2;
		S2 = S2l > S2 ? S2l : S2;
	}
	else
	{
		// Si hay rarefacción transónica aplico estimadores HLL para evitar violación entropía.
		if ( S1l < 0 && S1r > 0 ) 
		{
			S1 = S1l < S1 ? S1l : S1;
			S1 = S1r < S1 ? S1r : S1;
			// if (S1_old == S1 && S2_old == S2)	used++;
			// S1_old = S1;
			// S2_old = S2;
		}
		else if ( S2l < 0 && S2r > 0 )
		{
			S2 = S2r > S2 ? S2r : S2;
			S2 = S2l > S2 ? S2l : S2;
			// if (S1_old == S1 && S2_old == S2)	used++;
			// S1_old = S1;
			// S2_old = S2;
		}
	}


	if (S1 >= 0) 
	{
		flujo(n_L, F);
		if (!minus)	F[1] += cte; //Resto vector términos fuente (solo tiene segunda componente)
	}
	else if (S2 <= 0)
	{
		flujo(n_R, F);
		if (minus)	F[1] += -cte; //Resto vector términos fuente
	}
	else
	{
		// Calculo S1*S2*H, donde H es el vector del artículo
		aux = (S1+S2)*(S1+S2);
		S1S2_H[0] = -cte*( Gamma );
		S1S2_H[2] = -cte*( Gamma*(3-Gamma)*aux*0.125 - S1*S2 )/(Gamma-1);

		flujo(n_R, F_R);
		flujo(n_L, F);

		aux = 1./(S2 - S1);
		F[0] = ( S2*F[0] - S1*F_R[0] + S1*S2*(n_R.U1 - n_L.U1) )*aux;
		F[1] = ( S2*F[1] - S1*F_R[1] + S1*S2*(n_R.U2 - n_L.U2) )*aux;
		F[2] = ( S2*F[2] - S1*F_R[2] + S1*S2*(n_R.U3 - n_L.U3) )*aux;
		
		F[0] += - S1S2_H[0]*aux;
		F[2] += - S1S2_H[2]*aux;
		if (minus)	F[1] += S1*cte*aux;
		else		F[1] += S2*cte*aux;
	}
}

inline void flujoHLLC(const Nodo& n_L, const Nodo& n_R, const double& v_prom, const double& a_prom, std::array<double, 3>& F)
{
	static double S1, S2, S_star;			// Signal speeds
	static double S1l, S1r, S2l, S2r;			
	static double aux;

	// Autovalores Roe como estimadores de velocidades máxima y mínima de señal
	S1 = v_prom - a_prom;
	S2 = v_prom + a_prom;

	S1l = n_L.v - n_L.a;
	S1r = n_R.v - n_R.a;
	S2l = n_L.v + n_L.a;
	S2r = n_R.v + n_R.a;

	// Saco velocidades máxima y mínima de señal
	S1 = S1l < S1 ? S1l : S1;
	S1 = S1r < S1 ? S1r : S1;

	S2 = S2r > S2 ? S2r : S2;
	S2 = S2l > S2 ? S2l : S2;
	
	if (S1 >= 0) 
	{
		flujo(n_L, F);
	}
	else if (S2 <= 0)
	{
		flujo(n_R, F);
	}
	else
	{
		S_star = ( n_R.P - n_L.P + n_L.U2*(S1 - n_L.v) - n_R.U2*(S2 - n_R.v) )/( n_L.U1*(S1 - n_L.v) - n_R.U1*(S2 - n_R.v) );
		
		if (S_star >= 0)
		{
			flujo(n_L, F);
			
			aux = n_L.U1*(S1 - n_L.v)/(S1 - S_star);
			F[0] += S1 * ( aux - n_L.U1);
			F[1] += S1 * ( aux*S_star - n_L.U2);
			F[2] += S1 * ( aux*( n_L.E/n_L.rho + (S_star - n_L.v)*( S_star + n_L.P/( n_L.rho*(S1 - n_L.v) ) ) ) - n_L.U3);
		}
		else
		{
			flujo(n_R, F);

			aux = n_R.U1*(S2 - n_R.v)/(S2 - S_star);
			F[0] += S2 * ( aux - n_R.U1);
			F[1] += S2 * ( aux*S_star - n_R.U2);
			F[2] += S2 * ( aux*( n_R.E/n_R.rho + (S_star - n_R.v)*( S_star + n_R.P/( n_R.rho*(S2 - n_R.v) ) ) ) - n_R.U3);
		}
	}
}

inline void flujoHLLCS(const Nodo& n_L, const Nodo& n_R, const double& v_prom, const double& H_prom, const double& a_prom, std::array<double, 3>& F, const bool& minus)
{
	static std::array<double, 3> S1S2_H;
	static double S1, S2, S_star, S_star_num, S_star_den;					// Signal speeds
	static double S1_bis, S2_bis;
	static double deltaA, P_prom;
	static double aux, cte;

	//static int used{ 0 };
	//static double S1_old, S2_old;
	static double S1l, S1r, S2l, S2r;
	//static double N_usos{ 10./(areaTubo(1+s.delta_x) - s.Amin) };
	

	deltaA = n_R.A - n_L.A; // Variación de área en i+1/2
	P_prom = 0.5 * (n_R.P + n_L.P);
	cte = P_prom*deltaA;

	// Autovalores Roe como estimadores de velocidades máxima y mínima de señal
	S1 = v_prom - a_prom;
	S2 = v_prom + a_prom;

	S1l = n_L.v - n_L.a;
	S1r = n_R.v - n_R.a;
	S2l = n_L.v + n_L.a;
	S2r = n_R.v + n_R.a;

	if ( cte == 0 ) // Si no hay término fuente, utilizo siempre los estimadores de HLL
	{
		S1 = S1l < S1 ? S1l : S1;
		S1 = S1r < S1 ? S1r : S1;

		S2 = S2r > S2 ? S2r : S2;
		S2 = S2l > S2 ? S2l : S2;
	}
	else
	{
		// Si hay rarefacción transónica aplico estimadores HLL para evitar violación de la entropía.
		// if ( used < N_usos )
		// {
		if ( S1l < 0 && S1r > 0 ) 
		{
			S1 = S1l < S1 ? S1l : S1;
			S1 = S1r < S1 ? S1r : S1;
			// if (S1_old == S1 && S2_old == S2)	used++;
			// S1_old = S1;
			// S2_old = S2;
		}
		else if ( S2l < 0 && S2r > 0 )
		{
			S2 = S2r > S2 ? S2r : S2;
			S2 = S2l > S2 ? S2l : S2;
			// if (S1_old == S1 && S2_old == S2)	used++;
			// S1_old = S1;
			// S2_old = S2;
		}
		//}
	}


	if (S1 >= 0) 
	{
		flujo(n_L, F);
		if (!minus)
		{
			F[1] += cte;
		}
	}
	else if (S2 <= 0)
	{
		flujo(n_R, F);
		if (minus)
		{
			F[1] += -cte;
		}
	}
	else
	{
		// Me construyo el vector H_L * S1*S2
		S1_bis = n_L.v - n_L.a;
		S2_bis = n_L.v + n_L.a;
		aux = (S1_bis+S2_bis)*(S1_bis+S2_bis);
		S1S2_H[0] = -cte*( Gamma );
		S1S2_H[2] = -cte*( Gamma*(3-Gamma)*aux - 8*S1_bis*S2_bis )/(Gamma-1)*0.125;


		// Calculo S_star+ y veo si es mayor que 0
		S_star_num = ( n_R.P - n_L.P + S2*n_L.U1*(S1 - n_L.v) - S1*n_R.U1*(S2 - n_R.v) );
		S_star_den = n_L.U1*(S1 - n_L.v) - n_R.U1*(S2 - n_R.v);
		S_star = S2*(S_star_num + S1S2_H[0]) / ( S2*S_star_den + S1S2_H[0] ); // Multiplico todo por S2 por si fuese 0 que no divida por 0

		if (S_star >= 0)
		{
			
			flujo(n_L, F);
			if (!minus)
			{
				F[1] += P_prom * deltaA;
			}

			aux = (n_L.U1*(S1 - n_L.v) + S1S2_H[0]/S2)/(S1 - S_star);
			F[0] += S1*aux - S1S2_H[0]/S2 - S1*n_L.U1;
			F[1] += S1*aux*S_star - S1*n_L.U2;
			F[2] += S1*aux*( n_L.E/n_L.rho + (S_star - n_L.v)*( S_star + n_L.P/( n_L.rho*(S1 - n_L.v) ) ) ) - S1S2_H[2]/S2 - S1*n_L.U3;
		}
		else
		{
			// Me construyo el vector H_R
			S1_bis = n_R.v - n_R.a;
			S2_bis = n_R.v + n_R.a;
			aux = (S1_bis+S2_bis)*(S1_bis+S2_bis);
			S1S2_H[0] = -cte*( Gamma );
			S1S2_H[2] = -cte*( Gamma*(3-Gamma)*aux - 8*S1_bis*S2_bis )/(Gamma-1)*0.125;

			// S_star-
			S_star = S1*(S_star_num + S1S2_H[0]) / ( S1*S_star_den + S1S2_H[0] ); //Multiplico todo por S1 por si es 0, no dividir para 0
			
			flujo(n_R, F);
			if (minus)
			{
				F[1] += -P_prom * deltaA;
			}

			aux = (n_R.U1*(S2 - n_R.v) - S1S2_H[0]/S1)/(S2 - S_star);
			F[0] += S2*aux + S1S2_H[0]/S1 - S2*n_R.U1;
			F[1] += S2*aux*S_star - S2*n_R.U2;
			F[2] += S2*aux*( n_R.E/n_R.rho + (S_star - n_R.v)*( S_star + n_R.P/( n_R.rho*(S2 - n_R.v) ) ) ) + S1S2_H[2]/S1 - S2*n_R.U3;
		}
	}
}

inline void SolverMacCormack(const Nodo& n_L, const Nodo& n_R, //i, i+1
					Gradiente& dU, const bool& pred=false)
{
	static double deltaA{};
	static std::array<double, 3> F_L, F_R;

	deltaA = n_R.A - n_L.A; // Variación de área
	flujo(n_L, F_L);
	flujo(n_R, F_R);

	if ( pred == true )	// Predictor
	{
		dU = { F_R[0] - F_L[0],
			   F_R[1] - F_L[1] - n_L.P * deltaA,	// Término fuente del área
			   F_R[2] - F_L[2] 
			  };
	}
	else 				// Corrector
	{
		dU = { F_R[0] - F_L[0],
			   F_R[1] - F_L[1] - n_R.P * deltaA,
			   F_R[2] - F_L[2]
			  };
	}
}

inline void SolverRoeRaw(const Nodo& n_L, const Nodo& n_R, //i, i+1 Sin la corrección de la entropía.
					Gradiente& dU_L, Gradiente& dU_R)
{
	static double lambda1, lambda2, lambda3;	// Autovalores
	static std::array<double, 3> e1, e2, e3;	// Autovectores
	static double alpha1, alpha2, alpha3;		// Wave-strenghts
	static double beta1, beta2, beta3;			// Coeficientes del término fuente
	static double v_prom, H_prom, a_prom, P_prom, deltaA;
	static double deltaU1, deltaU2, deltaU3;
	// Variables auxiliares
	static double aux;

	if ( n_L.rho*n_R.rho < TOL )
	{
		v_prom = H_prom = a_prom = 0.0;
	}
	else
	{
		promediosRoe(n_L, n_R, v_prom, H_prom, a_prom);
	}
	
	// Autovalores
	lambda1 = v_prom - a_prom;
	lambda2 = v_prom;
	lambda3 = v_prom + a_prom;

	// Autovectores
	e1[0] = 1.;      e1[1] = v_prom - a_prom;    e1[2] = H_prom - v_prom*a_prom;
	e2[0] = 1.;      e2[1] = v_prom;    		 e2[2] = 0.5*v_prom*v_prom;
	e3[0] = 1.;      e3[1] = v_prom + a_prom;    e3[2] = H_prom + v_prom*a_prom;
	
	deltaU1 = n_R.U1 - n_L.U1;
	deltaU2 = n_R.U2 - n_L.U2;
	deltaU3 = n_R.U3 - n_L.U3;

	//Coeficientes alfa
	if ( n_L.rho*n_R.rho < TOL )
	{
		alpha1 = alpha2 = alpha3 = 0.0;
	}
	else
	{
		aux = 2*H_prom - v_prom*v_prom;
		alpha1 = 1./(a_prom*aux) * (
				(H_prom*v_prom + 0.5*a_prom*v_prom*v_prom - 0.5*v_prom*v_prom*v_prom) * deltaU1 +
				(0.5*v_prom*v_prom - H_prom - a_prom*v_prom) * deltaU2 +
				a_prom * deltaU3  );
		alpha2 = 2./aux * (
				(H_prom - v_prom*v_prom) * deltaU1 + 
				v_prom * deltaU2 - 
				deltaU3  );
		alpha3 = 1./(a_prom*aux) * (
				(-H_prom*v_prom + 0.5*a_prom*v_prom*v_prom + 0.5*v_prom*v_prom*v_prom) * deltaU1 +
				(-0.5*v_prom*v_prom + H_prom - a_prom*v_prom) * deltaU2 +	
				a_prom * deltaU3  );
	}

	deltaA = n_R.A - n_L.A; //Variación de área
	P_prom = 0.5 * (n_R.P + n_L.P);
			
	//Coeficientes beta para la descomposición del término fuente S
	aux = (P_prom*deltaA/a_prom) / (2*H_prom-v_prom*v_prom);
	beta1 = -aux * (H_prom + v_prom*a_prom - 0.5*v_prom*v_prom);
	beta2 = aux * (2*v_prom*a_prom);
	beta3 = aux * (H_prom - v_prom*a_prom - 0.5*v_prom*v_prom);
	

	if (lambda1 < 0)
	{
		dU_L = { dU_L.p_U1 + (lambda1*alpha1 - beta1)*e1[0],
					dU_L.p_U2 + (lambda1*alpha1 - beta1)*e1[1],
					dU_L.p_U3 + (lambda1*alpha1 - beta1)*e1[2]
				};
	}
	else
	{
		dU_R = { dU_R.p_U1 + (lambda1*alpha1 - beta1)*e1[0],
					dU_R.p_U2 + (lambda1*alpha1 - beta1)*e1[1],
					dU_R.p_U3 + (lambda1*alpha1 - beta1)*e1[2]
				};
		}
	
	if (lambda2 < 0)
	{
		dU_L = { dU_L.p_U1 + (lambda2*alpha2 - beta2)*e2[0],
				 dU_L.p_U2 + (lambda2*alpha2 - beta2)*e2[1],
				 dU_L.p_U3 + (lambda2*alpha2 - beta2)*e2[2]
			   };
	}
	else
	{
		dU_R = { dU_R.p_U1 + (lambda2*alpha2 - beta2)*e2[0],
				 dU_R.p_U2 + (lambda2*alpha2 - beta2)*e2[1],
				 dU_R.p_U3 + (lambda2*alpha2 - beta2)*e2[2]
			   };
	}

	if (lambda3 < 0)
	{
		dU_L = { dU_L.p_U1 + (lambda3*alpha3 - beta3)*e3[0],
					dU_L.p_U2 + (lambda3*alpha3 - beta3)*e3[1],
					dU_L.p_U3 + (lambda3*alpha3 - beta3)*e3[2]
				};
	}
	else
	{
		dU_R = { dU_R.p_U1 + (lambda3*alpha3 - beta3)*e3[0],
					dU_R.p_U2 + (lambda3*alpha3 - beta3)*e3[1],
					dU_R.p_U3 + (lambda3*alpha3 - beta3)*e3[2]
				};
	}
}

inline void SolverRoe(const Nodo& n_L, const Nodo& n_R, //i, i+1
					Gradiente& dU_L, Gradiente& dU_R)
{
	static double lambda1, lambda2, lambda3;	// Autovalores
	static std::array<double, 3> e1, e2, e3;	// Autovectores
	static double alpha1, alpha2, alpha3;		// Wave-strenghts
	static double beta1, beta2, beta3;			// Coeficientes del término fuente
	static double lb1l, lb1r, lb3l, lb3r;
	static double v_prom, H_prom, a_prom, P_prom, deltaA;
	static double deltaU1, deltaU2, deltaU3;
	// Variables auxiliares
	static double aux;

	if ( n_L.rho*n_R.rho < TOL )
	{
		v_prom = H_prom = a_prom = 0.0;
	}
	else
	{
		promediosRoe(n_L, n_R, v_prom, H_prom, a_prom);
	}
	
	// Autovalores
	lb1l = n_L.v - n_L.a;
	lb1r = n_R.v - n_R.a;
	lambda1 = v_prom - a_prom;

	lambda2 = v_prom;

	lb3l = n_L.v + n_L.a;
	lb3r = n_R.v + n_R.a;
	lambda3 = v_prom + a_prom;

	// Autovectores
	e1[0] = 1.;      e1[1] = v_prom - a_prom;    e1[2] = H_prom - v_prom*a_prom;
	e2[0] = 1.;      e2[1] = v_prom;    		 e2[2] = 0.5*v_prom*v_prom;
	e3[0] = 1.;      e3[1] = v_prom + a_prom;    e3[2] = H_prom + v_prom*a_prom;
	
	deltaU1 = n_R.U1 - n_L.U1;
	deltaU2 = n_R.U2 - n_L.U2;
	deltaU3 = n_R.U3 - n_L.U3;

	//Coeficientes alfa
	if ( n_L.rho*n_R.rho < TOL )
	{
		alpha1 = alpha2 = alpha3 = 0.0;
	}
	else
	{
		aux = 2*H_prom - v_prom*v_prom;
		alpha1 = 1./(a_prom*aux) * (
				(H_prom*v_prom + 0.5*a_prom*v_prom*v_prom - 0.5*v_prom*v_prom*v_prom) * deltaU1 +
				(0.5*v_prom*v_prom - H_prom - a_prom*v_prom) * deltaU2 +
				a_prom * deltaU3  );
		alpha2 = 2./aux * (
				(H_prom - v_prom*v_prom) * deltaU1 + 
				v_prom * deltaU2 - 
				deltaU3  );
		alpha3 = 1./(a_prom*aux) * (
				(-H_prom*v_prom + 0.5*a_prom*v_prom*v_prom + 0.5*v_prom*v_prom*v_prom) * deltaU1 +
				(-0.5*v_prom*v_prom + H_prom - a_prom*v_prom) * deltaU2 +	
				a_prom * deltaU3  );
	}

	deltaA = n_R.A - n_L.A; //Variación de área
	P_prom = 0.5 * (n_R.P + n_L.P);
			
	//Coeficientes beta para la descomposición del término fuente S
	aux = (P_prom*deltaA/a_prom) / (2*H_prom-v_prom*v_prom);
	beta1 = -aux * (H_prom + v_prom*a_prom - 0.5*v_prom*v_prom);
	beta2 = aux * (2*v_prom*a_prom);
	beta3 = aux * (H_prom - v_prom*a_prom - 0.5*v_prom*v_prom);
	
	if (lb1l < 0 && lb1r > 0)
	{ 
		aux = lb1l*(lb1r-lambda1)/(lb1r-lb1l); //TÉRMINO CORRECCIÓN ENTROPÍA
		dU_L = { dU_L.p_U1 + (aux*alpha1 - beta1)*e1[0],
				 dU_L.p_U2 + (aux*alpha1 - beta1)*e1[1],
				 dU_L.p_U3 + (aux*alpha1 - beta1)*e1[2]
			   };

		aux = lb1r*(lambda1-lb1l)/(lb1r-lb1l);
		dU_R = { dU_R.p_U1 + (aux*alpha1)*e1[0],
				 dU_R.p_U2 + (aux*alpha1)*e1[1],
				 dU_R.p_U3 + (aux*alpha1)*e1[2]
			   };
	}
	else
	{
		if (lambda1 < 0)
		{
			dU_L = { dU_L.p_U1 + (lambda1*alpha1 - beta1)*e1[0],
					 dU_L.p_U2 + (lambda1*alpha1 - beta1)*e1[1],
					 dU_L.p_U3 + (lambda1*alpha1 - beta1)*e1[2]
				   };
		}
		else
		{
			dU_R = { dU_R.p_U1 + (lambda1*alpha1 - beta1)*e1[0],
					 dU_R.p_U2 + (lambda1*alpha1 - beta1)*e1[1],
					 dU_R.p_U3 + (lambda1*alpha1 - beta1)*e1[2]
				   };
		}
	}
	
	if (lambda2 < 0)
	{
		dU_L = { dU_L.p_U1 + (lambda2*alpha2 - beta2)*e2[0],
				 dU_L.p_U2 + (lambda2*alpha2 - beta2)*e2[1],
				 dU_L.p_U3 + (lambda2*alpha2 - beta2)*e2[2]
			   };
	}
	else
	{
		dU_R = { dU_R.p_U1 + (lambda2*alpha2 - beta2)*e2[0],
				 dU_R.p_U2 + (lambda2*alpha2 - beta2)*e2[1],
				 dU_R.p_U3 + (lambda2*alpha2 - beta2)*e2[2]
			   };
	}

	if (lb3l < 0 && lb3r > 0)
	{
		aux = lb3l*(lb3r-lambda3)/(lb3r-lb3l);
		dU_L = { dU_L.p_U1 + (aux*alpha3)*e3[0],
				 dU_L.p_U2 + (aux*alpha3)*e3[1],
				 dU_L.p_U3 + (aux*alpha3)*e3[2]
			   };
		
		aux = lb3r*(lambda3-lb3l)/(lb3r-lb3l);
		dU_R = { dU_R.p_U1 + (aux*alpha3 - beta3)*e3[0],
				 dU_R.p_U2 + (aux*alpha3 - beta3)*e3[1],
				 dU_R.p_U3 + (aux*alpha3 - beta3)*e3[2]
			   };
	}
	else
	{
		if (lambda3 < 0)
		{
			dU_L = { dU_L.p_U1 + (lambda3*alpha3 - beta3)*e3[0],
					 dU_L.p_U2 + (lambda3*alpha3 - beta3)*e3[1],
					 dU_L.p_U3 + (lambda3*alpha3 - beta3)*e3[2]
				   };
		}
		else
		{
			dU_R = { dU_R.p_U1 + (lambda3*alpha3 - beta3)*e3[0],
					 dU_R.p_U2 + (lambda3*alpha3 - beta3)*e3[1],
					 dU_R.p_U3 + (lambda3*alpha3 - beta3)*e3[2]
				   };
		}
	}
}

inline void SolverHLL(const Nodo& n_L, const Nodo& n_C, const Nodo& n_R, //i, i+1
    	                Gradiente& dU)
{
	static int calls{ 0 };
	static std::array<double, 3> F_L, F_R;	// Flujos
	static const double CteA{ 4.*(s.Amax-s.Amin)/(s.L*s.L) };
	static double deltaA;
	// Guardo todos, consume más RAM pero evito calcularlos dos veces (en una llamada y la siguiente)
	static std::vector<double> v_prom( s.N_nodos+1 ), H_prom( s.N_nodos+1 ), a_prom( s.N_nodos+1 );

	// Calculo los promedios de Roe como estimadores para velocidad señales HLLC
	if (calls == 0)  promediosRoe(n_L, n_C, v_prom[calls], H_prom[calls], a_prom[calls]);
	flujoHLL(n_L, n_C, v_prom[calls], a_prom[calls], F_L);  //F_i-1/2
	//Sumo una llamada, calculo los otros promedios y flujo
	calls++;
	promediosRoe(n_C, n_R, v_prom[calls], H_prom[calls], a_prom[calls]);
	flujoHLL(n_C, n_R, v_prom[calls], a_prom[calls], F_R); //F_i+1/2

	deltaA =  ( CteA * (2*n_C.x - s.L) ) * s.delta_x;

	dU = { ( F_R[0] - F_L[0] ),
		   ( F_R[1] - F_L[1] ) - n_C.P * deltaA,
		   ( F_R[2] - F_L[2] )
		 };
	if (calls == s.N_nodos) calls = 0;
}

inline void SolverHLLS(const Nodo& n_L, const Nodo& n_C, const Nodo& n_R, //i, i+1
    	                Gradiente& dU)
{
	static int calls{ 0 };
	static std::array<double, 3> F_L{}, F_R{};	// Flujos
	// Guardo todos, consume más RAM pero evito calcularlos dos veces (en una llamada y la siguiente)
	static std::vector<double> v_prom( s.N_nodos+1 ), H_prom( s.N_nodos+1 ), a_prom( s.N_nodos+1 );

	// Calculo los promedios de Roe como estimadores para velocidad señales HLLC
	if (calls == 0) promediosRoe(n_L, n_C, v_prom[calls], H_prom[calls], a_prom[calls]);
	flujoHLLS(n_L, n_C, v_prom[calls], H_prom[calls], a_prom[calls], F_L, false);  //F_i-1/2
	//Sumo una llamada, calculo los otros promedios y flujo
	calls++;
	promediosRoe(n_C, n_R, v_prom[calls], H_prom[calls], a_prom[calls]);
	flujoHLLS(n_C, n_R, v_prom[calls], H_prom[calls], a_prom[calls], F_R, true); //F_i+1/2

	dU = { F_R[0] - F_L[0],
		   F_R[1] - F_L[1],
		   F_R[2] - F_L[2]
		 };
	if (calls == s.N_nodos) calls = 0;
}

inline void SolverHLLC(const Nodo& n_L, const Nodo& n_C, const Nodo& n_R, //i, i+1
    	                Gradiente& dU)
{
	static int calls{ 0 };
	static std::array<double, 3> F_L, F_R;	// Flujos
	static const double CteA{ 4.*(s.Amax-s.Amin)/(s.L*s.L) };
	static double deltaA;
	// Guardo todos, consume más RAM pero evito calcularlos dos veces (en una llamada y la siguiente)
	static std::vector<double> v_prom( s.N_nodos+1 ), H_prom( s.N_nodos+1 ), a_prom( s.N_nodos+1 );

	// Calculo los promedios de Roe como estimadores para velocidad señales HLLC
	if (calls == 0) promediosRoe(n_L, n_C, v_prom[calls], H_prom[calls], a_prom[calls]);
	flujoHLLC(n_L, n_C, v_prom[calls], a_prom[calls], F_L);
	//Sumo una llamada, calculo los otros promedios y flujo
	calls++;
	promediosRoe(n_C, n_R, v_prom[calls], H_prom[calls], a_prom[calls]); 
	flujoHLLC(n_C, n_R, v_prom[calls], a_prom[calls], F_R);

	deltaA = n_R.A - n_L.A;//( CteA * (2*n_C.x - s.L) ) * s.delta_x; //Variación de área *** INFLUYE EN RESULTADO FINAL LA FORMA EN LA QUE EXPRESAR ESTO ***

	dU = { ( F_R[0] - F_L[0] ),
		   ( F_R[1] - F_L[1] ) - n_C.P * deltaA, //n_C.rho*R*n_C.T * deltaA,
		   ( F_R[2] - F_L[2] )
		 };
	if (calls == s.N_nodos) calls = 0;
}

inline void SolverHLLCS(const Nodo& n_L, const Nodo& n_C, const Nodo& n_R, //i, i+1
    	                Gradiente& dU)
{
	static int calls{ 0 };
	static std::array<double, 3> F_L, F_R;	// Flujos
	// Guardo todos, consume más RAM pero evito calcularlos dos veces (en una llamada y la siguiente)
	static std::vector<double> v_prom( s.N_nodos+1 ), H_prom( s.N_nodos+1 ), a_prom( s.N_nodos+1 );

	// Calculo los promedios de Roe como estimadores para velocidad señales HLLC
	if (calls == 0) promediosRoe(n_L, n_C, v_prom[calls], H_prom[calls], a_prom[calls]);
	flujoHLLCS(n_L, n_C, v_prom[calls], H_prom[calls], a_prom[calls], F_L, false);
	//Sumo una llamada, calculo los otros promedios y flujo
	calls++;
	promediosRoe(n_C, n_R, v_prom[calls], H_prom[calls], a_prom[calls]); 
	flujoHLLCS(n_C, n_R, v_prom[calls], H_prom[calls], a_prom[calls], F_R, true);

	dU = { F_R[0] - F_L[0],
		   F_R[1] - F_L[1],
		   F_R[2] - F_L[2]
		 };
	if (calls == s.N_nodos) calls = 0;
}

void Malla::calcularNodos(Frontera& in_flow, Frontera& out_flow, const METHOD method)
{
	static std::vector<Gradiente> dU( s.N_nodos );
	static double del_t_x{};
	del_t_x = delta_t/delta_x;

	for (int i{ 0 }; i < s.N_nodos; i++)
	{
		dU[i] = {0, 0, 0};
	}

	switch (method) // Se elige qué solver utilizar
	{
	case METHOD::MACCORMACK:
		
		static std::vector<Nodo> n_pred( s.N_nodos );
		static std::vector<Gradiente> dU_old( s.N_nodos );
		
		n_pred = n;
		
		for (int i{ 0 }; i < s.N_nodos-1; i++) 			// Paso predictor
		{ 
		    SolverMacCormack(n_pred[i], n_pred[i+1],    // Nodos de las celdas i, i+1
		           				 dU[i], true);    		// Incrementos de U en i, indico que es paso predictor
		}
		dU_old = dU;

		for (int i{ 0 }; i < s.N_nodos-1; i++)
		{
		    n_pred[i].U1 -= dU[i].p_U1*del_t_x;
		    n_pred[i].U2 -= dU[i].p_U2*del_t_x;
		    n_pred[i].U3 -= dU[i].p_U3*del_t_x;
		}
		in_flow.actualizaFrontera(n_pred[0], n_pred[0 + 1], n_pred[0 + 2], del_t_x, method);
		out_flow.actualizaFrontera(n_pred[s.N_nodos-1], n_pred[s.N_nodos-1 - 2], n_pred[s.N_nodos-1 - 1], del_t_x, method);

		for (int i{ 0 }; i < s.N_nodos-1; i++) 			// Paso corrector
		{   
		    SolverMacCormack(n_pred[i], n_pred[i+1],
		           				 dU[i+1]);
		}
		
		for (int i{ 0 }; i < s.N_nodos-1; i++)
		{
		    dU[i] = {0.5 * (dU[i].p_U1 + dU_old[i].p_U1), 
		             0.5 * (dU[i].p_U2 + dU_old[i].p_U2),
		             0.5 * (dU[i].p_U3 + dU_old[i].p_U3)};
		}

		break;

	case METHOD::ROE:

		for (int i{ 0 }; i < s.N_nodos-1; i++)
		{
			SolverRoe(n[i], n[i+1],   			// Nodos de las celdas i, i+1
					 dU[i], dU[i+1]); 			// Incrementos de U en i, i+1
		}
		
		break;
	
	case METHOD::HLL:
		// PARECE MÁS ESTABLE CON EL HLL QUE CON EL MODIFICADO
		for (int i{ 1 }; i < s.N_nodos-1; i++)
        { 
            SolverHLLS(n[i-1], n[i], n[i+1],	// Nodos de las celdas i, i+1
                   	    	  dU[i]);    		// Incremento de U en i
        }
		
		break;

	case METHOD::HLLC:

		for (int i{ 1 }; i < s.N_nodos-1; i++)
        { 
            SolverHLLCS(n[i-1], n[i], n[i+1],	// Nodos de las celdas i, i+1
                   	    	   dU[i]);    		// Incremento de U en i
        }
		
		break;
	
	default:
		std::cout << "\n\nNO METHOD SELECTED. ABORTING EXECUTION.\n\n";
		exit(1);
	}

	// Ya he calculado todos los incrementos de U -> Avanzo el paso y saco valores finales.
	for (int i{ 0 }; i < s.N_nodos; i++)
	{
		n[i].U1 -= dU[i].p_U1*del_t_x;	// Signo menos porque dU/dt = -dF/dx (dentro de F van tb términos fuente)
		n[i].U2 -= dU[i].p_U2*del_t_x;
		n[i].U3 -= dU[i].p_U3*del_t_x;
		//if (PAUSA) std::cout << "\ni = " << i << ": " << dU[i].p_U1 << " " << dU[i].p_U2 << " " << dU[i].p_U3;
	}
	
	in_flow.actualizaFrontera(n[0], n[0 + 1], n[0 + 2], del_t_x, method);
	out_flow.actualizaFrontera(n[s.N_nodos-1], n[s.N_nodos-1 - 2], n[s.N_nodos-1 - 1], del_t_x, method);
	
	for (int i{ 0 }; i < s.N_nodos; i++)
	{
		n[i].rho = n[i].U1/n[i].A;
		n[i].v = n[i].U2/n[i].U1;
		n[i].E = n[i].U3/n[i].A;
		n[i].P = (Gamma-1)*(n[i].E - 0.5*n[i].rho*n[i].v*n[i].v);
		n[i].T = n[i].P/(n[i].rho*R);
		n[i].H = (n[i].E + n[i].P)/n[i].rho;
		n[i].a = sqrt(Gamma*R*n[i].T);
		n[i].M = n[i].v/n[i].a;
	}
	// Busco el paso de tiempo que usar y actualizo el tiempo
	tiempo += delta_t;
	delta_t = buscaDelta_t();
}

void Malla::escribirDatos(std::ofstream& f, const int precision)
{
	using namespace std;
	static bool primera_llamada{ false };

	if (primera_llamada == false) 
	{
		f << fixed << setprecision(precision);
		f << "Tiempo\t \n";
		f << "x\t U1\t U2\t U3\t rho\t v\t T\t P\t a\t A\t H\t E\t M\t \n";
	};
	primera_llamada = true;
	
	f << tiempo << "\n";
	
	for (const Nodo& nodo : n)
	{
		f << left << nodo.x << setw(2) << "\t ";
		f << left << nodo.U1 << setw(2) << "\t ";
		f << left << nodo.U2 << setw(2) << "\t ";
		f << left << nodo.U3 << setw(2) << "\t ";
		f << left << nodo.rho << setw(2) << "\t ";
		f << left << nodo.v << setw(2) << "\t ";
		f << left << nodo.T << setw(2) << "\t ";
		f << left << nodo.P << setw(2) << "\t ";
		f << left << nodo.a << setw(2) << "\t ";
		f << left << nodo.A << setw(2) << "\t ";
		f << left << nodo.H << setw(2) << "\t ";
		f << left << nodo.E << setw(2) << "\t ";
		f << left << nodo.M << setw(2) << "\t ";
		f << "\n";
	};
	f << "\n";
}

void Malla::sondear(Sonda& sonda, const int numero_nodo)
{
	sonda.punto_sondeo.resize(sonda.llamadas+1); 		// Añado una componente a los vectores y anoto datos
	sonda.tiempos.resize(sonda.llamadas+1);

	sonda.punto_sondeo[sonda.llamadas] = n[numero_nodo]; //	Copia el nodo
	sonda.tiempos[sonda.llamadas] = tiempo;              //	Copia el tiempo
	sonda.llamadas++;                                    //	Pasa al siguiente instante de sondeo
}