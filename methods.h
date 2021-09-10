#ifndef METHODS_H
#define METHODS_H

#include "constants.h"
#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>

extern double areaTubo(const double x);
extern bool PAUSA;

enum class METHOD
{
    MACCORMACK,
    ROE,
    HLL,
    HLLC,
};

class Settings
{
public:

    friend class Malla;
    friend class frontera;
    friend class sonda;

    double CFL;
    double delta_x;
    int N_nodos;
    double simTime;
    int intervalo_extraer;
    bool is_sonda;
    int N_nodo_sonda;

    bool entrada_libre;
    bool salida_libre;
    int N_puntos_definidos;
    std::vector<int> nodos_definidos;
    std::vector<double> P;
    std::vector<double> T;
    std::vector<double> v;
    std::vector<double> rho;

    double Amax;
    double Amin;
    double L;

    METHOD method = METHOD::ROE;

    Settings()
    {
        int aux{};
        double aux2{};
        FILE *setup;	// Archivo de configuración
        setup = fopen("inputs.in", "r");
        //LECTURA DE DATOS DE ENTRADA
        fscanf(setup,"%*s");
        fscanf(setup,"%*s %*s %lf", &CFL);
        fscanf(setup,"%*s %*s %lf", &delta_x);
        fscanf(setup,"%*s %*s %d", &N_nodos);
        fscanf(setup,"%*s %*s %lf", &simTime);
        fscanf(setup,"%*s %*s %d", &intervalo_extraer);
        fscanf(setup,"%*s %*s %d", &aux);
        is_sonda = static_cast<bool>(aux);
        fscanf(setup,"%*s %*s %d", &N_nodo_sonda);
        fscanf(setup,"%*s");
        fscanf(setup,"%*s %*s %d", &aux);
        entrada_libre = static_cast<bool>(aux);
        fscanf(setup,"%*s %*s %d", &aux);
        salida_libre = static_cast<bool>(aux);
        fscanf(setup,"%*s %*s %d", &N_puntos_definidos);
        N_puntos_definidos--;
        fscanf(setup,"%*s %*s");
        nodos_definidos.resize(N_puntos_definidos+1);
        for (int i{ 0 }; i<N_puntos_definidos; i++)
        {
            fscanf(setup," %d", &aux);
            nodos_definidos[i] = aux;
        }
        nodos_definidos[N_puntos_definidos] = N_nodos-1;
        P.resize(N_puntos_definidos+1);
        fscanf(setup,"%*s %*s");
        for (int i{ 0 }; i<N_puntos_definidos+1; i++)
        {
            fscanf(setup," %lf", &aux2);
            P[i] = aux2;
        }
        T.resize(N_puntos_definidos+1);
        fscanf(setup,"%*s %*s");
        for (int i{ 0 }; i<N_puntos_definidos+1; i++)
        {
            fscanf(setup," %lf", &aux2);
            T[i] = aux2;
        }
        v.resize(N_puntos_definidos+1);
        fscanf(setup,"%*s %*s");
        for (int i{ 0 }; i<N_puntos_definidos+1; i++)
        {
            fscanf(setup," %lf", &aux2);
            v[i] = aux2;
        }
        fscanf(setup,"%*s");
        fscanf(setup,"%*s %*s %lf", &Amax);
        fscanf(setup,"%*s %*s %lf", &Amin);
        fscanf(setup,"%*s");
        fscanf(setup,"%*s %*s %d", &aux);

        fclose(setup);

        rho.resize(N_puntos_definidos+1);
        for (int i{0}; i < N_puntos_definidos; i++)
        {
            rho[i] = P[i]/(R*T[i]);
        }
        L = (N_nodos-1)*delta_x;
        method = static_cast<METHOD>(aux);
    }

};

extern Settings s;

struct Nodo // {frontera, entrada, x, A, U1, U2, U3, rho, v, T, P, a, H, E, M}
{
    bool frontera{};    // ¿El nodo es frontera?
    bool entrada{};     // Si es frontera, ¿es entrada?
    double x{};         // Posición
    double A{};         // Área

    // VARIABLES CONSERVADAS                
    double U1{};        // rho*A
    double U2{};        // rho*A*v
    double U3{};        // E*A
    // VARIABLES PRIMITIVAS
    double rho{};       // Densidad
    double v{};         // Velocidad
    double T{};         // Temperatura
    double P{};         // Presión
    double a{};         // Velocidad del sonido
    double H{};         // Entalpía total
    double E{};         // Energía total
    double M{};         // Mach
    
};

struct Gradiente
{
    double p_U1{};
    double p_U2{};
    double p_U3{};
};

inline void promediosRoe(const Nodo& n_L, const Nodo& n_R, double& v_prom, double& H_prom, double& a_prom);

inline void flujo(const Nodo& n, std::array<double, 3>& F);
inline void flujoHLL(const Nodo& n_L, const Nodo& n_R, const double& v_prom, const double& a_prom, std::array<double, 3>& F);
inline void flujoHLLS(const Nodo& n_L, const Nodo& n_R, const double& v_prom, const double& H_prom, const double& a_prom, std::array<double, 3>& F, const bool& minus);
inline void flujoHLLC(const Nodo& n_L, const Nodo& n_R,
            const double& v_prom, const double& a_prom, std::array<double, 3>& F);
inline void flujoHLLCS(const Nodo& n_L, const Nodo& n_R,
            const double& v_prom, const double& H_prom, const double& a_prom, std::array<double, 3>& F, const bool& minus);

inline void SolverMacCormack(const Nodo& n_L, const Nodo& n_R,
                    Gradiente& dU, const bool& pred);

inline void SolverRoeRaw(const Nodo& n_L, const Nodo& n_R, //i, i+1 Sin la corrección de la entropía.
					Gradiente& dU_L, Gradiente& dU_R);

inline void SolverRoe(const Nodo& n_L, const Nodo& n_R,
                    Gradiente& dU_L, Gradiente& dU_R);

inline void SolverHLL(const Nodo& n_L, const Nodo& n_C, const Nodo& n_R,
                    Gradiente& dU);

inline void SolverHLLS(const Nodo& n_L, const Nodo& n_C, const Nodo& n_R, //i, i+1
                    Gradiente& dU);

inline void SolverHLLC(const Nodo& n_L, const Nodo& n_C, const Nodo& n_R,
                    Gradiente& dU);

inline void SolverHLLCS(const Nodo& n_L, const Nodo& n_C, const Nodo& n_R,
                    Gradiente& dU);

class Timer
{
private:
	// Type aliases to make accessing nested type easier
	using clock_t = std::chrono::high_resolution_clock;
	using second_t = std::chrono::duration<double, std::ratio<1> >;
	
	std::chrono::time_point<clock_t> m_beg;
 
public:
	Timer() : m_beg(clock_t::now())
	{
	}
	
	void reset()
	{
		m_beg = clock_t::now();
	}
	
	double elapsed() const
	{
		return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
	}
};

class Sonda
{
public:
    friend class Malla;
private:
    int llamadas;
	std::vector<Nodo> punto_sondeo{};
	std::vector<double> tiempos{}; 

public:
    Sonda():llamadas{ 0 } {}

    void escribirSondeo(std::ofstream& f);

};

class Frontera
{
private:
    const Settings s{};
    const bool entrada{};   // Indica si es frontera de entrada o de salida
    const bool libre{};     // Indica si la salida o entrada son libres o se fijan variables
    
    Nodo n_ghost{};

public:

    friend class Malla;

    Frontera() = default;
    Frontera(const bool& def_entrada, const Settings& set, bool def_libre):
    s{ set }, entrada{ def_entrada }, libre{def_libre}
    {
        Frontera();
        if (entrada == true)
        {
            n_ghost.rho = s.rho[0]; 
            n_ghost.v = s.v[0];
            n_ghost.T = s.T[0]; 
            n_ghost.P = s.P[0];
            n_ghost.x = 0;
        }
        else
        {
            n_ghost.rho = s.rho[s.N_puntos_definidos-1]; 
            n_ghost.v = s.v[s.N_puntos_definidos-1];
            n_ghost.T = s.T[s.N_puntos_definidos-1]; 
            n_ghost.P = s.P[s.N_puntos_definidos-1];
            n_ghost.x = s.delta_x*s.N_nodos;
        }
        
        n_ghost.A = areaTubo(n_ghost.x);
        // Primer paréntesis es Cv, además se E_k = 0, porque en reservorio v=0 ESTO ES UN POCO DELICADO
        n_ghost.H = (R/(Gamma-1))*n_ghost.T + n_ghost.P/n_ghost.rho; //+ 0.5*n_ghost.v*n_ghost.v; 
        n_ghost.E = n_ghost.P/(Gamma-1); //+ 0.5*n_ghost.rho*n_ghost.v*n_ghost.v;
        n_ghost.a = sqrt(Gamma*n_ghost.P/n_ghost.rho);
        n_ghost.M = n_ghost.v/n_ghost.a;
                    
        n_ghost.U1 = n_ghost.rho*n_ghost.A;
        n_ghost.U2 = n_ghost.rho*n_ghost.A*n_ghost.v;
        n_ghost.U3 = n_ghost.E*n_ghost.A;
    }

    // Interpolación lineal
    inline void actualizaFrontera(Nodo& n0, const Nodo& n1, const Nodo& n2, const double& del_t_x, const METHOD method = METHOD::ROE);
    // Interpolación cuadrática
    //inline void actualizaFrontera(Nodo& n0, const Nodo& n1, const Nodo& n2, const Nodo& n3, const METHOD method = METHOD::ROE);
};

class Malla
{
private:

    const Settings s;
    std::vector<Nodo> n;            // Nodos de la malla (centros de las celdas)
    const double delta_x{};         // Paso de malla (equiespaciado)
    double delta_t{};               // Paso de tiempo
    double tiempo{};                // Tiempo actual simulado de la malla

    inline double buscaDelta_t()    // Busca el valor máximo de dt en base a condición CFL
    {
        static double min{};
        static double aux{};
		// En teoría aquí es del máximo autovalor, que he supuesto que siempre es v+a, pero no tiene porque
		min = fabs( 1./(n[0].a + n[0].v) );
        
		for (const Nodo& x : n)
        {
            aux = fabs( 1./(x.a + x.v) );
            if (aux<min)
                min = aux;
        };
        
        return fabs( s.CFL * delta_x * min );
    }

public:

    friend class Frontera;
    friend class Sonda;

    Malla(const Settings& set) : // La malla se construye vacía, pero nodos ya tienen sus posiciones.
    s{ set }, delta_x{ s.delta_x }, tiempo{ 0 }
    {   
        n.resize(s.N_nodos);
        for (int i{ 0 }; i<s.N_nodos; i++)
        {
            n[i].x = delta_x * i;
            n[i].frontera = false;
            n[i].entrada = false;
        }
        n[0].frontera = n[s.N_nodos-1].frontera = true;
        n[0].entrada = true;
    }

	void inicializar();
    // Si es inline parece que no puede usarse de fuera del cpp methods.cpp
    // NOTA: Parece que los algoritmos son más estables si no se desarrolla el álgebra y se dejan los términos 
    // separados por paréntesis (flujo) + (término fuente)

	void calcularNodos(Frontera& in_flow, Frontera& out_flow, const METHOD method = METHOD::ROE);   // FUNCIÓN PRINCIPAL

	void escribirDatos(std::ofstream& f, const int precision);  // Input: fichero y el número de decimales con el que se mostrarán los números

	void sondear(Sonda& sonda, const int numero_nodo);          // Registra datos y los guarda en la sonda

	inline double tiempoTrans() { return tiempo; }
};

#endif