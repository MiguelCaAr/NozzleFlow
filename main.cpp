#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "constants.h"
#include "methods.h"

Settings s{};
bool PAUSA = false;

double areaTubo(const double x)
{
	return ( s.Amin + 4.*(s.Amax-s.Amin)/(s.L*s.L) * (x*x - x*s.L + 0.25*s.L*s.L) );
}

int main()
{
	Timer t;
	Malla malla{ s };
	Frontera in_flow{ true, s, s.entrada_libre }, out_flow{ false, s, s.salida_libre };
	Sonda sonda{};
	std::ofstream salida; // Fichero de salida
	int pasos_sim{ 0 };


	printf("\n>> Simulation setup loaded");
	printf("\n   Simulation time: %.5lf s", s.simTime);
	printf("\n   CFL: %.2lf", s.CFL);
	printf("\n   Number of cells: %d", s.N_nodos);
	printf("\n   dif_x: %.5lf m", s.delta_x);
	printf("\n   Pipe length: %.5lf m", s.L);
	switch (s.method)
	{
		case METHOD::MACCORMACK:
			std::cout << "\n   Using MacCormack's solver";
			break;
		case METHOD::ROE:
			std::cout << "\n   Using Roe's solver";
			break;
		case METHOD::HLL:
			std::cout << "\n   Using HLL solver";
			break;
		case METHOD::HLLC:
			std::cout << "\n   Using HLLC solver";
			break;
		default:
			std::cout << "\n   No valid method selected, aborting";
			exit(1);
			break;
	}


	malla.inicializar();

	salida.open("salida.dat");

    malla.escribirDatos(salida, 8);
	if (s.is_sonda == true) malla.sondear(sonda, s.N_nodo_sonda);

	std::cout << "\n>> Starting simulation";

	while ( malla.tiempoTrans() <= s.simTime ) // NÚCLEO DEL PROGRAMA
	{
		malla.calcularNodos(in_flow, out_flow, s.method);

		pasos_sim++;
	
		if ( ( pasos_sim % s.intervalo_extraer ) == 0)
        {
            malla.escribirDatos(salida, 8);
            if (s.is_sonda == true) malla.sondear(sonda, s.N_nodo_sonda); 
        };
	}

	malla.escribirDatos(salida, 8);
	if (s.is_sonda == true) malla.sondear(sonda, s.N_nodo_sonda); 

	salida.close();
	if (s.is_sonda == true) 
	{
		salida.open("sonda.dat");
		sonda.escribirSondeo(salida);
		salida.close();
	}

	std::cout << "\n\n>> Simulation completed!";
	std::cout << "\n   Time steps simulated: " << pasos_sim;
	std::cout << "\n   Time elapsed: " << t.elapsed() << " seconds";
	std::cout << "\n\n";

	return 0;
}