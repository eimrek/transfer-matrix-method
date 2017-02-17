

#include <vector>
#include <algorithm>
#include <cmath> // M_PI
#include <complex>
#include <iostream>
#include <array>
#include <cstdio>
#include <chrono>

namespace constants {
	const double hbar = 6.582119514e-16; // eV*s
	const double mass_e = 5.6856300620910916e-30; // eV*s^2/nm^2
	const double q_e = 1.6021766208e-19; // C
	const double eps0 = 1.4185972717563562e-39; // C^2/(eV*nm)

	const double pi = 3.14159265358979323846;

	const double const_c = hbar*hbar/(2*mass_e); // eV*nm^2
}

double barrier_rectangle(double x) {
	if (x > 1.0 && x < 2.0) return 1.0;
	if (x >= 2.0) return 0.0;
	return 0.0;
}


double min_pot = -4.0;
/** x - nm; f - V/nm */
double barrier_schottky_nordheim(double x, double f) {
	double b = 11.5;  // mu + phi; eV
	double a = -f; // q_e*V/nm = eV/nm
	double c = -constants::q_e*constants::q_e/(16*constants::pi*constants::eps0); // eV*nm

	double x_start = 1.0;
	double x1 = x_start - (-b+std::sqrt(b*b-4*a*c))/(2*a);
	double x2 = x1 + (-b-std::sqrt(b*b-4*a*c))/(2*a);

	//std::cout << x1 << " " << x2 << std::endl;

	if (x < x_start) return 0.0;
	double pot = b + a*(x-x1) + c/(x-x1);
	if (pot < min_pot) return min_pot;
	return pot;

}

std::complex<double> calculate_k(double energy, double position, double (*barrier)(double x)) {

	double potential = (*barrier)(position);
	double kin_energy = energy-potential;
	if (kin_energy == 0) kin_energy = 1e-14;
	return std::sqrt(std::complex<double>(kin_energy/constants::const_c, 0));
}

double calculate_transmission(double energy, double xmin, double xmax, int num_regions, double (*barrier)(double x)) {

	double dx = (xmax-xmin)/num_regions;

	// start the transfer matrix as identity
	// Based on the refs: transfer_matrix = (T1 T3)
	//										(T2 T4)
	std::array<std::complex<double>, 4> transfer_matrix = {1., 0., 0., 1.};

	// potential and k in first region
	std::complex<double> k1 = calculate_k(energy, 0.5*dx, barrier);

	std::complex<double> k_first = k1;
	std::complex<double> k2;

	// starting from the 2nd region, loop over all transitions
	for (int i = 1; i < num_regions; i++) {
		// potential and k in the next region of the transition
		k2 = calculate_k(energy, (i+0.5)*dx, barrier);

		// discontinuity matrix elements
		std::complex<double> dm1 = 0.5*(1. + k2/k1);
		std::complex<double> dm2 = 0.5*(1. - k2/k1);

		//std::cout << dm1 << " " << dm2 << std::endl;

		// multiply the accumulated transfer matrix with the current discontinuity matrix
		std::complex<double> tm0 = transfer_matrix[0];
		std::complex<double> tm1 = transfer_matrix[1];
		transfer_matrix[0] = tm0*dm1 + transfer_matrix[2]*dm2;
		transfer_matrix[1] = tm1*dm1 + transfer_matrix[3]*dm2;
		transfer_matrix[2] = tm0*dm2 + transfer_matrix[2]*dm1;
		transfer_matrix[3] = tm1*dm2 + transfer_matrix[3]*dm1;

		// propagation matrix elements
		std::complex<double> pm1 = std::exp(-std::complex<double>(0., 1.)*k2*dx);
		std::complex<double> pm2 = std::exp( std::complex<double>(0., 1.)*k2*dx);

		//std::cout << pm1 << " " << pm2 << std::endl;

		// multiply the accumulated transfer matrix with the current propagation matrix
		transfer_matrix[0] = transfer_matrix[0]*pm1;
		transfer_matrix[1] = transfer_matrix[1]*pm1;
		transfer_matrix[2] = transfer_matrix[2]*pm2;
		transfer_matrix[3] = transfer_matrix[3]*pm2;

		// move onto the next transition
		k1 = k2;
	}

	/*
	std::cout << std::endl;
	for (int j = 0; j < 4; j++) {
		std::cout << transfer_matrix[j] << std::endl;
	}
	std::cout << std::endl;
	*/

	double transmission = k2.real()/k_first.real() *
						  1./(std::abs(transfer_matrix[0])*std::abs(transfer_matrix[0]));
	double reflection = (std::abs(transfer_matrix[2])*std::abs(transfer_matrix[2])) /
						(std::abs(transfer_matrix[0])*std::abs(transfer_matrix[0]));
	/*
	printf("T:   %10.3e\n", transmission);
	printf("R:   %10.3e\n", reflection);
	printf("T+R: %10.3e\n", transmission + reflection);
	*/
	return transmission;
}

int main() {
	//double transmission =  calculate_transmission(0.3);
	//printf("T:   %10.3e\n", transmission);

	/*
	auto barrier = [](double x) {return barrier_schottky_nordheim(x, 1.0, -5.0);};

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 5; i++) {
		double xmax = 14.0+i*2.0;
		printf("%.1f %.10e\n", xmax, calculate_transmission(10.0, 0.0, xmax, 2000000, barrier));
	}
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	*/
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	FILE * outf = fopen("./data/outf.txt", "w");
	for (int i = 0; i < 1500; i++) {
		double energy = 0.01*i;
		auto barrier = [](double x) {return barrier_schottky_nordheim(x, 1.0);};
		printf("%.5f %.10e\n", energy, calculate_transmission(energy, 0.0, 20.0, 20000, barrier));
		fprintf(outf, "%.5f %.10e\n", energy, calculate_transmission(energy, 0.0, 20.0, 20000, barrier));
	}
	fclose(outf);
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()/1e6 << std::endl;

	FILE * file = fopen("./data/test.txt", "w");
	for (int i = 0; i < 400; i++) {
		double x = i*0.1;
		fprintf(file, "%.2f %.5f\n", x, barrier_schottky_nordheim(x, 1.0));
	}
	fclose(file);

	return 0;
}
