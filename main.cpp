

#include <vector>
#include <algorithm>
#include <cmath> // M_PI
#include <complex>
#include <iostream>
#include <array>
#include <cstdio>
#include <chrono>
#include <functional>
#include <iomanip>
#include <sstream>

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

void barrier_schottky_nordheim_roots(double f, double phi, double mu, double& x1_out, double& x2_out) {
	double a = -f; // q_e*V/nm = eV/nm
	double b = mu + phi;  // mu + phi; eV
	double c = -constants::q_e*constants::q_e/(16*constants::pi*constants::eps0); // eV*nm
	x1_out = (-b+std::sqrt(b*b-4*a*c))/(2*a);
	x2_out = (-b-std::sqrt(b*b-4*a*c))/(2*a);
}

/** x - nm; f - V/nm; min_pot - eV */
double barrier_schottky_nordheim(double x, double f, double phi = 4.5, double mu = 7.0, double min_pot = 0.0) {

	double a = -f; // q_e*V/nm = eV/nm
	double b = mu + phi;  // mu + phi; eV
	double c = -constants::q_e*constants::q_e/(16*constants::pi*constants::eps0); // eV*nm

	double x1, x2;
	barrier_schottky_nordheim_roots(f, phi, mu, x1, x2);

	if (x < 0.0) return 0.0;
	double pot = b + a*(x+x1) + c/(x+x1);
	if (pot < min_pot) return min_pot;
	return pot;
}

std::complex<double> calculate_k(double energy, double position, std::function<double(double)> barrier) {

	double potential = barrier(position);
	double kin_energy = energy-potential;
	if (kin_energy == 0) kin_energy = 1e-14;
	return std::sqrt(std::complex<double>(kin_energy/constants::const_c, 0));
}

double calculate_transmission(double energy, double xmin, double xmax,
							  int num_regions, std::function<double(double)> barrier) {

	double dx = (xmax-xmin)/num_regions;

	// start the transfer matrix as identity
	// Based on the refs: transfer_matrix = (T1 T3)
	//										(T2 T4)
	std::array<std::complex<double>, 4> transfer_matrix = {1., 0., 0., 1.};

	// potential and k in first region
	std::complex<double> k1 = calculate_k(energy, xmin + 0.5*dx, barrier);

	std::complex<double> k_first = k1;
	std::complex<double> k2;

	// starting from the 2nd region, loop over all transitions
	for (int i = 1; i < num_regions; i++) {
		// potential and k in the next region of the transition
		k2 = calculate_k(energy, xmin + (i+0.5)*dx, barrier);

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

	//double reflection = (std::abs(transfer_matrix[2])*std::abs(transfer_matrix[2])) /
	//					(std::abs(transfer_matrix[0])*std::abs(transfer_matrix[0]));
	/*
	printf("T:   %10.3e\n", transmission);
	printf("R:   %10.3e\n", reflection);
	printf("T+R: %10.3e\n", transmission + reflection);
	*/
	return transmission;
}

template <typename T>
std::string to_string_prec(const T a_value, const int n = 1) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(n) << a_value;
    return out.str();
}

int main() {

	double field = 4.0;
	double phi = 4.5;
	double mu = 60.0;

	double emin = 50.0;
	double emax = 70.0;

	std::string params = "_f" + to_string_prec(field) + "_p" + to_string_prec(phi)
						 + "_m" + to_string_prec(mu) + "_em" + to_string_prec(emax);
	std::string tunnel_file = "./data/tunnel_data"+params+".txt";
	std::string potential_file = "./data/potential"+params+".txt";

	// Roots of the potential:
	double x1, x2;
	barrier_schottky_nordheim_roots(field, phi, mu, x1, x2);

	// The coordinate system is chosen such that the first root is x=0
	// and thus, the second root is at x2-x1
	// corresponding integration limits
	double xmin = - 1.0;
	double xmax = x2 - x1 + 1.0;

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	FILE * outf = fopen(tunnel_file.c_str(), "w");
	for (double energy = emin; energy <= emax; energy += 0.01) {
		auto barrier = [&field, &phi, &mu](double x) {return barrier_schottky_nordheim(x, field, phi, mu);};
		double transm = calculate_transmission(energy, xmin, xmax, 40000, barrier);

		printf("%.5f %.10e\n", energy, transm);
		fprintf(outf, "%.5f %.10e\n", energy, transm);
	}
	fclose(outf);
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()/1e6 << std::endl;

	FILE * file = fopen(potential_file.c_str(), "w");
	for (double x = xmin; x < xmax; x+=0.01) {
		fprintf(file, "%.2f %.5f\n", x, barrier_schottky_nordheim(x, field, phi, mu));
	}
	fclose(file);

	return 0;
}
