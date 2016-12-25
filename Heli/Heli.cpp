#include <cmath>
#include "Heli.h"

#define MAX(x,y) (x < y ? y : x)

const double PI = 4 * atan(1.0);

Heli::Heli(double X0[], int Xn, double U0[], int Un, double t0){

	X = new double[Xn];
	this->Xn = Xn;
	U = new double[Un];
	this->Un = Un;
	t = t0;

	for (int i = 0; i < Xn; i++) {
		X[i] = X0[i];
	}
	for (int i = 0; i < Un; i++) {
		U[i] = U0[i];
	}
}

Heli::~Heli(){
	if (X) delete X;
	if (U) delete U;
}
// dt = step time; Ttr = tail rotor thrust; Tmr = Main rotor thrust; *q = output q[0]:x_fly and q[1]:
void Heli::sim_step(double dt, double Ttr, double Tmr, double d_lat, double d_lon, double *q_array, bool collided) {
	double u, v, w, p, q, r, phi, theta, psi; // state variables

	t += dt;
	// unpacking states
	u = X[0];
	v = X[1];
	w = X[2];
	p = X[3];
	q = X[4];
	r = X[5];
	phi = X[6];
	theta = X[7];
	psi = X[8];


	// define paramaters 
	// many unused due to physics simplification
	double g = 9.81; // gravity
	double m = 8.2; // mass
	double rho_a = 1.225; // air density
	double u_w = 0.0; // x wind velocity
	double v_w = 0.0; // y wind velocity
	double w_w = 0.0; // z wind velocity
	double S_vf = 0.012; // Vertical fin area
	double S_ht = 0.01; // Horizontal stabilizer area
	double htr = 0.08; // T.R. z-distance from C.G.
	double k_mr_x = -0.3;
	double l_tr = 0.91; // T.R. x-distance from C.G.
	double l_ht = 0.71;
	double S_fus_x = 0.1;
	double S_fus_y = 0.22;
	double S_fus_z = 0.15;
	double u_a = u - u_w;
	double v_a = v - v_w;
	double w_a = w - u_w;
	double V_inf = sqrt(u_a*u_a + v_a*v_a); //
	double Ixx = 0.18;
	double Iyy = 0.34;
	double Izz = 0.28;
	double Xmr = -Tmr*k_mr_x*sin(theta); // Main rotor X component of force
	double Ymr = Tmr*k_mr_x*sin(phi)*cos(theta); // Main rotor Y component of force
	double Zmr = -Tmr*cos(phi)*cos(theta); // Main rotor Z component of force
	double Ytr = 0.0; //
	double Ltr = Ytr*htr;
	double Ntr = -Ttr*l_tr; //
	double Qe = r*0.5; // Engine moment
	double Xfus = u < 0.0 ? 0.5*rho_a*S_fus_x*pow(5.0*u_a, 2) : -0.5*rho_a*S_fus_x*pow(5.0*u_a, 2); // limit u
	double Yfus = v < 0.0 ? 0.5*rho_a*S_fus_y*pow(5.0*v_a, 2) : -0.5*rho_a*S_fus_y*pow(5.0*v_a, 2); // limit v
	double Nfus = r < 0.0 ? 0.5*rho_a*S_vf*pow(50.0*r,2) : -0.5*rho_a*S_vf*pow(50.0*r, 2); // limit r
	double Lfus = p < 0.0 ? 0.5*rho_a*S_ht*pow(10.0*p, 2) : -0.5*rho_a*S_ht*pow(10.0*p, 2); // Limit p
	double Mfus = q  < 0.0 ? 0.5*rho_a*S_ht*pow(5.0*q, 2) : -0.5*rho_a*S_ht*pow(5.0*q, 2); // Limit q
	double Zfus = w < 0.0 ? 0.5*rho_a*S_fus_z*pow(10.0*w_a, 2) : -0.5*rho_a*S_fus_z*pow(10.0*w_a, 2); // Limit w
	double Lst = -phi*10.0; // Pseudo roll control (prevents excessive roll angles)
	double Mst = -theta*10.0; // Pseudo pitch control (prevents excessive pitch angles)
	double Yvf = 0.0; // Vertical fin Y force component 
	double Lvf = 0.0; // Vertical fin moment about X component 
	double Nvf = 0.0; // Vertical fin moment about Z component 
	double Mht = 0.0;
	double Zht = 0.0;



	// pack inputs
	U[0] = Tmr;
	U[1] = Ttr;

	// derivatives 
	double u_dot, v_dot, w_dot, p_dot, q_dot, r_dot, phi_dot, theta_dot, psi_dot;

	// calculate derivatives 
	if (collided) {
		u_dot = 0.0;
		v_dot = 0.0;
		w_dot = 0.1*g*cos(phi)*cos(theta) + (Zmr + Zfus + Zht) / m < 0.0 ? 0.1*g*cos(phi)*cos(theta) + (Zmr + Zfus + Zht) / m : 0.0;
		p_dot = 0.0;
		q_dot = 0.0;
		r_dot = 0.0;
		phi_dot = 0.0;
		theta_dot = 0.0;
		psi_dot = 0.0;

		w = w + w_dot*dt < 0.0 ? w + w_dot*dt : 0.0;
		u = 0.0;
		v = 0.0;
		psi = psi;
		theta = theta;
		phi = phi;

		X[0] = u;
		X[1] = v;
		X[2] = w;
		X[3] = p;
		X[4] = q;
		X[5] = r;
		X[6] = phi;
		X[7] = theta;
		X[8] = psi;
		
	}
	else {
		// Physics simplified due to conflicting signs causing variables to go to infinity nearly instantly 
		u_dot = 0.1*g*sin(theta) + (Xmr + Xfus) / m;  //v*r - w*q - g*sin(theta) + (Xmr + Xfus) / m; 
		v_dot = -0.06*g*sin(phi)*cos(theta) + (Ymr + Yfus + Ytr + Yvf) / m;  //w*p - u*r + g*sin(phi)*cos(theta) + (Ymr + Yfus + Ytr + Yvf) / m; 
		w_dot = 0.1*g*cos(phi)*cos(theta) + (Zmr + Zfus + Zht) / m; //u*p - v*p + g*cos(phi)*cos(theta) + (Zmr + Zfus + Zht) / m; 
		p_dot = (d_lat + Lfus) / Ixx;//q*r*(Iyy-Izz)/Ixx + (Lmr + Lvf + Ltr)/Ixx; // roll accel
		q_dot = (d_lon + Mfus) / Iyy;//p*r*(Izz-Ixx)/Iyy + (Mmr + Mht) / Iyy;  // pitch accel
		r_dot = (-Qe + Nvf + Ntr + Nfus) / Izz; //p*q*(Ixx-Iyy)/Izz + (-Qe + Nvf + Ntr + Nfus) / Izz; // yaw accel
		phi_dot = p + Lst; // roll vel
		theta_dot = q + Mst; // pitch vel
		psi_dot = r; // yaw vel

		X[0] += u_dot*dt;
		X[1] += v_dot*dt;
		X[2] += w_dot*dt;
		X[3] += p_dot*dt;
		X[4] += q_dot*dt;
		X[5] += r_dot*dt;
		X[6] = abs(X[6] + phi_dot*dt) <= 0.005 ? 0.0 : X[6] + phi_dot*dt; // To help it go back to natural state faster
		X[7] = abs(X[7] + theta_dot*dt) <= 0.001 && abs(X[7] + theta_dot*dt) >= -0.001 ? 0.000 : X[7] + theta_dot*dt; // To help it go back to natural state faster
		X[8] += psi_dot*dt;
	}

	q_array[0] = X[0]; //u
	q_array[1] = X[1]; //v
	q_array[2] = X[2]; //w
	q_array[3] = X[6]; //phi
	q_array[4] = X[7]; //theta
	q_array[5] = X[8]; //psi

}

