
#pragma once

#ifndef Heli_H
#define Heli_H

class Heli {
	double *X; //
	double *U; //
	int Xn, Un; // Size of arrays X and U
	double t; // sim time

public:

	Heli( double X0[], int Xn, double U[], int Un, double t0);
	~Heli();
	void sim_step(double dt, double Ttr, double Tmr, double d_lat, double d_lon, double *q, bool collided);

};
#endif
//Heli::Heli(char *mesh_filename, double X0[], int Xn, double U0[], int Un, double t0);

//Heli::~Heli();
// dt = step time; Ttr = tail rotor thrust; Tmr = Main rotor thrust; *q = output q[0]:x_fly and q[1]:
//void Heli::sim_step(double dt, double Ttr, double Tmr, double *q);
