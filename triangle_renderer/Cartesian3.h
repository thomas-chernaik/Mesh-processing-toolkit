///////////////////////////////////////////////////
//
//	Hamish Carr
//	January, 2018
//
//	------------------------
//	Cartesian3.h
//	------------------------
//	
//	A minimal class for a point in Cartesian space
//	
///////////////////////////////////////////////////

#ifndef CARTESIAN3_H
#define CARTESIAN3_H

#include <iostream>

// the class - we will rely on POD for sending to GPU
class Cartesian3
	{ // Cartesian3
	public:
	// the coordinates
	float x, y, z;

	// constructors
	Cartesian3();
	Cartesian3(float X, float Y, float Z);
	
	// equality operator
	bool operator ==(const Cartesian3 &other);

	// addition operator
	Cartesian3 operator +(const Cartesian3 &other);

	// subtraction operator
	Cartesian3 operator -(const Cartesian3 &other);
	
	// multiplication operator
	Cartesian3 operator *(float factor);

	// division operator
	Cartesian3 operator /(float factor);

	// crossproduct routine
	Cartesian3 cross(const Cartesian3 &other);
	
	// routine to find the length
	float length();
	
	// normalisation routine
	Cartesian3 normalise();

	}; // Cartesian3

// stream output
std::ostream & operator << (std::ostream &outStream, Cartesian3 value);
		
#endif