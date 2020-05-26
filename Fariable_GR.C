
//This is a C code for calculating the motion of the sun earth mercury system as two seperate 2-body problems.
//Written by Daniel Foulds-Holt 
//df00177@surrey.ac.uk

//== INCLUDES ==
#include <stdio.h>
#include <math.h>
#include <string.h>

//== GLOBALS ==
double G = 1.0; //Gravitational constant
double c; // Set later in code
double pi = 3.14159265359;

int TimeLength = 1000;

//== STRUCTURES ==

// Data structure encompassing all information about a planet (, sun, or COM) 
struct Status{

	//Empty constructor
	Status() {
		x    = 0.0;
		v_x  = 0.0;
		y    = 0.0;
		v_y  = 0.0;
		z    = 0.0;
		v_z  = 0.0;
		Mass = 0.0;
		TimeElapsed = 0.0;
		No_Orbits = 0;
		e = 0;
		a = 0;
	}

	//Full constructor
	Status(double xx, double yy, double zz, double vx, double vy, double vz, double M, double T) {
		x    = xx;
		v_x  = vx;
		y    = yy;
		v_y  = vy;
		z    = zz;
		v_z  = vz;
		Mass = M;
		TimeElapsed = T;
		No_Orbits = 0;
		e = 0;
		a = 0;
	}

	//Prints Status to Screen
	void Print() {
		printf("\nStatus:\n");
		printf("\tx  = %.4E\n", x);
		printf("\tvx = %.4E\n",v_x);
		printf("\ty  = %.4E\n", y);
		printf("\tvy = %.4E\n",v_y);
		printf("\tz  = %.4E\n", z);
		printf("\tvz = %.4E\n",v_z);
		printf("\tMass = %.4E\n", Mass);
		printf("\tTimeElapsed = %.4E\n", TimeElapsed);
		printf("\n");
	}


	//Variables ====================================

	double x;
	double v_x;
	double y;
	double v_y;
	double z;
	double v_z;
	double Mass;
	double TimeElapsed;
	int No_Orbits;
	double e;
	double a;
};

// Data structure containing the initial conditions of a planet - used for passing back and forth
struct InitialConditions{

	//Constructor
	InitialConditions(double M, double A, double P, double E, double V_a, double V_p) {
		Mass = M;
		Ap   = A;
		Pe   = P;
		Ecc  = E;
		Va   = V_a;
		Vp   = V_p;
	}

	//Variables ===========

	double Mass;
	double Ap;
	double Pe;
	double Ecc;
	double Va;
	double Vp;
};

// Data structure for holding vectors 
struct Vect{

	//Empty constructor
	Vect() {
		x    = 0.0;
		y    = 0.0;
		z    = 0.0;
	}

	//Constructor
	Vect(double X_in, double Y_in, double Z_in){
		x = X_in;
		y = Y_in;
		z = Z_in;
	}

	void Print(){
		printf("[%f, %f, %f]\n",x,y,z);
	}

	double Length(){
		double LSqr = pow(x,2.0) + pow(y,2.0) + pow(z,2.0);
		double L = pow(LSqr, 0.5);
		return L;
	} 

	//Variables ===========

	double x; 
	double y;
	double z;
};

// This is the semimajor axis returning thing
struct SemiEcc{

	SemiEcc(double Ain, double Ein){
		A = Ain;
		E = Ein;
	}

	double A;
	double E;
};

// Parameters to be returned after calculating the orbit
struct OrbitalParameters{

	OrbitalParameters(Vect Ecc1, Vect Ecc2, SemiEcc SemEcc1, SemiEcc SemEcc2){
		EL[0] = Ecc1.Length();
		EL[1] = Ecc2.Length();
		EPhi[0] = atan(Ecc1.y / Ecc1.x);
		EPhi[1] = atan(Ecc2.y / Ecc2.x);
		DeltaPhi = EPhi[1] - EPhi[0];
		E[0] = SemEcc1.E;
		E[1] = SemEcc2.E;
		A[0] = SemEcc1.A;
		A[1] = SemEcc2.A;
	}

	double E[2];
	double A[2];
	double EL[2];
	double EPhi[2];
	double DeltaPhi;

};

//== VECTOR FUNCTIONS ==

//Allows for the cross product to be calculated between two vectors
Vect CrossProduct(Vect A, Vect B){
	Vect C = Vect();
	C.x = A.y * B.z - A.z * B.y;
	C.y = A.x * B.z - A.z * B.x;
	C.z = A.x * B.y - A.y * B.x;
	return C;
}

//Allows for the dot product to be calculated between two vectors
double DotProduct(Vect A, Vect B){
	double C = A.x * B.x + A.y * B.y + A.z * B.z;
	return C;
}

//Allows for two vectors to be added 
Vect VectAdd(Vect A, Vect B, double factor){
	Vect C = Vect();
	C.x = A.x  + (factor * B.x);
	C.y = A.y  + (factor * B.y);
	C.z = A.z  + (factor * B.z);
	return C;
}

//Allows for a vector to be scaler
Vect VectScale(Vect A, double factor){
	Vect C = Vect();
	C.x = (factor * A.x);
	C.y = (factor * A.y);
	C.z = (factor * A.z);
	return C;
}


//== FUNCTIONS ==

// Outputs the initial conditions to a file
//int OutputIC(InitialConditions IC, char* Filename){ \\ <--- For the old initial conditions file
int OutputIC(Status Planet1, Status Planet2, char* Filename){

	FILE *fp;
	fp = fopen(Filename, "w+");

	fprintf(fp, "%f, ", Planet1.x);
	fprintf(fp, "%f, ", Planet1.y);
	fprintf(fp, "%f, ", Planet1.z);
	fprintf(fp, "%f, ", Planet1.v_x);
	fprintf(fp, "%f, ", Planet1.v_y);
	fprintf(fp, "%f " , Planet1.v_z);
	fprintf(fp, "\n");

	fprintf(fp, "%f, ", Planet2.x);
	fprintf(fp, "%f, ", Planet2.y);
	fprintf(fp, "%f, ", Planet2.z);
	fprintf(fp, "%f, ", Planet2.v_x);
	fprintf(fp, "%f, ", Planet2.v_y);
	fprintf(fp, "%f " , Planet2.v_z);
	fprintf(fp, "\n");

	fclose(fp);
	return 1;
}

// Prints the initial conditions to the screen 
int DisplayIC(InitialConditions Earth, InitialConditions Merc) {

	printf("\n == Initial Conditions ==\n\n");

	printf("Earth:\n");
	printf("\tMass           = %.4E Solar Massess\n", Earth.Mass);
	printf("\tApoapsis       = %.4E AU\n", Earth.Ap);
	printf("\tPeriapsis      = %.4E AU\n", Earth.Pe); 
	printf("\tEccentricity   = %.4E\n", Earth.Ecc);
	printf("\n");
	printf("\tApoapsis Vel.  = %.4E Units\n", Earth.Va);
	printf("\tPeriapsis Vel. = %.4E Units\n", Earth.Vp);
	printf("\n");

	printf("Mercury:\n");
	printf("\tMass           = %.4E Solar Massess\n", Merc.Mass);
	printf("\tApoapsis       = %.4E AU\n", Merc.Ap);
	printf("\tPeriapsis      = %.4E AU\n", Merc.Pe);
	printf("\tEccentricity   = %.4E\n", Merc.Ecc);
	printf("\n");
	printf("\tApoapsis Vel.  = %.4E Units\n", Merc.Va);
	printf("\tPeriapsis Vel. = %.4E Units\n", Merc.Vp);
	printf("\n");

	return 1;
}

// Finds the Centre Of Mass (COM) for a "Sun-Planet" system and outputs the Status 
Status FindCOM(Status Sun, Status Planet) {

	Status Output;
	
	Output.x = (Planet.Mass * Planet.x + Sun.Mass * Sun.x) / (Planet.Mass + Sun.Mass);
	Output.y = (Planet.Mass * Planet.y + Sun.Mass * Sun.y) / (Planet.Mass + Sun.Mass);
	Output.z = (Planet.Mass * Planet.z + Sun.Mass * Sun.z) / (Planet.Mass + Sun.Mass);

	return Output;
}

// Repositions the planets (& Sun) to the COM reference frame
void MoveToCOM(Status *Sun, Status *Planet) { 

	Status COM = FindCOM(*Sun, *Planet);
	//COM.Print();

	double RM  = (Planet->Mass * Sun->Mass) / (Planet->Mass + Sun->Mass);

	// Position =========================
	double Temp;
	//Planet
	Temp = Planet->x - COM.x;
	Planet->x = (RM * Temp) / Planet->Mass;
	Temp = Planet->y - COM.y;
	Planet->y = (RM * Temp) / Planet->Mass;
	Temp = Planet->z - COM.z;
	Planet->z = (RM * Temp) / Planet->Mass;

	//Sun
	Temp = Sun->x - COM.x;
	Sun->x = - (RM * Temp) / Sun->Mass;
	Temp = Sun->y - COM.y;
	Sun->y = - (RM * Temp) / Sun->Mass;
	Temp = Sun->z - COM.z;
	Sun->z = - (RM * Temp) / Sun->Mass;

	// Velocity =========================

	//Planet
	Temp = Planet->v_x;
	Planet->v_x = - (Temp * RM) / Planet->Mass; 
	Temp = Planet->v_y;
	Planet->v_y = - (Temp * RM) / Planet->Mass; 
	Temp = Planet->v_z;
	Planet->v_z = - (Temp * RM) / Planet->Mass; 
	
	//Sun
	Temp = Sun->v_x;
	Sun->v_x = - (Temp * RM) / Sun->Mass; 
	Temp = Sun->v_y;
	Sun->v_y = - (Temp * RM) / Sun->Mass; 
	Temp = Sun->v_z;
	Sun->v_z = - (Temp * RM) / Sun->Mass; 

	// ==================================
}

// Returns the coefficant A
double Return_AGR(Vect R, Vect V, double m1, double m2, int* Order, double c){

	double M = m1 + m2;
	double r_dot =  DotProduct(R, V) / R.Length();
	double nu = m1 * m2 / pow(M,2.0);

	double RL = R.Length();
	double VL = V.Length();

	double a0, a1, a2, a2_1, a2_2, a2_3, a52;

	a0 = 1.0;

	a1 = -1.5 * pow(r_dot,2.0) * nu + (1 + 3 * nu) * pow(VL,2.0) - 2 * (2 + nu) * G * M / RL;

	a2_1 = (15.0 / 8.0) * pow(r_dot, 4.0) * nu * (1 - 3*nu) + 3 * pow(r_dot,2.0) * nu * pow(VL, 2.0) * (2 * nu - 1.5);
	a2_2 = nu * pow(VL,4.0)  * (3 - 4 * nu) + pow((G * M / RL),2.0) * (9 + nu * 87 / 4); 
	a2_3 = (G * M / RL) * (-2.0 * pow(r_dot,2.0) * (1.0 + pow(nu,2.0)) - 25.0 * r_dot * nu - 13.0 * nu * pow(VL,2.0) / 2.0 );
	a2 = a2_1 + a2_2 + a2_3;

	a52 = -8.0 / 5.0 * G * M / RL * nu * r_dot * (17.0 * G * M / (RL * 3.0) - 3 * pow(VL, 2.0));

	double Result = Order[0] * a0 + Order[1] * a1 * pow(c, -2.0) + Order[2] * a2 * pow(c, -4.0) + Order[3] * a52 * pow(c, -5.0);
	return Result;
}

// Returns the coefficant B
double Return_BGR(Vect R, Vect V, double m1, double m2, int* Order, double c){	
	double M = m1 + m2;
	double r_dot =  DotProduct(R, V) / R.Length();
	double nu = m1 * m2 / pow(M,2.0);

	double RL = R.Length();
	double VL = V.Length();

	double b0, b1, b2, b2_1, b2_2, b52;

	b0 = 0.0;
	
	b1 = -2.0 * (2.0 - nu) * r_dot;

	b2_1 = 3.0 * pow(r_dot, 3.0) * nu * (1.5 + nu) - r_dot * nu * pow(RL,2.0) * (15.0 / 2.0 + 2.0 * nu);
	b2_2 = G * M * r_dot / RL * (2.0 + 41.0 * nu / 2.0 + 4.0 * pow(nu,2.0));
	b2 = b2_1 + b2_2;

	b52 = 8.0 * G * M * nu / (5.0 * RL) * (3.0 * G * M / RL + pow(VL,2.0));

	double Result = Order[0] * b0 + Order[1] * b1 * pow(c, -2.0) + Order[2] * b2 * pow(c, -4.0) + Order[3] * b52 * pow(c, -5.0);
	return Result;
}

// Performs an iteration of the Euler method & updates planet positions
void MovePlanets(Status *Sun, Status *Planet, int*GRCorrections, double c, double AP) {

	double NewVX, NewVY, NewVZ;	

	double M = Sun->Mass + Planet->Mass;
	//Vect R = Vect((Sun->x - Planet->x), (Sun->y - Planet->y), (Sun->z - Planet->z));
	Vect R = Vect(Planet->x, Planet->y, Planet->z);
	Vect V = Vect(Planet->v_x, Planet->v_y, Planet->v_z);
	double RL = pow(R.Length(), 3.0);

	Vect N = VectScale(R, 1.0/R.Length());

	double Temp = RL / (G * M);
	double TimePeriod = AP * pow(Temp, 0.5);
	Planet->TimeElapsed += TimePeriod;

	double tempy = Planet->y;

	Planet->a = 1.0 / (2.0 / R.Length() - pow(V.Length(),2.0) / (G * M));
	Vect H = CrossProduct(R,V);
	Planet->e = pow((1.0 - pow(H.Length(), 2.0) / (G * M * Planet->a)), 0.5);

	// ============================

	double A_GR, B_GR;
	A_GR = Return_AGR(R,V,Planet->Mass,Sun->Mass,GRCorrections,c);
	B_GR = Return_BGR(R,V,Planet->Mass,Sun->Mass,GRCorrections,c);	

	// ============================

	//Acceleration - Pt. 1
	//GR Correction
	double GRFact = -G * M / pow(R.Length(),2.0);
	Vect tempV = VectAdd(VectScale(N, A_GR), V, B_GR);
 	Vect A = VectScale(tempV, GRFact);

	//Velocity - Pt. 1
	NewVX = Planet->v_x + (A.x * TimePeriod / 2.0);
	NewVY = Planet->v_y + (A.y * TimePeriod / 2.0);
	NewVZ = Planet->v_z + (A.z * TimePeriod / 2.0);

	//Position
	Planet->x = Planet->x + (NewVX * TimePeriod);
	Planet->y = Planet->y + (NewVY * TimePeriod);
	Planet->z = Planet->z + (NewVZ * TimePeriod);

	//Acceleration - Pt. 2
	//GR Correction
	R = Vect(Planet->x, Planet->y, Planet->z);
	V = Vect(Planet->v_x, Planet->v_y, Planet->v_z);
	N = VectScale(R, 1.0/R.Length());
	GRFact = -G * M / pow(R.Length(),2.0);
	tempV = VectAdd(VectScale(N, A_GR), V, B_GR);
 	A = VectScale(tempV, GRFact);
	
	//Velocity - Pt. 2
	Planet->v_x = NewVX + (A.x * TimePeriod / 2.0);
	Planet->v_y = NewVY + (A.y * TimePeriod / 2.0);
	Planet->v_z = NewVZ + (A.z * TimePeriod / 2.0);

	if((Planet->y < 0.0) and (tempy >= 0.0)) Planet->No_Orbits++;
}


// Outputs Planetry data to a file that is passed to function
void OutputData(Status Planet, FILE* OutputFile){

	fprintf(OutputFile, "%f, ", Planet.TimeElapsed);
	fprintf(OutputFile, "%f, ", Planet.x);
	fprintf(OutputFile, "%f, ", Planet.y);
	fprintf(OutputFile, "%f, ", Planet.z);
	fprintf(OutputFile, "%f, ", Planet.v_x);
	fprintf(OutputFile, "%f, ", Planet.v_y);
	fprintf(OutputFile, "%f " , Planet.v_z);
	fprintf(OutputFile, "\n");

}

// Calculate Eccentricity Vector & output as a "vector" 
Vect CalculateEccVector(Status Planet, Status Sun){
	Vect R = Vect(Planet.x, Planet.y, Planet.z);
	Vect V = Vect(Planet.v_x, Planet.v_y, Planet.v_z);

	Vect H = CrossProduct(R, V);

	double factor = G * (Sun.Mass + Planet.Mass);
	Vect Temp = VectScale(CrossProduct(V,H), factor);

	factor = - 1.0 / R.Length();
	Vect Result = VectAdd(Temp, R, factor);
	
	return Result;
}

// Return the Semi-Major axis
SemiEcc ReturnSemiMajor_Ecc(Status Planet, Status Sun){
	double M = Sun.Mass + Planet.Mass;
	Vect R = Vect(Planet.x, Planet.y, Planet.z);
	Vect V = Vect(Planet.v_x, Planet.v_y, Planet.v_z);
	Vect H = CrossProduct(R, V);
	
	double temp = (2.0 / R.Length()) - (pow(V.Length(),2.0) / (G * M));
	double A = 1 / temp;

	// Return the eccentricisty
	
	temp = 1 - (pow(H.Length(), 2.0) / (G * M * A));
	double E = pow(temp,0.5);
	SemiEcc Res = SemiEcc(A,E);

	return Res;
}

//Compute Mass Units
double ComputeMassUnits(double InputMass){
	double mpc = 3.08568 * pow(10.0,13.0); // m
	double yrs = 365.25 * 24.0 * 60.0 * 60.0; // s
	double G_ms = 6.67430 *pow(10.0,-11.0); // m3 kg–1 s–2

	double M_Units = pow(mpc,3.0) / (G_ms * pow(yrs,2.0)); // in kg
	double M_Convert = 1.98855 * pow(10.0,30.0) / M_Units;
	double res = InputMass * M_Convert;
	return res;
}

// Compute Velocity Units
double ComputeVelUnits(double InputVel){
	double mpc = 3.08568 * pow(10.0,13.0); // m
	double yrs = 365.25 * 24.0 * 60.0 * 60.0; // s
	double G_ms = 6.67430 *pow(10.0,-11.0); // m3 kg–1 s–2

	double M_Units = pow(mpc,3.0) / (G_ms * pow(yrs,2.0)); // in kg
	double V_Convert = pow((G_ms * M_Units / mpc), 0.5);
	double res = InputVel / V_Convert;
	return res;
}

// Compute Initial conditions
Status ComputeInitConditions(double M1, double MS, double a, double e, double vfact){
	double Ap = a * (1.0 + e);
	double Pe = a * (1.0 - e);
	double Temp = (G * (M1 + MS)) / a;
	double Vp = pow(Temp * ((1 + e) / (1 - e)), 0.5);
	double Va = pow(Temp * ((1 - e) / (1 + e)), 0.5);

	Status Planet = Status(Ap * vfact, 0.0, 0.0, 0.0, Va * vfact, 0.0, M1, 0.0);	

	return Planet;
}

// Version 2

OrbitalParameters IntOrbits2(Status BH1, Status BH2, double c, int* GRCorrections, double TimeLength, double OP, double AP, char* Filename){

	// Open & Name Files
	FILE *OutputFile_BH1;
	char title[32];
	strcpy(title,Filename);
	strcat(title, "_BH1Orbit.txt");
	OutputFile_BH1 = fopen(title, "w+");

	FILE *OutputFile_BH2;
	char title2[32];
	strcpy(title2,Filename);
	strcat(title2, "_BH2Orbit.txt");
	OutputFile_BH2 = fopen(title2, "w+");

	FILE *OutputFile_EA;
	OutputFile_EA = fopen("binary.circ.orb", "w+");

	printf("Performing integration...\n");
	printf("\n");

	//Get Parameters
	Vect Ecc1 = CalculateEccVector(BH2, BH1);
	SemiEcc SemEcc1 = ReturnSemiMajor_Ecc(BH2, BH1);

	// == DO THE LOOP ==
	// =================
	for(int t = 0; t < TimeLength; t++) {
		//BH1.Print();
		//MovePlanets(&BH2, &BH1, GRCorrections, c, AP);
		//BH1.Print();
		OutputData(BH1, OutputFile_BH1);
		MovePlanets(&BH1, &BH2, GRCorrections, c, AP);
		//BH2.Print();
		fprintf(OutputFile_EA, "%.15f, %.15f, %.15f\n", BH2.TimeElapsed, BH2.e, BH2.a);
		OutputData(BH2, OutputFile_BH2);
		if(t % 10000 == 0.0){
			printf("Time Elapsed = %f / %f\r", BH2.TimeElapsed, OP);
			fflush(stdout);
		}
		if(BH2.TimeElapsed >= OP){
			printf("Orbital period reached! T = %f\n", BH2.TimeElapsed);
			break;
		}
	}  
	// =================
	// =================

	printf("\n");
	printf("Saving data to '%s' and '%s'.\n", title, title2);
	printf("\n");

	//Close Files
	fclose(OutputFile_BH1);
	fclose(OutputFile_BH2);
	fclose(OutputFile_EA);

	//Get Parameters
	Vect Ecc2 = CalculateEccVector(BH2, BH1);
	SemiEcc SemEcc2 = ReturnSemiMajor_Ecc(BH2, BH1);

	//Return Orbital Parameters
	OrbitalParameters Param = OrbitalParameters(Ecc1, Ecc2, SemEcc1, SemEcc2);
	return Param;
}

// Semi-Major theory through Peter's equations
void SemiMajorTheory(double c, double M1, double M2, double a_0, double lim){

	double Beta = (pow(G, 3.0) / pow(c, 5.0)) * M1 * M2 * (M1 + M2);
	double A;	

	FILE* File;
	File = fopen("peters.orb", "w+");

	for(int ti=0; ti < lim; ti++){ 
		A = pow(pow(a_0,4.0) - (256 * Beta * ti / 5.0), 0.25);
		fprintf(File, "%f, %f\n", (ti / 100.0), A);
	}

	fclose(File);

}


//== PROGRAM ==

int main() {

	// ============= INTRO =============
	printf("\n");
	printf("General Relativity - Coursework\n");
	printf("===============================\n");
	printf("\n");
	printf("Submitted by: Daniel Foulds-Holt\n");
	printf("         URN: 6424523\n");
	printf("        Date: 26/05/20\n");
	printf("\n");
	printf("This code is written to perform simulations of the two body\n");
	printf("problem using a Leapfrog integrator with variable time-step.\n");
	printf("It also can apply Post-Newtonian corrections of the order\n");
	printf("PN1, PN2, and PN2.5 .\n");
	printf("\n");
	printf("===============================\n");
	printf("\n");
	// =================================
	

	// ============= Q. I  =============
	//Reframe BH mass in new units [M]:
	double M1 = ComputeMassUnits(7 * pow(10,6)); //[M]
	double M2 = ComputeMassUnits(7 * pow(10,3)); //[M]
	double a = 0.15; //mpc 
	double e = 0.0;
	double M = M1 + M2;
	double mu = M1 * M2 / M;

	//Print to screen
	printf("The binary black hole system of two black holes:\n");
	printf("\n");
	printf("M1 = %f [M]\n", M1);
	printf("M2 = %f [M]\n", M2);
	printf("Eccentricity    = %f\n", e);
	printf("Semi-major axis = %f mpc\n", a);
	printf("\n");
	printf("The initial conditions have been saved to the file \n");
	printf("'binary.ecc.init'.\n");
	printf("\n");
	printf("===============================\n");
	printf("\n");
	
	//Get Initial conditions
	Status BH2 = ComputeInitConditions(mu, M, a, e, 1.0);
	Status BH1 = Status(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, M, 0.0);
	
	//Correct positions
	MoveToCOM(&BH1, &BH2);

	//BH2.Print();	

	//Output IC to file
	OutputIC(BH1, BH2,(char*)"binary.circ.init");
	// =================================


	// ============= Q. J  =============
	double c = ComputeVelUnits(299792458);
	double T_p = 600.0;
	int GRCorrections[] = {1, 0, 0, 1};
	double TimeLength = 100000000;
	double AP = 0.001;
	//AP = 0.00001; // For accuracy

	OrbitalParameters Params = IntOrbits2(BH1, BH2, c, GRCorrections, TimeLength, T_p, AP, (char*)"Meaningless2");
	// =================================


	// ============= Q. K  =============

	SemiMajorTheory(c, M1, M2, a, 60000);
	printf("Saved all the a data to.... \n");

	// =================================

}












