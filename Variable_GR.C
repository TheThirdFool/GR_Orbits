
//This is a C code for calculating the motion of the sun earth mercury system as two seperate 2-body problems.
//Written by Daniel Foulds-Holt 
//df00177@surrey.ac.uk

//== INCLUDES ==
#include <stdio.h>
#include <math.h>

//== GLOBALS ==
double A_Earth = 1.0; //AU
double E_Earth = 0.0167; //Eccentricity 

double A_Merc  = 0.387098;
double E_Merc  = 0.205635;

double M_Sun   = 1.0; //Solar Masses
double M_Earth = 3.0024584  * pow(10,-6);
double M_Merc  = 1.65956463 * pow(10,-7);

double G = 1.0; //Gravitational constant

int TimeLength = 10000;
//double AccuracyParam = 0.1; // Time interval - For varying to find energy drift
double AccuracyParam = 0.001; // For accuracy


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

//== FUNCTIONS ==

// Outputs the initial conditions to a file
//int OutputIC(InitialConditions IC, char* Filename){ \\ <--- For the old initial conditions file
int OutputIC(Status Planet, char* Filename){

	FILE *fp;
	fp = fopen(Filename, "w+");

//	fprintf(fp, "\n == Initial Conditions ==\n\n");
//	fprintf(fp, "Mass           = %f Solar Massess\n", IC.Mass);
//	fprintf(fp, "Apoapsis       = %f AU\n", IC.Ap);
//	fprintf(fp, "Periapsis      = %f AU\n", IC.Pe);
//	fprintf(fp, "Eccentricity   = %f\n", IC.Ecc);
//	fprintf(fp, "\n");
//	fprintf(fp, "Apoapsis Vel.  = %f Units\n", IC.Va);
//	fprintf(fp, "Periapsis Vel. = %f Units\n", IC.Vp);
//	fprintf(fp, "\n");

	fprintf(fp, "%f, ", Planet.x);
	fprintf(fp, "%f, ", Planet.y);
	fprintf(fp, "%f, ", Planet.z);
	fprintf(fp, "%f, ", Planet.v_x);
	fprintf(fp, "%f, ", Planet.v_y);
	fprintf(fp, "%f " , Planet.v_z);
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

// Returns the total energy of the planet
double GetEnergy(Status Planet) {
	double TotalVSq = pow(Planet.v_x, 2.0) + pow(Planet.v_y, 2.0) + pow(Planet.v_z, 2.0);
	
	double Temp = pow(Planet.x, 2.0) + pow(Planet.y, 2.0) + pow(Planet.z, 2.0);
	double R = pow(Temp, 0.5);

	double KE = 0.5 * TotalVSq;
	double PE = G * Planet.Mass / R;

	return (KE - PE);
} 

// Performs an iteration of the Euler method & updates planet positions
void MovePlanets(Status *Sun, Status *Planet) {

	double NewVX, NewVY, NewVZ;	
	double ax, ay, az;

	double M = Sun->Mass + Planet->Mass;
	double R_ns = pow(Planet->x,2.0) + pow(Planet->y, 2.0) + pow(Planet->z, 2.0); 
	double R = pow(R_ns, 1.5);

	double Temp = R / (G * M);
	double TimePeriod = AccuracyParam * pow(Temp, 0.5);
	Planet->TimeElapsed += TimePeriod;

	//Acceleration - Pt. 1
	ax = -((G * M * Planet->x) / R);
	ay = -((G * M * Planet->y) / R);
	az = -((G * M * Planet->z) / R);
 
	//Velocity - Pt. 1
	NewVX = Planet->v_x + (ax * TimePeriod / 2.0);
	NewVY = Planet->v_y + (ay * TimePeriod / 2.0);
	NewVZ = Planet->v_z + (az * TimePeriod / 2.0);

	//Position
	Planet->x = Planet->x + (NewVX * TimePeriod);
	Planet->y = Planet->y + (NewVY * TimePeriod);
	Planet->z = Planet->z + (NewVZ * TimePeriod);

	//Acceleration - Pt. 2
	ax = -((G * M * Planet->x) / R);
	ay = -((G * M * Planet->y) / R);
	az = -((G * M * Planet->z) / R);

	//Velocity - Pt. 2
	Planet->v_x = NewVX + (ax * TimePeriod / 2.0);
	Planet->v_y = NewVY + (ay * TimePeriod / 2.0);
	Planet->v_z = NewVZ + (az * TimePeriod / 2.0);

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

//Outputs the Energy data to file
void OutputEnergyData(Status Planet, FILE* OutputFile_Energy, double TotalE, double TimePeriod){

	double FractionalE;
	TotalE = TotalE - GetEnergy(Planet);
	FractionalE = TotalE / GetEnergy(Planet);
	
	printf("\tMercury II = %.15f\n", TotalE);
	//Print to file the TimePeriod and the change in Total and Fractional Energy
	fprintf(OutputFile_Energy, "%.15f, %.15f, %.15f\n", TimePeriod, TotalE, FractionalE);
	//printf("%.15f, %.15f, %.15f\n", TimePeriod, TotalE, FractionalE);
}


//== PROGRAM ==

int main() {

	//Initial conditions
	Status Sun_E = Status(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, M_Sun, 0.0);
	Status Sun_M = Status(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, M_Sun, 0.0);

	//Earth Initial Conditions ================================
	double Ap_Earth = A_Earth * (1.0 + E_Earth);
	double Pe_Earth = A_Earth * (1.0 - E_Earth);
	double Temp = (G * (M_Earth + M_Sun)) / A_Earth;
	double Vp_Earth = pow(Temp * ((1 + E_Earth) / (1 - E_Earth)), 0.5);
	double Va_Earth = pow(Temp * ((1 - E_Earth) / (1 + E_Earth)), 0.5);

	Status Earth = Status(Ap_Earth, 0.0, 0.0, 0.0, Va_Earth, 0.0, M_Earth, 0.0);	

	//Mercury Initial Conditions ==============================
	double Ap_Merc = A_Merc * (1.0 + E_Merc);
	double Pe_Merc = A_Merc * (1.0 - E_Merc);
	Temp = (G * (M_Merc + M_Sun)) / A_Merc;
	double Vp_Merc = pow(Temp * ((1 + E_Merc) / (1 - E_Merc)), 0.5);
	double Va_Merc = pow(Temp * ((1 - E_Merc) / (1 + E_Merc)), 0.5);

	Status Merc = Status(Ap_Merc, 0.0, 0.0, 0.0, Va_Merc, 0.0, M_Merc, 0.0);


	//Declare Variables =======================================
	InitialConditions IC_Earth = InitialConditions(M_Earth, Ap_Earth, Pe_Earth, E_Earth, Va_Earth, Vp_Earth);
	InitialConditions IC_Merc = InitialConditions(M_Merc, Ap_Merc, Pe_Merc, E_Merc, Va_Merc, Vp_Merc);

	DisplayIC(IC_Earth, IC_Merc);

	//This is for the old initial conditions!
//	OutputIC(IC_Earth, (char*)"Earth.ini");
//	OutputIC(IC_Merc, (char*)"Mercury.ini");

	OutputIC(Earth, (char*)"Earth.ini");
	OutputIC(Merc, (char*)"Mercury.ini");


	//Move to COM frame =======================================
	//MoveToCOM(&Sun_E, &Earth);
	//MoveToCOM(&Sun_M, &Merc);


	//Start Iterator ==========================================
	FILE *OutputFile_E;
	OutputFile_E = fopen("EarthOrbit.txt", "w+");

	FILE *OutputFile_M;
	OutputFile_M = fopen("MercuryOrbit.txt", "w+");


	// In order to vary the accuracy parameter, uncomment all the code in both 'ENERGY' sections.	
	//ENERGY ===================================================
	//FILE *OutputFile_Energy;
	//OutputFile_Energy = fopen("EnergyShift.txt", "w+");

	  //For finding the change in Energy with a changing timestamp
    //double TotalEM = 0.0;

	//for(int i = 0; i < 100; i++){
	  //Finds the total Energy in the starting positions
	//TotalEM = GetEnergy(Merc);
	//==========================================================

	// == DO THE LOOP ==
	for(int t = 0; t < TimeLength; t++) {
		// Calculate velocity and update position
		//Earth ============
		MovePlanets(&Sun_E, &Earth);
		//MoveToCOM(&Sun_E, &Earth);
		OutputData(Earth, OutputFile_E);

		//Mercury ==========
		MovePlanets(&Sun_M, &Merc);
		//MoveToCOM(&Sun_E, &Merc);
		OutputData(Merc, OutputFile_M);

		if(Merc.TimeElapsed > 1.5) break;
	
	}  

	//ENERGY ===================================================
	  //Output energy differences to file
	//OutputEnergyData(Merc, OutputFile_Energy, TotalEM, AccuracyParam);

	  //Increment TimePeriod and reset planetry positions
	//AccuracyParam = AccuracyParam - 0.001;
	//Earth = Status(Ap_Earth, 0.0, 0.0, 0.0, Va_Earth, 0.0, M_Earth, 0.0);	
	//Merc = Status(Ap_Merc, 0.0, 0.0, 0.0, Va_Merc, 0.0, M_Merc, 0.0);
	//TotalEM = 0.0;
	//}
	//fclose(OutputFile_Energy);
	//==========================================================

	fclose(OutputFile_E);
	fclose(OutputFile_M);

	printf("\n DONE!\n");

}












