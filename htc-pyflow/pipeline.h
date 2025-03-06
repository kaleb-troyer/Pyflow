#pragma once

#ifndef PIPELINE
#define PIPELINE

struct results {
	double Reynolds; 
	double FrictionFactor; 
	double Nusselt; 
	double HeatTransferCoefficient; 
	double PressureDrop; 

	results() {
		Reynolds = Nusselt = FrictionFactor
			= HeatTransferCoefficient = PressureDrop = 0.0;
	}; 
};

struct dimensions {
	double hydraulic_diameter; 
	double major_diameter; 
	double minor_diameter; 
	double height; 
	double width; 
	double length;
	double relRough; 

	dimensions() {
		hydraulic_diameter = major_diameter = minor_diameter
			= height = width = length = relRough = 0.0;
	}; 
};

struct fluid {
	double prandtl; 
	double density; 
	double viscosity; 
	double conductivity; 
	double velocity; 

	fluid() {
		prandtl = density = viscosity
			= conductivity = velocity = 0.0; 
	}; 
};

class FlowCalculator {
private:
	// Warning Flags enumeration
	enum WarningFlags {
		WARN_A = 1 << 0, // invalid geometry code
		WARN_B = 1 << 1, // Reynolds number too small or negative (< 0.001)
		WARN_C = 1 << 2, // Reynolds number too large (> 5E6)
		WARN_D = 1 << 3, // Prandtl number outside of valid range
		WARN_E = 1 << 4, // L/D ratio must be greater than 1 
		WARN_F = 1 << 5, // relRough must be < 0.05 and positive 
		WARN_G = 1 << 6, // the Graetz number is outside the acceptable range
	};

	// For bitmask tracking every warning raised
	int warnings; 
	
	// Data structures
	results s_results;
	dimensions s_dims;
	fluid s_fluid;	

	// Private core functions
	void pipe_laminar_flow(double* dims);
	void duct_laminar_flow(double* dims);
	void annulus_laminar_flow(double* dims);
	void turbulent_flow(double* dims);
	void laminar_flow(int code, double* dims); 
	void transitional_flow(int code, double* dims); 
	void pressure_drop(); 

public:
	FlowCalculator(double Ac, double Dh, double roughness, double velocity,
		double Prandtl, double density, double viscosity, double conductivity);
	~FlowCalculator() {};
	
	bool updater(double Ac, double Dh, double roughness, double velocity,
		double Prandtl, double density, double viscosity, double conductivity);

	bool calculate(unsigned int code, double* dims);

	// Accessor functions to get results
	double getReynolds() const;
	double getNusselt() const;
	double getFrictionFactor() const;
	double getHeatTransferCoefficient() const;
	double getPressureDrop() const; 

	// Warning processing functions
	void setWarning(WarningFlags flag) {
		warnings |= flag; 
	} 

	bool isWarning(WarningFlags flag) const {
		return (warnings & flag) != 0; 
	}

	int getWarnings() const {
		return warnings; 
	}
};

#endif // PIPELINE


