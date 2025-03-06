
#include "pch.h"
#include "pipeline.h"
#include <cmath>

FlowCalculator::FlowCalculator(
	double Ac, double Dh, double roughness, double velocity,
	double Prandtl, double density, double viscosity, double conductivity
) {
	updater(Ac, Dh, roughness, velocity,Prandtl, density, viscosity, conductivity); 
}

bool FlowCalculator::updater(
	double Ac, double Dh, double roughness, double velocity,
	double Prandtl, double density, double viscosity, double conductivity
) {

	// Reset warning flags
	warnings = 0; 

	// Set fluid properties
	s_fluid.prandtl = Prandtl;
	s_fluid.density = density;
	s_fluid.velocity = velocity;
	s_fluid.viscosity = viscosity;
	s_fluid.conductivity = conductivity;

	// Set fixed dimensions
	s_dims.hydraulic_diameter = Dh;
	s_dims.relRough = roughness;

	// Calculate Reynolds number (based on velocity and fixed properties)
	s_results.Reynolds = density * Dh * velocity / viscosity;
	return true; 
}

bool FlowCalculator::calculate(
	unsigned int code, double* dims
) {	

	if (s_results.Reynolds < 0.001) {
		setWarning(WARN_B);
	} else if (s_results.Reynolds > 5e6) {
		setWarning(WARN_C);
	}

	if (s_results.Reynolds > 3000) { // turbulent flow
		turbulent_flow(dims);
	} else if (s_results.Reynolds <= 2300) { // laminar flow
		laminar_flow(code, dims); 
	} else {
		transitional_flow(code, dims); 
	}

	pressure_drop();
	return true; 
}; 

void FlowCalculator::laminar_flow(
	int code, double* dims
) {

	switch (code) {
	case 10:
		pipe_laminar_flow(dims);
		break;
	case 20:
		duct_laminar_flow(dims);
		break;
	case 30:
		annulus_laminar_flow(dims);
		break;
	default: 
		setWarning(WARN_A); 
	}
}

void FlowCalculator::turbulent_flow(double* dims) {
	// Use dims to update dimensions for turbulent flow.
	s_dims.length = dims[2];

	// Checking for violations on parameter ranges
	if (s_fluid.prandtl < 0.004 || s_fluid.prandtl > 2000) {
		setWarning(WARN_D); 
	}

	if (s_dims.length / s_dims.hydraulic_diameter <= 1) {
		setWarning(WARN_E); 
	}

	if (s_dims.relRough > 0.05 || s_dims.relRough < 0.0) {
		setWarning(WARN_F); 
	}

	// friction factor
	double ff_fd = 0.0; 
	if (s_dims.relRough > 1.0E-5) {
		// Zigrang and Sylvester
		//double A0 = 2.00;
		//double A1 = 7.54;
		//double A2 = 5.02;
		//double A3 = 13.0;
		//double A4 = A0 * s_dims.relRough / (A1 * s_dims.hydraulic_diameter);
		//double ff_fd = std::pow(-A0 * std::log10((A4 - (A2 / s_results.Reynolds)
		//	* std::log10((A4 + (A3 / s_results.Reynolds))))), -A0);

		double A0 = -2.000;
		double A1 =  3.710;
		double A2 =  1.975;
		double A3 =  3.930;
		double A4 =  1.092;
		double A5 =  7.627;
		double A6 =  395.9;
		// Offor and Alabi, Advances in Chemical Engineering and Science, 2016, 6, 237-245
		ff_fd = std::pow(A0 * std::log10((s_dims.relRough / A1) - (A2 / s_results.Reynolds)
			* std::log(std::pow(s_dims.relRough / A3, A4) + (A5 / (s_results.Reynolds + A6)))), A0);
	} else {

		double A0 = -0.001570232;
		double A1 =  0.394203137;
		double A2 =  2.534153311;
		// Li, Seem, and Li, A New Explicity Equation for Accurate Friction Factor Calculation for Smooth Tubes, 2011
		ff_fd = std::pow(A0 / std::log(s_results.Reynolds) + A1
			/ std::pow(std::log(s_results.Reynolds), 2.0) + A2 / std::log(s_results.Reynolds), 3.0) * 4.0; 
	}

	// Nusselt number
	double N1 = (ff_fd / 8.0) * (s_results.Reynolds - 1000.0) * s_fluid.prandtl;
	double N2 = 1.0 + 12.7 * (std::pow(s_fluid.prandtl, 2.0 / 3.0) - 1.0) * std::sqrt(ff_fd / 8.0);
	double Nu_fd = N1 / N2;
	
	// Accounting for low prandtl numbers in Nusselt calculation
	if (s_fluid.prandtl < 0.5) {
		double Nu_low_prandtl = 4.8 + 0.0156 * std::pow(s_results.Reynolds, 0.85) * std::pow(s_fluid.prandtl, 0.93); 
		if (s_fluid.prandtl < 0.1) {
			Nu_fd = Nu_low_prandtl; 
		} else {
			Nu_fd = Nu_low_prandtl + (s_fluid.prandtl - 0.1) * (Nu_fd - Nu_low_prandtl) / 0.4; 
		}
	}

	// Accounting for developing flow
	s_results.FrictionFactor = ff_fd * (1.0 + std::pow(s_dims.hydraulic_diameter / s_dims.length, 0.7));
	s_results.Nusselt = Nu_fd * (1.0 + std::pow(s_dims.hydraulic_diameter / s_dims.length, 0.7));

	// Heat Transfer Coefficient
	s_results.HeatTransferCoefficient = s_results.Nusselt * s_fluid.conductivity / s_dims.hydraulic_diameter;

}

void FlowCalculator::transitional_flow(int code, double* dims) {
	
	// calculating laminar results
	laminar_flow(code, dims); 
	double Nu_la = s_results.Nusselt; 
	double ff_la = s_results.FrictionFactor; 
	double ht_la = s_results.HeatTransferCoefficient; 

	// calculating turbulent results
	turbulent_flow(dims);

	// interpolating the solution using Reynolds number
	double factor = (s_results.Reynolds - 2300) / (3000 - 2300);
	s_results.Nusselt = Nu_la + factor * (s_results.Nusselt - Nu_la); 
	s_results.FrictionFactor = ff_la + factor * (s_results.FrictionFactor - Nu_la); 
	s_results.HeatTransferCoefficient = ht_la + factor * (s_results.HeatTransferCoefficient - ht_la); 
	
}

void FlowCalculator::pipe_laminar_flow(double* dims) {
	// Use dims to update dimensions for pipe flow.
	s_dims.major_diameter = dims[0];
	s_dims.length = dims[2];

	// Non-dimensional numbers for Nusselt and friction factor calculations
	double Gz = s_dims.hydraulic_diameter * s_results.Reynolds * s_fluid.prandtl / s_dims.length; 
	double ZH = 1.0 / Gz; 
	double ZM = ZH * s_fluid.prandtl; 

	// Checking for violations on parameter ranges
	if (s_fluid.prandtl < 0.1) {
		setWarning(WARN_D); 
	}

	if (ZH < 1.0E-6) {
		setWarning(WARN_G); 
	}

	// friction factor
	double A0 = 64.0; 
	double A1 = 3.44; 
	double A2 = 1.25; 
	double A3 = 0.00021; 
	double A4 = 4.00; 
	double ff_fd = A0 / s_results.Reynolds; 
	double Lp = s_dims.length / (s_dims.hydraulic_diameter * s_results.Reynolds); 
	double F1 = A1 / std::sqrt(Lp);
	double F2 = (A2 / (A4 * Lp)) + (A0 / A4) - (A1 / std::sqrt(Lp)); 
	double F3 = 1 + (0.00021 / std::pow(Lp, 2)); 
	s_results.FrictionFactor = (A4 / s_results.Reynolds) * (F1 + (F2 / F3)); 

	// Nusselt Number
	//double Nu_fd_T = 3.66; // constant temperature boundary condition (lower bound)
	//double Nu_fd_H = 4.36; // constant heat flux boundary condition (upper bound)

	//double B0 = 0.049; 
	//double B1 = 0.020; 
	//double B2 = 1.120; 
	//double B3 = 0.065;
	//double B4 = 0.700; 
	//double N1 = (B0 + (B1 / s_fluid.prandtl)) * std::pow(Gz, B2); 
	//double N2 = 1.0 + (B3 * std::pow(Gz, B4));
	
	// Bennet, T.D., Journal of Heat Transfer, 2020
	double Nu_fd_T = (std::pow(5.001 / std::pow(ZH, 1.119) + 136.0, 0.2978) - 0.6628)
		/ (std::tanh(2.444 * std::pow(ZM, 1.0/6.0) * (1 + 0.565 * std::pow(ZM, 1.0/3.0)))); 
	double Nu_fd_H = (std::pow(6.562 / std::pow(ZH, 1.137) + 220.4, 0.2932) - 0.5003)
		/ (std::tanh(2.530 * std::pow(ZM, 1.0 / 6.0) * (1 + 0.639 * std::pow(ZM, 1.0 / 3.0))));

	/*
	The average of the upper and lower bounds on the Nusselt number is taken, 
	as opposed to letting the user select their boundary conditon. Because these
	are limits on the heat transfer coefficient, the actual condition is likely
	found somewhere in between. 
	*/
	//s_results.Nusselt = (N1 / N2) + (Nu_fd_T + Nu_fd_H) / 2.0;
	s_results.Nusselt = (Nu_fd_T + Nu_fd_H) / 2.0; 

	// Heat Transfer Coefficient
	s_results.HeatTransferCoefficient = s_results.Nusselt * s_fluid.conductivity / s_dims.hydraulic_diameter; 

};

void FlowCalculator::duct_laminar_flow(double* dims) {
	// Use dims to update dimensions for duct flow.
	s_dims.height = dims[0];
	s_dims.width  = dims[1];
	s_dims.length = dims[2];

	// Implementation would go here.
	// For now, you might set some default or computed values.
}

void FlowCalculator::annulus_laminar_flow(double* dims) {
	// Use dims to update dimensions for annulus flow.
	s_dims.major_diameter = dims[0];
	s_dims.minor_diameter = dims[1];
	s_dims.length = dims[2];

	// Implementation would go here.
	// For now, you might set some default or computed values.
}

void FlowCalculator::pressure_drop() {
	/*
	If a length is not specified, the pressure drop is not calculated
	and set as -1.0. Otherwise, the friciton factor is used to estimate
	the pressure drop across the length. 
	*/
	if (std::isinf(s_dims.length)) {
		s_results.PressureDrop = -1.0; 
	} else {
		s_results.PressureDrop = s_results.FrictionFactor * (s_dims.length / s_dims.hydraulic_diameter)
			* 0.5 * s_fluid.density * std::pow(s_fluid.velocity, 2.0);
	}
}; 

double FlowCalculator::getReynolds() const{
	return s_results.Reynolds; 
};

double FlowCalculator::getNusselt() const{
	return s_results.Nusselt; 
};

double FlowCalculator::getFrictionFactor() const{
	return s_results.FrictionFactor; 
};

double FlowCalculator::getHeatTransferCoefficient() const {
	return s_results.HeatTransferCoefficient;
};

double FlowCalculator::getPressureDrop() const {
	return s_results.PressureDrop;
};

