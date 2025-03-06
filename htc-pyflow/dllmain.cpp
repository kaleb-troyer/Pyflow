// dllmain.cpp : Defines the entry point for the DLL application.
#include "pch.h"
#include "pipeline.h"

#ifdef _MSC_VER
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif

// Expose a C-style API for use from Python.
extern "C" {

    // Create a new FlowCalculator instance.
    DLL_EXPORT FlowCalculator* CreateFlowCalculator(
        double Ac, double Dh, double roughness, double velocity,
        double Prandtl, double density, double viscosity, double conductivity
    ) {
        return new FlowCalculator(Ac, Dh, roughness, velocity, Prandtl, density, viscosity, conductivity);
    }

    // Delete an existing FlowCalculator instance.
    DLL_EXPORT void DeleteFlowCalculator(FlowCalculator* calc) {
        delete calc;
    }

    // Update flow calculator instance parameters. 
    DLL_EXPORT bool updater(
        FlowCalculator* calc, double Ac, double Dh, double roughness, double velocity,
        double Prandtl, double density, double viscosity, double conductivity
    ) {
        if (!calc) return false; 
        return calc->updater(Ac, Dh, roughness, velocity, Prandtl, density, viscosity, conductivity);
    }

    // Call calculate on the FlowCalculator instance. The dims array is passed
    // directly to the method that matches the chosen flow code.
    DLL_EXPORT bool core(FlowCalculator* calc, unsigned int code, double* dims) {
        if (!calc) return false;
        return calc->calculate(code, dims);
    }

    // Getter functions to retrieve results.
    DLL_EXPORT double getReynolds(FlowCalculator* calc) {
        return calc ? calc->getReynolds() : 0.0;
    }

    DLL_EXPORT double getNusselt(FlowCalculator* calc) {
        return calc ? calc->getNusselt() : 0.0;
    }

    DLL_EXPORT double getFrictionFactor(FlowCalculator* calc) {
        return calc ? calc->getFrictionFactor() : 0.0;
    }

    DLL_EXPORT double getHeatTransferCoefficient(FlowCalculator* calc) {
        return calc ? calc->getHeatTransferCoefficient() : 0.0;
    }

    DLL_EXPORT double getPressureDrop(FlowCalculator* calc) {
        return calc ? calc->getPressureDrop() : 0.0;
    }

    DLL_EXPORT int getWarnings(FlowCalculator* calc) {
        return calc ? calc->getWarnings() : 0;
    }

} // extern "C"

// Standard DLL entry point.
BOOL APIENTRY DllMain( 
    HMODULE hModule,
    DWORD  ul_reason_for_call,
    LPVOID lpReserved
) {

    switch (ul_reason_for_call) {
        case DLL_PROCESS_ATTACH:
        case DLL_THREAD_ATTACH:
        case DLL_THREAD_DETACH:
        case DLL_PROCESS_DETACH:
            break;
    }

    return TRUE;
}

