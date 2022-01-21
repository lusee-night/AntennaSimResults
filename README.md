# AntennaSimResults

Sweep Parameter (Discrete Simulation Results):
Frequency 20 - 50 MHz, increments of 1 MHz;
Antenna Angle 15 - 60 deg, increments of 15 deg;

## Constant Parameters:
Antenna Length 6 m; 
Antenna Radius 20 cm; 
Distance of Antenna to center (offset) 25 cm;
LuSEE Box 1x1x1 m^3;
Lander Radius 1.5 m; 
Lander Thickness 0.5 m; 
Lander Height above Lunar Surface 1 m; 
Lunar Surface Roughness 10 cm;

## Further Sim Details: 
Excitation Impedance 50 ohm;
Freq_low 10 MHz (sim box); Freq_high 50 MHz (mesh size of radiation boundary); Mesh Factor 10;
Solution Frequency 100 MHz (mesh size of interior objects);
Far Field Radiation Boundary for all sides of simulation volume except -z; Lunar Soil simulated via an Impedance Layer.
