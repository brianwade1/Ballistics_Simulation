# 6-DOF Ballistics Simulation

> MATLAB Simulation of a six degree of freedom (6-DOF) ballistic flight of a 120mm mortar round.

![Sample Output](/images/output.png)

![Sample Output](/images/output_console.png)

## Files

The simulation is setup, run, and controlled from the Mortar_Sim.m file. This main file uses the ODE45 solver in MATLAB with the equations of motion (EoM) specified for each time step in the EoM.m file. The file reads input data from a standard atmosphere table and a table of aerodynamic coefficients vs. Mach number (McCoy, 1998). The 6-DOF model was developed based on the methodology in McCoy, chapter 9.

1. [Mortar_Sim.m](Mortar_Sim.m) - Main program
2. [EoM.m](EoM.m) - Equations of Motion
3. [std_atm.csv](std_atm.csv) - Table of the standard atmosphere
4. [Aerodynamic_Char_120mm_Mortar.xlsx](Aerodynamic_Char_120mm_Mortar.xlsx) - Aerodynamic coefficients of the 120mm for different Mach number. Data extracted from McCoy, 1998 on page 220.

## Setup and Run

The simulation is setup from within the main program ([Mortar_Sim.m](Mortar_Sim.m)). The simulation parameters are set in lines 68-81. These lines set the initial conditions of the mortar round as it exits the mortar tube and enters free flight.

![initial_conditions](/images/initial_conditions.png)

The first variable (Vo_set) is the initial velocity at the barrel exit in meters per second. The next variable (phi_0_set) is the vertical angle of departure with respect to the inertial frame in degrees (this is the angle of the mortar tub measured from the horizontal ground plane). The third variable (theta_0_set) is the horizontal angle of departure in degrees (this is the angle of the mortar tub measured from a vertical plan perpendicular to the ground and along the x-axis).

The next two variables (w_z_0_set and w_y0_set) are the initial pitch and yaw of the mortar round as it exists the barrel both in radians per second. The sign convention is positive for a nose-up or left-yaw rate.

The third set of variables (alpha_0_set and beta_0_set) are the pitch and yaw angles of the round at the barrel exit. These are measures of the angular difference between the mortar body axis and the inertial (gun tube) axis. Both angles are in degrees. These use the same sign convention as above where a nose-up or left-yaw are positive. Note that these angles establish an initial off-axis flight of the round.

The final set of variables (x_0, y_0, and z_0) are the x, y, and z locations of the end of the mortar barrel (where the round enters free flight). These are all measured in meters.

## MATLAB version and add-ons

This simulation was built using MATLAB 2019b. No other MATLAB add ons should be needed to run the simulation.

## References

1. McCoy, RL, Modern Exterior Ballistics: The Launch and Flight Dynamics of Symmetric Projectiles, Schiffer Military History, Atglen, PA, 1998.
