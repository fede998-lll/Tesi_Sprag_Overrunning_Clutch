# Sprag Overrunning Clutch: Torsional Dynamics Simulation

This project is part of a thesis focusing on the torsional dynamic behavior of a mechanical system featuring a sprag overrunning clutch. The MATLAB code provided simulates the system's dynamic behavior, calculates key parameters, and visualizes results for analysis.

## Project Overview

The main goal of this project is to study the dynamic response of a system that incorporates a sprag overrunning clutch under various conditions. The simulation accounts for:
- Angular positions and velocities of key components (input/output shafts, inner/outer races, and sprags).
- Damping coefficients and stiffness of the system.
- Geometric deformations and forces in the clutch mechanism.
- Freewheeling, engaged states, and transitions between these states of the sprag clutch.

## Code Description

The repository includes the following MATLAB scripts:

### Core Files
1. **`main.m`**:
   - Entry point for the simulation.
   - Sets up initial conditions, system parameters, and time-stepping for the simulation.
   - Implements a Runge-Kutta 4th-order integration loop to compute the system's response over time.

2. **`system_eq.m`**:
   - Defines the system of equations governing the dynamics of the shafts, races, and clutch.

3. **`clutch.m`**:
   - Implements the behavior of the sprag overrunning clutch in different states (freewheeling, engagement, fully engaged).
   - Computes forces, angular accelerations, and transitions between states.

4. **`calculate_damping.m`**:
   - Calculates the damping coefficients for the system based on stiffness, inertia, and a specified damping ratio.

5. **`computing_NBI.m`**:
   - Computes the normal force (`NBI`) on the inner race and the geometric parameters of the clutch.

6. **`objective_function.m`**:
   - Evaluates the relative angular displacement between the inner and outer races based on normal force and geometric properties.

### Visualization
1. **`plotter.m`**:
   - Generates plots for angular positions, velocities, and the behavior of the sprags over time.

