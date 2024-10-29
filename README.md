# 2D Heat Conduction Simulation

This MATLAB code simulates transient heat conduction in functionally graded materials circualr plate using the Generalized Differential Quadrature (GDQ) and Crank-Nicolson (CN) methods.

## Purpose
The script models temperature changes over time in 2D circular geometry, with customizable boundary conditions and thermal properties that vary by material type and temperature.

## Usage
Run the code by specifying the mesh (`nr`, `nz`), time step (`dt`), and material properties (e.g., `a`, `b`, `hh`, `kisi`). Set boundary conditions for each edge.

## Inputs and Outputs
- **Inputs**: `nr`, `nz`, `dt`, `a`, `b`, `hh`, `kisi`, `right_edge`, `top_edge`, `bottom_edge`, `left_edge`, `partial_ratio_top`, `partial_ratio_bottom`
- **Outputs**: `TT` (temperature matrix), `tt` (time vector), and material-specific outputs.
