# Super-Fast Mesh-Free 2D Transient Heat Conduction Simulation in Circular Plate

This MATLAB code simulates transient heat conduction in circular functionally graded material (FGM) plates using the Generalized Differential Quadrature (GDQ) and Crank-Nicolson (CN) methods.

## Purpose
The script models temperature distribution over time within a 2D circular FGM plate. It incorporates material properties that vary with temperature and can handle partial thermal loading, an essential feature for real-life applications.

## Key Features
- Super-Fast
- Customizable Boundary Conditions: The code supports boundary conditions such as insulation, convection, prescribed temperature, and heat flux on any edge within the domain. Each condition can be applied partially across boundaries to accurately simulate real-world thermal loading.
- Variable Material Properties: Thermal conductivity, specific heat, Young’s modulus, and other material properties change based on material type and temperature, allowing for a more accurate FGM representation.

## Usage
To run the simulation, specify the following:

- Mesh details: Number of nodes along the radius (`nr`) and thickness (`nz`).
- Time details: Time step (`dt`) and simulation duration.
- Material and geometric properties: Parameters such as a, b, hh, and kisi.
- Boundary conditions: Set `right_edge`, `top_edge`, `bottom_edge`, `left_edge` for each boundary with options for full or partial application using `partial_ratio_top` and `partial_ratio_bottom`.

## Inputs and Outputs
- **Inputs**: `nr`, `nz`, `dt`, `a`, `b`, `hh`, `kisi`, `right_edge`, `top_edge`, `bottom_edge`, `left_edge`, `partial_ratio_top`, `partial_ratio_bottom`
- **Outputs**: `TT` (temperature matrix), `tt` (time vector), and material-specific outputs.
