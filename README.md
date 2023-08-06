# MSAT
MSAT aims to provide a helpful tool to quantitatively analyze the shock instability problem for the schemes with three-point stencils. Detailed information about the matrix stability analysis method can be found in the following: 
 - [M. Dumbser, J.-M. Moschetta, J. Gressier, A matrix stability analysis of the carbuncle phenomenon, Journal of Computational Physics 197 (2) (2004) 647–670](https://doi.org/10.1016/j.jcp.2003.12.013)
 - [W.Ren,W.Xie,Y.Zhang,H.Yu,Z.Tian,Numerical stability analysis of shock- capturing methods for strong shocks i:second-order muscl schemes (2023)](https://arxiv.org/abs/2305.03281)

## Authors
- Weijie Ren <rwj@nudt.edu.cn>
- Wenjia Xie <xiewenjia@nudt.edu.cn>
- Ye Zhang <zhangye17a@nudt.edu.cn>
- Hang Yu <yuhang18@gfkd.edu.cn>
- Zhengyu Tian <tianzhengyu_kd@163.com>

## Usage
Following the next step to use MSAT:
 1. Clone the directory with <https://github.com/JameRwj/MSAT_SoftwareX.git>
 2. Change the settings in *Settings.dat*
 3. Generate the computational grid and export to *Grid.dat*
 4. Execute *Main Program* to run the software
 
Note that MSAT has windows version and Linux version. In Windows, this software is tested with Visual Studio 2022 and Intel oneAPI 2023. And in Linux, MSAT is tested on Ubuntu system with gfortran 9.4.0, and lapack 3.8.0.

## Settings.dat
The main settings of the MSAT are stored in *Settings.dat*, including:

- The reconstruction method used in MAST:
		1. MUSCL;
		2. ROUND.

 - The limiter function used in the MUSCL approach. The correspondence between the number and the limiter function is:
		 1. Superbee limiter;
		 2. van Leer limiter;
		 3. van Albada limiter;
		 4. minmod limiter;
		 5. the limiter proposed by Xi Deng in  [X.Deng, A unificd framework for non-linear reconstruction schemes in a compact stencil. Part 1: Beyond second order, Journal of Computational Physics 481 (2023) 112052](https://doi.org/10.1016/j.jcp.2023.112052) ;
	And it will be first-order accurate if choose 0.
	
 - The Riemann solver used in MSAT. The correspondence between the number and the Riemann solver is:
		 1. Roe solver;
		 2. HLLC solver;
		 3. HLL solver;
		 4. van Leer solver;
		 5. AUSM+ solver;
		 6. SLAU solver;
		 7. HLLE solver;
		 8. HLLEM solver.
- Test case:
		1. 2D normal shock
		2. other test cases
	
 - Numerical shock structure $\varepsilon$ of the 2D normal shock, which is between 0 and 1.
 - Mach number of the 2D normal shock.
 - The method of Initializing the 2D normal shock. If it is 1, the 2D flow field is initialized by the Rankine-Hugoniot conditions. And the 2D flow field will be initialized by projecting the steady flow field from 1D computation onto the 2D domain.
 - The iteration steps of the 1D computation.

## Grid.dat
The coordinates of the grid nodes are stored in *Grid.dat*. Note that the origin of the coordinates is the lower-left corner of the grid. There are two parts in *Grid.dat*. The first part specifies the number of grid nodes in the  x  and  y  dimensions. The second part is divided into three columns: the first column contains the  x  coordinates of each grid point, the second column contains the  y  coordinates, and the third column displays the  z  values, which should be 0 due to the two-dimensional nature of the computation.

## InitialFlow
This folder contains the initial flow field, including *InitialFlow_rho.dat*, *InitialFlow_u.dat*, *InitialFlow_v.dat*, and *InitialFlow_p.dat*. Those files will be used to initialize the flow field if other test cases are used. Note that other programs should obtain the four files.

## Results
The results of MSAT will be stored in the file of *Results*, including the flow field, result and residual of the 1D computation (if have), scatters of all eigenvalues, and the eigenvectors of the most unstable eigenvalue if it exceeds 0.

Note that the preceding introduction is mainly about the Windows version, and the Linux version has the similar usage.
