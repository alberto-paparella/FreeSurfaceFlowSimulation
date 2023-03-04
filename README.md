# Free Surface Flow Simulation on Structured Grids

Aim of this project is the simulation of a fluid dynamic phenomena, that is the free surface flow on structured grids. Thanks to this simulation, we are capable of observe the variation over the depth of the water in each point of the considered domain during time, starting from a given condition, as well as the evolution of speed in each point among both axis.

The simulation can be executed making exploiting any number of processors given that the computation matrix is divisible for such number. This parallelisation allows to reduce the computation time, as proven by the speed-up analysis, efficiency analysis and Kuck function analysis at the varying of the number of used processors.

## Conditions

- The domain is the square: $\Omega=[-0.5;0.5]\times[-0.5;0.5]$
- The initial condition is pressure $\eta=1+e^{ -\frac{1}{2s^{2}}(x^{2}+y^{2})}$ and $s=0.1$
- Speeds $u(i+\frac{1}{2},j)$ and $v(i,j+\frac{1}{2})$ are defined at the cell interfaces of the grid, while pressure $\eta$ is defined at the baricenters of the cells
- $\Delta t$ represents the temporal step
  
## Equations and numeric method

Schema is based on a semi-implicit temporal discretization and a finite difference spatial approximation.

Total water depth at the interface, also defined as draught, is expressed as:
- $H_{i+\frac{1}{2},j}^{n}=\max(0,h_{i+\frac{1}{2},j}^{n}+\eta_{i,j}^{n};h_{i+\frac{1}{2},j}^{n}+\eta_{i+1,j}^{n})$

- $H_{i,j+\frac{1}{2}}^{n}\max(0,h_{i,j+\frac{1}{2}}^{n}+\eta_{i,j}^{n};h_{i,j+\frac{1}{2}}^{n}+\eta_{i,j+1})$

Convective terms $Fu_{i+\frac{1}{2},j}^{n}$ and $Fv_{i,j+\frac{1}{2}}^{n}$ have been evaluated using an explicit upwind method.

The following finite differences discretization has been considered for the equations governing pressure $\eta$ (simulation is evaluated on a flat seabed, hence $\eta$ is equal to the height of the fluid) and speeds $u$ and $v$ on axis $x$ and $y$, respectively:
- $u_{i+\frac{1}{2};j}^{n+1}=Fu_{i+\frac{1}{2};j}^{n}-g*\frac{\Delta t}{\Delta x}*(\eta_{i+1;j}^{n+1}-\eta_{i;j}^{n+1})$
- $v_{i;j+\frac{1}{2}}^{n+1}=Fv_{i;j+\frac{1}{2}} ^{n}-g*\frac{\Delta t}{\Delta y}*(\eta_{i;j+1}^{n+1}-\eta_{i;j}^{n+1})$
- $\eta_{i;j}^{n+1}=\eta_{i;j}^{n}\frac{\Delta t}{\Delta x}*(H_{i+\frac{1}{2};j}^nu_{i+\frac{1}{2};j}^{n+1}-H_{i-\frac{1}{2};j}^{n}u_{i-\frac{1}{2};j}^{n})$
- $\frac{\Delta t}{\Delta y}*( H_{i;j+\frac{1}{2}}^{n}v_{i;j+\frac{1}{2}}^{n+1}-H_{i;j-\frac{1}{2}}^{n}v_{i;j-\frac{1}{2}}^{n})$

## Linear system

The previous equations can be evaluated using the following pentadiagonal linear system, which can be solved using the conjugate gradient method:

- $-g\frac{\Delta t^{2}}{\Delta x^{2}}H_{i-\frac{1}{2},j}^{n}\eta_{i-1,j}^{n+1}-g\frac{\Delta t^{2}}{\Delta y^{2}}H_{i,j-\frac{1}{2}}^{n}\eta_{i,j-1}^{n+1}+[1+g\frac{\Delta t^{2}}{\Delta x^{2}}(H_{i+\frac{1}{2},j}^{n}+H_{i-\frac{1}{2},j}^{n})+g\frac{\Delta t^{2}}{\Delta y^{2}}(H_{i,j+\frac{1}{2}}^{n}+H_{i,j-\frac{1}{2}}^{n})]\eta_{i,j}^{n+1}-g\frac{\Delta t^{2}}{\Delta x^{2}}H_{i+\frac{1}{2},j}^{n}\eta_{i+1,j}^{n+1}-g\frac{\Delta t^{2}}{\Delta y^{2}}H_{i,j+\frac{1}{2}}^{n}\eta_{i,j+1}^{n+1} = b_{i}^{n}$
- $b_{i}^{n} = \eta_{i,j}^{n} -  \frac{\Delta t^{2}}{\Delta x^{2}}[(HFu)_{i+\frac{1}{2},j}-(HFu)_{i-\frac{1}{2},j}^{n}] -  \frac{\Delta t^{2}}{\Delta y^{2}}[(HFv)_{i,j+\frac{1}{2}}-(HFv)_{i,j+\frac{1}{2}}^{n}]$

## MPI parallelization

Code parallelization has been realized using a message passing paradigm with distributed memory through MPI directives. Hence, matrices containing $\eta$, $u$ and $v$, of dimensions **IMAX**$\times$**IMAX**, **(IMAX+1)**$\times$**JMAX** and **IMAX**$\times$**(JMAX+1)** respectively, have been distributed among processors, vertically subdividing them into **nCPU** rectangles with the addition of two rows for each submatrix serving as ghost cells and storing information about the border rows of the adjacent matrices, to the enhance simulation accuracy. This exchange is realized using **MPI_SEND** and **MPI_RECEIVE** directives.

Information among processors about the state of the whole system is synchronized at given times (e.g., at the time of writing the results) using support matrices $\eta_1$, $u_1$ and $v_1$. Hence, writing is always done by rank 0 CPU in three moments during the execution: at the initial state $t=0.00$, at the intermidiate state $t=0.05$ and at the final state $t=0.10$. This operation is realized making use of the **MPI_ALLGATHER** directive and is done in a cycle, so that the rows of matricex $\eta_1$, $u_1$ and $v_1$ can be used as buffers.

## Results plots

Results given a $120\times 120$ matrix. Simulation has been executed using 2 cpus.

$t=0.00$

![t=0.00](https://github.com/alberto-paparella/FreeSurfaceFlowSimulation/blob/main/images/test_2_cpu_120_120_matrix_0000_0000.PNG?raw=true)

$t=0.05$

![t=0.05](https://github.com/alberto-paparella/FreeSurfaceFlowSimulation/blob/main/images/test_2_cpu_120_120_matrix_0003_0000.PNG?raw=true)

$t=0.10$

![t=0.10](https://github.com/alberto-paparella/FreeSurfaceFlowSimulation/blob/main/images/test_2_cpu_120_120_matrix_0008_0000.PNG?raw=true)

# Results

The following tables show the results of some simulations regarding the computation times while changing the number of used CPUs. The average of these results will be used for the speed-up, efficiency and Kuck function analysis. The simulations have been executed on an **i5-8250u** CPU (4 cores, 8 threads, 1.6 GHZ).

Results on a $200\times 200$ grid:
| CPU | Time 1 (s) | Time 2 (s) |Time 3 (s)  |Time 4 (s)  | Avg (s)     |
| --- | ---------- | ---------- | ---------- | ---------- | ----------- |
| 1   | 12.7693350 | 13.8580796 | 13.4519254 | 13.1758675 | 13.31380188 |
| 2   | 7.8881675  | 6.8482584  | 6.6088351  | 6.6232503  | 6.992127825 |
| 4   | 4.2107446  | 4.1933532  | 4.0723427  | 4.0085329  | 4.12123888  |

Results on a $400\times 400$ grid:
| CPU | Time 1 (s)  | Time 2 (s) |Time 3 (s)   |Time 4 (s)  | Avg (s)     |
| --- | ----------- | ---------- | ----------- | ---------- | ----------- |
| 1   | 109.2646735 | 95.8739504 | 95.853.359  | 95.7647304 | 99.18909755 |
| 2   | 73.5225888  | 50.7483651 | 50.3488792  | 51.9943347 | 56.90354195 |
| 4   | 30.8085830  | 34.3894880 | 35.2220634  | 38.6311879 | 34.76283058 |

## Speed-up

Speed-up is calculated as:
- $S(p)=\dfrac{T_{(p=1)}}{T(p)}$

where $T_{(p=1)}$ is the average computation time using **1 single CPU** (serial execution) and $T(p)$ is the average computation time using **p CPUs**.

![speedup](https://github.com/alberto-paparella/FreeSurfaceFlowSimulation/blob/main/images/speedup.png?raw=true)

## Efficiency

Efficiency is calculated as:
- $E(p)=\dfrac{S(p)}{p}$

i.e., the ratio between the speed-up and the relative number of used CPUs.

![efficiency](https://github.com/alberto-paparella/FreeSurfaceFlowSimulation/blob/main/images/efficiency.png?raw=true)

## Kuck function

Kuck function is evaluated as:
- $K(p)=S(p)*E(p)$

The number of CPUs **p\*** such that $K(p)$ is at its maximum peek is the best tradeoff between speed-up and efficiency.

![kuck](https://github.com/alberto-paparella/FreeSurfaceFlowSimulation/blob/main/images/kuck.png?raw=true)
