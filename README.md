# Heat2d Stencil DOALL parallelism
 This project aim to compare the serial and parallel execution of Heat2d stencil algorithm. 
# Problem Statement
 Write Heat2D stencil algorithm that utilizes parallel programming DOALL parallelism in OpenMP in C. Compare the performance with sequential ordering. 
# How to run Code
 ```bash run_script.sh ```
 
The execution command generates an ‘output.txt’ file, which records the execution times for both serial and parallel computations, along with their respective parameters.
Upon convergence, the resulting matrices are stored in ‘ser_<Num of Itr>.txt’ for serial execution and ‘par_<Num of Iterations>.txt’ for parallel execution.
Additionally, the final heat maps are visualized using GNU plot files, named ‘serial_<Num of Itr>.gp’ and ‘parallel_<Num of Itr>.gp’ for serial and parallel computations respectively. 
# Visualisation
``` gnuplot --persist <filename>.gp  ```

Use this command to Visualize the gnu plot maps. Replace <filename> with the actual filename of gnu plot heatmap.

<img src="https://github.com/pavansai444/Heat2d-Stencil-DOALL-parallelism/blob/main/Visualisation_img.png" alt="Visualisation" width=60% height=auto>

# Outputs
| N        | Sequential time | Parallel time |
|----------|-----------------|---------------|
| 5        | 0.000020        | 0.004029      |
| 10       | 0.000125        | 0.001270      |
| 50       | 0.022399        | 0.060640      |
| 100      | 0.202321        | 0.534386      |
| 500      | 40.586021       | 35.525296     |
| 1000     | 164.265472      | 157.228200    |
| 5000     | 4021.034153     | 3312.206442   |
| 10000    | text1           | text2         |
