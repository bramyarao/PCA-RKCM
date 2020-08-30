# Reduced Order Modeling using Principal Component Analysis
## This work is copyrighted.

### Introduction
<div style="text-align: justify"> 
In this work model order reduction is proposed for the strong form <a href="https://doi.org/10.1002/num.20539" target="blank">Reproducing Kernel Collocation Method</a> (RKCM), which is a meshfree numerical method based on points as discretization in the domain.
<br/>
<br/>
The main aim of model order reduction is to find a lower dimension approximation of a full model solution by projecting it onto a lower dimensional space. Proper Orthogonal Decomposition (POD), also known as Principal Component Analysis (PCA), is one of the most popular methods used to construct this projection operator. In this method the proper orthogonal modes of a system of equations are constructed and then these modes are truncated as required to construct the lower dimensional approximation. Model reduction for RKCM is done using a Least Squares Galerkin projection to project the over-determined system of equations. For the theory, full formulation and details of examples given below, refer to <a href="LatexFormulation/main.pdf" target="blank">this pdf</a>. The coding is done in MATLAB. 
</div>

<br/>
### Numerical Examples
<div style="text-align: justify"> 
Examples from elasto-dynamics are given below which show the comparison between the full model and the reduced model using 'r' number of POD modes.
</div>

### Example 1: 1D Wave equation
<div style="text-align: justify"> 
In this simple example a single POD is enough to capture the deformation of the full scale solution.
</div>

<p align="center">
<video width="500" controls preload>
  <source src="videos/Wave.mp4">
</video>
</p>

### Example 2: 1D Bi-material bar
<div style="text-align: justify">
In this example, the bar is made up of two materials in either half, which have different modulus of elasticity. More POD modes are required in order to obtain a reduced solution which captures the deformation of the full scale solution.
</div>

<p align="center">
<video width="500" controls preload>
  <source src="videos/Bimaterial.mp4">
</video>
</p>
