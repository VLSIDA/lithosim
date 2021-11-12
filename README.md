# lithosim

This is a very basic lithography simulation and pixel-based OPC tool.

## Simulation

The simulation uses an analytical model similar to [A. Poonawala,
P. Milanfar, “A Pixel-Based Regularization Approach to Inverse
Lithography”,Microelectronic Engineering, 84 (2007)
pp. 2837–2852](https://users.soe.ucsc.edu/~milanfar/publications/journal/Microelectronic_Final.pdf).

## OPC

The OPC just does a simulated annealing algorithm to minimize error
between the target mask and the simulated mask. It does not convert
the pixel-based mask into a manufacturable mask.

## Example Results

<img src="examples/tiny-mask-90nm.jpg" alt="Mask (target)" width="100"/>
<img src="examples/tiny-aerial-90nm.jpg" alt="Aerial Image" width="100"/>
<img src="examples/tiny-contours-90nm.jpg" alt="Contours" width="100"/>

<img src="examples/tiny-opc-90nm.jpg" alt="OPC Mask" width="100"/>
<img src="examples/tiny-opc-aerial-90nm.jpg" alt="OPC Aerial Image" width="100"/>
<img src="examples/tiny-opc-contours-90nm.jpg" alt="OPC Contours" width="100"/>
