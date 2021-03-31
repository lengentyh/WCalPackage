# basis_rotation_five
The program rotate the spatial coordinates of the orbitals that represent hopping parameters.

## restriction
The orbitals should NOT be hybridized, nor should them be spin polarized.   

## Getting Started
1) load a ```ftn58sparse``` into workspace and fill it into
```
[new_ftn58sparse,trot]=rotateftn58_five(ftn58sparse,a,b,c,orbitype,orbital)
```
2) Determine the Euler angles ```a,b,c``` of the rotation and the orbitals (also their type, e.g. s,p,d ).
3) The orbitals should be in a continuous series in Wannier90
4) Also, these orbitals should be ordered in the following sequence
```
p:pz,px,py
d:dz2, dxz, dyz, dx2-y2, dxy
f:fz3, fxz2, fyz2, fz(x2-y2), fxyz, fx(x2-3y2), fy(3x2-y2)
```

## Other files and directory
1.) The ```debug_rotateftn58_five.m``` demonstrate the rotation
2.) The ```illustration ver5.pptx``` illustrate the principles of this program
3.) The ```Wannier_orbits.pptx``` illustrate the principle of finding the representation of the rotation matrix
4.) The ```./bandplot_concise``` and ```./ftn58_related``` are programs supporting ```debug_rotatftn58_five.m``` 