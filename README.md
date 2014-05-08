NBODY6 Input File
------------------

### Structure of Input File

`1` `10000.0`  
`n` `1` `ncrit` `nrand` `300` `1`  
`0.02` `0.02` `0.15` `dtadj` `dtout` `tcrit` `1.0E-04` `2.0` `meanm`  
`2` `0` `O03` `0` `O05` `0` `2` `O08` `0` `0`  
`0` `O12` `0` `O14` `1` `1` `0` `0` `3` `O20`  
`2` `O22` `O23` `0` `1` `2` `3` `4` `0` `1`  
`0` `1` `2` `2` `1` `0` `0` `2` `0` `3`  
`0` `0` `0` `0` `0` `0` `0` `0` `0` `0`  
`2.0E-6` `1.0E-4` `0.2` `1.0` `1.0E-6` `0.01`  
`2.3` `mmax` `mmin` `nbin` `0` `Z` `epoch` `dtplot`  
`virq` `0.0` `0.0` `0.0`  
`0.125`

### Explanation of Parameters

`n` -- Total number of systems in cluster (singles + binary systems).  _Suggestion: 3000_

`ncrit` -- Number of stars where simulation stops. For clusters with numbers of stars below this value, it is no longer considered a cluster. _Suggestion: 10_

`nrand` -- Random seed used for simulation. You (should) get reproducable results if you specify the same seed for each simulation.

`dtadj` -- Energy checking time-step (Myr). _Suggestion: 2.0_

`dtout` -- Output interval (Myr). Must be a multiple of `dtadj`. _Suggestion: 10.0_

`tcrit` -- Length of simulation (Myr). _Suggestion: 10.0_

`meanm` -- Average stellar mass within cluster (Solar Masses). _Suggestion: 0.6_

`O03` -- Position / Velocity output flag.
> `0` -- Do not output any position or velocity information.  
> `1` -- Output rows of `m` `x` `y` `z` `vx` `vy` `vz` for each output timestep.  

`O05` -- Initial cluster shape. _Suggestion: 1_
> `0` -- Uniform and isotropic sphere.  
> `1` -- Plummer sphere.

`O08` -- Primordial binaries flag. _Suggestion: 0_
> `0` -- No primordial binaries.  
> `1` -- Primordial binaries. Must add extra line to input file.

`O12` -- H-R diagram output flag. _Suggestion: 1_
> `0` -- No H-R digram output.  
> `1` -- Stellar parameters necessary for making H-R digram are output to fort.82 and fort.83

`O14` -- Galactic external force model. _Suggestion: 1_
> `-1` -- Cut-off galactic force.  
> `1` -- Linearized galactic force.  
> `2` -- Point-mass galaxy force.  
> `3` -- Generalized galactic force (point mass + disk + logarithmic halo). Must add extra line to input file.

`O20` -- Initial mass function. _Suggestion: 4_

`O22` -- Manual initial conditions flag. _Suggestion: 0_
> `0` -- NBODY6-specified initial conditions.  
> `2` -- Manual input of intial conditions: `m` `x` `y` `z` `vx` `vy` `vz` on fort.10. Input scaled.  
> `3` -- Manual input of intial conditions: `m` `x` `y` `z` `vx` `vy` `vz` on fort.10. Input not scaled.

`O23` -- Escaper output flag. _Suggestion: 2_
> `0` -- No output for escaping stars.  
> `1` -- Basic output to screen for escaping stars.  
> `2` -- Diagnostic output to ESC for escaping stars, and escape angles printed to screen.








