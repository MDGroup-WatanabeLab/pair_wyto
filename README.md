# pair_wyto
- Repository to develop code for Lammps pair style of extended SW potential for Si-O mixed system
  - The first version of the potential function was published in [1] and later revised in [2].

## Usage
### Build LAMMPS
1. Install the following two files in the `src` directory of LAMMPS.
   - `pair_wyto.cpp`
   - `pair_wyto.h`
2. Build LAMMPS.

### Simulation in LAMMPS
- In addition to the input script and structure file, prepare a `.wyto` file containing parameter information.
  - There is `SiO.wyto` in the `/test` directory which describes the parameters proposed in the paper.

## Test
- A test script for lammps is provided in the `/test` directory.

|File name|Description|
---|---
|`in.wyto_test`|Script for performing NVT MD of Si/SiO2 interface with 3 different time steps|
|`SiO.wyto`|Potential parameter file for SiO system|
|`SiO2_Si_interface`|Structure file of SiO2/Si interface|
|`check_econserve.py`|Python script to plot the Hamiltonian transitions in the NVT MD at each time step|

- In this test, please check that the Hamiltonian is conserved more accurately with shoter time steps.

## Credit
- Kentaro Hirai (Waseda University): Code development
- Takanobu Watanabe (Waseda University): Code debugging assistance

## Reference
[1] Watanabe, T., Fujiwara, H., Noguchi, H., Hoshino, T., & Ohdomari, I. (1999). Jpn. J. Appl. Phys., 38(4A), L366.   
[2] Watanabe, T., Yamasaki, D., Tatsumura, K., & Ohdomari, I. (2004). Appl. Surf. Sci, 234(1-4), 207-213.
