# oxi_si.py
- This program makes SiO<sub>2</sub>/Si interface with lammps (using ```mpirun``` command).
- Before running this program, make sure the ```SiO.wyto``` file is in the same directory and check variables in this program (```lmp_path``` and ```mpi_num```)   .

## Usage
Oxidation of a $\times$ b $\times$ c Si supercell begins by typing the following command :

```
python oxi_si.py -n "a b c"
```

Oxidation starts with the closest layer to the surface and ends when half of the Si has been oxidized.

## Documentation
### Initial structure
- Make a vacuum region at the top of a $\times$ b $\times$ c Si supercell
- Make the initial structure by the following steps :

```
1. Structure optimization (cg)
2. NVT MD at 300K (Timestep 1fs, 10000step)
3. Structure optimization (cg)
```

### Oxidation
- Oxidation proceeds repeating the following steps (layer by layer oxidation) :

```
1. Insert O atom at the midpoint of Si-Si bonds (Randomly select 10% of the insertable positions)
2. Structure optimization (cg)
3. NVT MD at 1000K (Timestep 0.1fs, 10000step)
4. Structure optimization (cg)
```

### Result
The results of oxidation are grouped in directories as follows.

```
.
├── si_init/
│   ├── log.lammps
│   ├── stable.data
│   └── stable.final
└── oxidation/
    ├── layer1/
    │   ├── 1_md/
    │   │   ├── log.lammps
    │   │   ├── md.final
    │   │   ├── oxi.traj
    │   │   └── stable.data
    │   ├── 1_stable_1st/
    │   │   ├── log.lammps
    │   │   ├── stable.data
    │   │   └── stable.final
    │   ├── 1_stable_2nd/
    │   │   ├── log.lammps
    │   │   ├── stable.data
    │   │   └── stable.final
    │   ├── 2_md
    │   ├── 2_stable_1st
    │   ├── 2_stable_2nd
    │   ├── 3_md
    │   ├── 3_stable_1st
    │   ├── 3_stable_2nd
    │   └── ...
    ├── layer2
    ├── layer3
    └── ...
```
- File  
```stalbe.data``` : Input file  
```stable.final``` : Output file  
```md.traj``` : Trajectory during NVT MD  

- Directory  
```Si_init``` : Contain the results of making the initial structure  
```oxidation``` : Contain the results of oxidation for each layer  
```stable_1st``` : Contain the results of structure optimization after inserting O atom at the midpoint of Si-Si bonds  
```md``` : Contain the results of NVT MD, after 1st structure optimization  
```stable_2nd``` : Contain the results of structure optimization after NVT MD  

## Credit
- Kotaro Takematsu (Waseda University) : Code development
- Kentaro Hirai (Waseda University) : Code debugging assistance
- Takanobu Watanabe (Waseda University) : Code debugging assistance
