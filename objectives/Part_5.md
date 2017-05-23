# Part 5: Running your first VASP calculation: Si relaxation
Now that you're basically familiar with the way Hipergator is set up, you can start submitting some of your own calculations.

-----------
## About VASP
The main software package our group uses to study materials is called [VASP](https://www.vasp.at/). It's not free, and it's not something you would install on your laptop. In fact, most users don't install it for themselves at all. Typically, one person in the group is in charge of installing various versions of VASP, and then other users just use the version he or she installed. In other words, the latest version of VASP is already installed for you, and you just need to know how to use it.

It's worth looking into more detail about what VASP actually does, but I'll try to describe it in simple terms here. VASP calculates the "ground state" of crystal structures. To do this, it needs to model every atom in a crystal structure, including all of its relevant electrons, and describe the way each atom interacts with each other atom. The electrons, which are quantum-mechanical particles, are tricky to model since most real materials have lots of them and we don't have an exact way to describe the way they interact with each other when there are a lot of them. So we make some approximations which we can actually calculate and which work pretty well at describing the electrons' behavior.

We often talk about using VASP to "relax" a structure. By this we mean that we're using VASP to calculate the atomic positions and lattice constants that give the structure its lowest possible energy. So even if we start out with a rough guess at what a crystal structure looks like, VASP can often relax it into the real crystal structure, provided the initial guess wasn't outrageously bad.

When you put all of VASP's capabilities together, it's a very reliable and useful way to study materials at the atomic level. VASP can give us insight into a material's mechanical properties, its electronic/ magnetic properties, and its thermodynamic properties. At different times, you will probably use VASP to study a material for each of those properties.

## Input files
VASP requires 4 input files to get started using these approximations to search for a material's ground state. These files are:
1. **POSCAR**. This file describes the crystal structure you're studying. An example POSCAR of silicon's crystal structure is given below:

```
Si2
1.0
2.734364 2.734364 0.000000
0.000000 2.734364 2.734364
2.734364 0.000000 2.734364
Si
2
direct
0.250000 0.250000 0.250000 Si
0.000000 0.000000 0.000000 Si
```
The first line in the POSCAR is a title, and can be anything you want. The second line is your lattice parameter, which is multiplied by all of your basis vectors on the next three lines. The fifth line has the names of all of the elements in your crystal structure in order. The sixth line tells the number of each kind of element. "direct" means that the coordinates of the atoms below are given with respect to the basis vectors. The lines after that are the coordinates of each atom.

2. **KPOINTS**. This file specifies the grid on which VASP will sample the electron density of your system. The main thing to know about this file is that the more k-points you include, the more accurate your calculation will be and the longer it will take. An example is given below:

```
KPOINTS
0
Gamma
10 10 10
```
Again, the first line is a title. The "0" on the second line means you want to generate the k-points automatically (this is almost always true). "Gamma" means to center the grid around the Gamma point, which is the origin ([0, 0, 0]) of your Brillouin zone. This is usually good to use, especially for hexagonal structures. The last line is the number of k-points you want along each dimension. So this file will make a grid of 10x10x10 k-points for sampling the electronic density of a structure.

3. **INCAR**. This file contains the instructions telling VASP what kind of calculation to run. There are many options that can be set here and I can't possibly go through all of them but they are online in the [VASP wiki site](https://cms.mpi.univie.ac.at/wiki/index.php/The_VASP_Manual) under INCAR-tags. Below is a very simple example INCAR for relaxing the Si structure above:

```
EDIFF = 1e-04
ENCUT = 350
IBRION = 2
ISIF = 3
ISMEAR = 0
NSW = 50
PREC = Accurate
SIGMA = 0.1
```

- **EDIFF**: The criteria for determining when your calculation has converged
- **ENCUT**: The cutoff energy for the plane-wave basis set. Higher = more accurate and takes longer
- **IBRION**: The algorithm used to relax your structure
- **ISIF**: How many degrees of freedom to use while relaxing your structure
- **ISMEAR**: How VASP should handle partially occupied electronic states
- **NSW**: The number of ionic steps your relaxation is allowed to run
- **PREC**: How precise your calculation will be
- **SIGMA**: Width in eV for smearing partial occupancies

It might also be helpful if you look up what each tag in that INCAR means in the VASP wiki page.

4. **POTCAR**. This file contains information about the core and valence electrons of each element in your crystal structure. You don't actually write POTCAR files- these files come with VASP and we're not allowed to distribute them so contact someone in the group about where to find them on Hipergator.

## Running a calculation
Make sure you're on your scratch filesystem (``cd /ufrc/hennig/your_username``) and then make a directory called ``Si_test``. Go into that directory and use vim or emacs to write the contents of the **POSCAR**, **KPOINTS** and **INCAR** above to those files in your directory. Then copy the Si POTCAR from the location you were told about on Hipergator to your directory. Your directory should now have the following contents:

```
POSCAR INCAR KPOINTS POTCAR
```
Then you can submit your job to the compute nodes with the submission script below:

```bash
#!/bin/bash
#SBATCH --job-name=Si_relax
#SBATCH -o out_%j.log
#SBATCH -e err_%j.log
#SBATCH --qos=hennig
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=1000mb
#SBATCH -t 10:00:00

cd $SLURM_SUBMIT_DIR

module load intel/2016.0.109
module load openmpi

mpirun /home/mashton/vasp.5.4.1/bin/vasp > job.log
```
So copy that to a file named "runjob", and then type ``sbatch runjob`` from your command line. Hipergator should respond with "Submitted batch job XXXXXXX" if everything worked okay.

Check on the status of your job periodically. Once it starts running, it won't take long to finish. The section below has some tips about understanding the output files VASP creates, so go through those after the job has finished.

## Output files
After your job has finished, you'll notice that there are about 15-20 new files in your ``Si_test`` directory. Don't worry about all of them, we'll just go over the most important ones.

**--- job.log ---**

This is the file we asked the compute node to write our job's main output to (look at the last line in your runjob file), so it contains a log of how the calculation ran. We should check this file to make sure the calculation converged properly.

There are a lot of lines and a lot of numbers in this file. The top of the file will look something like this:

```
running on   16 total cores
distrk:  each k-point on   16 cores,    1 groups
distr:  one band on    1 cores,   16 groups
using from now: INCAR
vasp.5.4.1 24Jun15 (build May 10 2016 14:40:36) complex

POSCAR found type information on POSCAR  Si
POSCAR found :  1 types and       2 ions
scaLAPACK will be used
```
giving you some basic information about the job. Below that might be some warnings- you can look at those for now but don't let them alarm you. Then you will see several blocks like this:

```
N       E                     dE             d eps       ncg     rms          rms(c)
Total vdW correction in eV:    5.44740944070352
DAV:   1     0.142822335994E+03    0.14282E+03   -0.13325E+04  7808   0.672E+02
DAV:   2     0.620650629458E+01   -0.13662E+03   -0.12063E+03  7568   0.152E+02
DAV:   3    -0.839822378505E+01   -0.14605E+02   -0.14484E+02  8464   0.499E+01
DAV:   4    -0.860551850367E+01   -0.20729E+00   -0.20725E+00  7312   0.733E+00
DAV:   5    -0.861013180243E+01   -0.46133E-02   -0.46132E-02  8848   0.120E+00    0.123E+01
Total vdW correction in eV:    5.46516819513896
RMM:   6    -0.726373637170E+01    0.13464E+01   -0.10124E+01  6491   0.153E+01    0.136E+01
Total vdW correction in eV:    5.47503716575538
RMM:   7    -0.690783802015E+01    0.35590E+00   -0.46268E-01  6983   0.384E+00    0.791E+00
Total vdW correction in eV:    5.49248204046353
RMM:   8    -0.686508522000E+01    0.42753E-01   -0.17091E-01  6384   0.240E+00    0.402E+00
```
Each row in this part of the output is VASP calculating the total energy of the system based on its most recent guess at the electronic wavefunctions. In other words, VASP is trying out different electronic configurations for your material, and if a new configuration gives it a lower energy, VASP takes that as the new configuration and keeps going so on and so on until the changes in the system's energy are lower than the EDIFF you specified in your INCAR file. DAV and RMM are abbreviations for the names of the algorithms being used to update the wavefunctions between steps. The column labeled "E" are the energies VASP calculated, and the column labeled "dE" are the differences between that electron configuration's energy and the one prior. Don't worry about the other columns for now.

If you keep scrolling down, you will see that these lines end when dE < EDIFF, and a line like this one is printed:

```
   1 F= -.68422849E+01 E0= -.68425880E+01  d E =-.684228E+01
```
VASP has found an electronic configuration that is sufficiently low in energy that it stopped, and now it needs to calculate the forces acting on each atom so it can properly move them around and do the whole thing over again. This is called an ionic step. After VASP has moved each atom according to the forces acting on it, it recalculates the wavefunctions based on those new atomic positions. This "self-consistent" cycle continues until the TOTAL energy changes between ionic steps are less than EDIFF * 10 (by default, but you can change this). After that criteria is met, you should see a line that says

```
 reached required accuracy - stopping structural energy minimisation
```
That means everything worked, your relaxation is done, and you have an optimized crystal structure. If you see any errors in your job.log file, you can ask someone else in the group about them- we have almost definitely run into the same error before. We're very good at running into errors.

**--- OUTCAR ---**

The OUTCAR file has pretty much all of the information about your calculation's results in a human-readable format. Almost everything in the OUTCAR file is also contained in other, more specific output files, but sometimes it's nice to have everything in one place.

**--- CHG, CHGCAR, and WAVECAR ---**

These files are not necessarily human readable, and are used to store the charge density and wavefunctions of your system in case you want to re-run the calculation in the future. These are typically the 3 largest files generated by VASP, so if you don't think you'll use them again you should probably delete them.

**--- CONTCAR ---**

This is exactly like the POSCAR, but it contains the optimized structure instead of your initial guess. Compare it with the POSCAR to see how much your structure changed during the relaxation.

**--- DOSCAR ---**

The DOSCAR contains the density of states for your material's electronic structure. For the calculation you just ran, the results in this file won't be extremely meaningful.

**--- OSZICAR ---**

This file is just like the job.log, but with a little less clutter (and, consequently, information).

**--- vasprun.xml ---**

The vasprun.xml is not designed for human-readability, but is very easy for a Python program like Pymatgen to parse. Like the OUTCAR, it contains all of the information about your calculation.
