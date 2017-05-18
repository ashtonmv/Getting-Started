# Part 5: Running your first VASP calculation
Now that you're basically familiar with the way Hipergator is set up, you can start submitting some of your own calculations.

-----------
## About VASP
The main software package our group uses to study materials is called [VASP](https://www.vasp.at/). It's not free, and it's not something you would install on your laptop. In fact, most users don't install it for themselves at all. Typically, one person in the group is in charge of installing various versions of VASP, and then other users just use the version he or she installed. In other words, the latest version of VASP is already installed for you, and you just need to know how to use it.

It's worth looking into more detail about what VASP actually does, but I'll try to describe it in simple terms here. VASP calculates the "ground state" of crystal structures. To do this, it needs to model every atom in a crystal structure, including all of its relevant electrons, and describe the way each atom interacts with each other atom. The electrons, which are quantum-mechanical particles, are tricky to model since most real materials have lots of them and we don't have an exact way to describe the way they interact with each other when there are a lot of them. So we make some approximations which we can actually calculate and which work pretty well at describing the electrons' behavior.

We often talk about using VASP to "relax" a structure. By this we mean that we're using VASP to calculate the atomic positions and lattice constants that give the structure its lowest possible energy. So even if we start out with a rough guess at what a crystal structure looks like, VASP can often relax it into the real crystal structure, provided the initial guess wasn't outrageously bad.

## Input files
VASP requires 4 input files to get started using these approximations to search for a material's ground state. These files are:
1. **POSCAR**. This file describes the crystal structure you're studying. An example POSCAR of silicon's crystal structure is given below:
```
Si2
1.0
-0.000000 2.734364 2.734364
2.734364 0.000000 2.734364
2.734364 2.734364 0.000000
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
NSW = 1
PREC = Accurate
SIGMA = 0.1
```
It's actually probably more helpful if you just look up what each tag in that INCAR means in the VASP wiki page.

4. **POTCAR**. This file contains information about the core and valence electrons of each element in your crystal structure. You don't actually write POTCAR files- these files come with VASP and we're not allowed to distribute them so contact someone in the group about where to find them on Hipergator.

## Running a calculation
Make sure you're on your scratch filesystem (``cd /ufrc/hennig/your_username``) and then make a directory called ``Si_test``. Go into that directory and use vim or emacs to write the contents of the **POSCAR**, **KPOINTS** and **INCAR** above to those files in your directory. Then copy the Si POTCAR from the location you were told about on Hipergator to your directory. Your directory should now have the following contents:
```
POSCAR INCAR KPOINTS POTCAR
```
Then you can submit your job to the compute nodes with the submission script below:
```
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

echo 'Done.'
```
So copy that to a file named "runjob", and then type ``sbatch runjob`` from your command line. Hipergator should respond with "Submitted batch job XXXXXXX" if everything worked okay.

Check on the status of your job periodically. Once it starts running, it won't take long to finish. The section below has some tips about understanding the output files VASP creates, so go through those after the job has finished.

## Output files
