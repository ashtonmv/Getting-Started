# Part 7: Using MPInterfaces to study 2D materials

Silicon was kind of a boring example, right? Chances are your research project will involve studying some kind of 2D materials, since they are a major focus in the Hennig group. Below you will get a chance to perform VASP calculations for a 2D material (of your choosing! See how not boring we're being now?) similar to what you did for silicon a moment ago.

-------------
## Choosing a 2D material
It's possible you already have a POSCAR file on hand for a 2D material, but if you don't you can download one from [MaterialsWeb](https://materialsweb.org/twod_materials), our group's website. If you're unsure which materials are in our database, enter "{}" into the search bar on MaterialsWeb and wait a second. It will give you a list of the hundreds of 2D materials that are in our database. For now, just pick one out that has 6 or fewer atoms, and 2 or fewer elements in its composition. Click on it and download its POSCAR file from its data page.

Upload that POSCAR file to your scratch filesystem on Hipergator:

~~~bash
$ scp Downloads/POSCAR your_username@hpg2.rc.ufl.edu:/ufrc/hennig/your_username/
~~~

## Relaxing the 2D material and getting its band structure
Then log onto Hipergator and navigate to your scratch directory. There, make a new directory for relaxing your 2D material. You might as well name it after the 2D material's formula. I'll use VSe2 (mw-64) as an example.

```shell
$ mkdir VSe2
$ mv POSCAR VSe2
$ cd VSe2
$ ipython
```
```python
In [1]: from mpinterfaces.mat2d.stability.startup import relax
In [2]: relax()
In [3]: exit
```

I'll mention here that the VASP execution statement in the runjob file is different for this material than it was for Si. Instead of ``..../bin/vasp``, it should be ``..../bin/vasp_noz`` or something similar. That's important because 2D materials should have vacuum above and below them in their POSCAR file. To visualize this, load the POSCAR file you downloaded in to VESTA or whatever crystal structure visualization software you have. Normally when VASP is relaxing a structure, it would collapse that vacuum. ``vasp_noz`` is a special version of VASP designed specifically not to collapse vacuum regions.

After your material has finished relaxing, check its job.log to make sure it converged okay. Technically, if you downloaded a structure from MaterialsWeb, it had already been relaxed by MPInterfaces before it was uploaded, so your relaxation should not take too long. If everything converged okay, then go ahead and run a band structure calculation just like you did for silicon before.

```shell
$ ipython
```
```python
In [1]: from mpinterfaces.mat2d.electronic_structure.startup import run_pbe_calculation
In [2]: run_pbe_calculation()
In [3]: exit
```

When it's done, plot the band structure:

```shell
$ ipython
```
```python
In [1]: from mpinterfaces.mat2d.electronic_structure.analysis import plot_band_structure
In [2]: plot_band_structure()
In [3]: exit
```

Is your material a metal or a semiconductor? Its band structure plot should look like the one on its data page on MaterialsWeb.

## Understanding thermodynamic stability
When we study 2D materials, it's important to keep in mind that they are actually not stable states of matter. They're never as stable as normal bulk materials, but if their energy is close enough to the bulk materials, then they can exist as metastable. This difference in energy between the 2D material and bulk materials is usually referred to as the 2D material's *hull distance*, for reasons I'll explain below. Afterward, you'll get to calculate the hull distance for your 2D material.

Open the [phase diagram app](https://materialsproject.org/#apps/phasediagram) on Materials Project. Enter the elements in your material. For me, I entered V and Se. After pressing "Generate", you should see a "v"-shaped plot. This plot is called the "convex hull" (for the V-Se system, in my case). The materials that lie on the convex hull are thermodynamically perfectly stable. To the right of the plot, click the "Unstable" tab. You'll probably see several light blue points show up above the convex hull. These materials are unstable, because their energies are higher than those of the stable phases on the convex hull. In real life they probably won't exist since their energy could be lowered by simply decomposing to one or more of the phases on the convex hull. The *hull distance* I referred to earlier is the length of the straight line downward from an unstable material to the convex hull. If this line is short enough, *i.e.* the material's energy is close enough to the real stable phases, it might exist as metastable. In the past we've determined that ~0.15-0.2 eV/atom hull distance is a reasonable cutoff for determining whether or not a 2D material might exist as metastable.

## Calculating the material's hull distance
To calculate our 2D material's hull distance, we need to download and relax the structures of the bulk materials below it on the convex hull diagram. MPInterfaces can do this automatically for us. Go into your 2D material's relaxation directory and run the following:

```shell
$ ipython
```
```python
In [1]: from mpinterfaces.mat2d.stability.analysis import get_competing_phases
In [2]: print(get_competing_phases())
[('VSe2', 'mp-694')]
In [3]: exit
```
That is the formula and Materials Project ID of the material directly below your 2D material on the convex hull. We need its energy, so let's download it and relax it. First make a directory called ``competing_phases`` in your top-level scratch directory (/ufrc/hennig/your_username/). Then go into that directory and make a directory named after the formula of the material (*e.g.* VSe2). To download its structure from the Materials Project, go into that new directory and do the following:

```shell
$ ipython
```
```python
In [1]: from pymatgen.matproj.rest import MPRester
In [2]: from mpinterfaces import CONFIG_FILE
In [3]: MPR = MPRester(CONFIG_FILE["mp_api"])
In [4]: MPR.get_structure_by_material_id("mp-694").to("POSCAR", "POSCAR")
In [5]: from mpinterfaces.mat2d.stability.startup import relax
In [6]: relax(dim=3)
In [7]: exit
```
Now your competing bulk phase is relaxing. Once it's done, you can compare its energy with the energy of your 2D material to get the hull distance by doing the following from within your 2D material's relaxation directory:

```shell
$ ipython
```
```python
In [1]: from mpinterfaces.mat2d.stability.analysis import get_hull_distance
In [2]: print(get_hull_distance())
0.09
In [3]: exit
```
That value is in eV/atom, so you can see that VSe2 is well within the 0.15 eV/atom cutoff. Compare your calculated hull distance with the Formation Energy listed on its data page on MaterialsWeb (formation energy is, in this case, just another name for hull distance).