# General Recipes in Python and Linux

This is a (non-comprehensive) list of recipes for performing tasks in Python and Linux that users in the group have found useful over the years.

## Linux recipes
**--- Create a shortcut command ---**

There are a few ways to do this, but `alias` is an easy one:

```
$ alias shortcut_name='actual command syntax'
```
e.g.

```
$ alias myq='squeue -u your_username'
$ alias devnode='module load ufrc; srundev --time=2:00:00'
```

**--- Set an evironment variable ---**

`export` behaves similarly to `alias` Use all capitals when defining variables in bash.

```
$ export VARIABLE_NAME=some_value
```
e.g.

```
export PYTHONPATH=/home/your_username/software/MPInterfaces:/home/your_username/software/Pymatgen
export HDGN='your_username@hydrogen_ip_address'
```

**--- Check the size of a directory ---**

```
$ du -h  # Recursively lists the size of every directory in the current directory
$ du -h --max-depth=N  # Change "N" to be however deep you want the recursion to run
```

**--- Check the status of the group's resources on Hipergator ---**

```
$ module load ufrc
$ slurmInfo -g hennig
```


## Python Recipes
**--- Make a list of all directories in the current one ---**

```python
import os
directories = [d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]
```

**--- Perform Linux-like commands ---**

```python
import os
os.mkdir(dir_name)  # Create a directory
os.chdir(dir_name)  # Go into a directory
os.system('literal command to execute')  # e.g. os.system('rm WAVECAR')
```

**--- Use [matplotlib](https://matplotlib.org/) to make plots ---**

```python
import matplotlib
matplotlib.use("Agg")  # This is important
import matplotlib.pyplot as plt


ax = plt.figure().gca()
ax.plot([1, 2, 3, 4, 5], [1, 4, 9, 16, 25])
plt.savefig('quadratic.pdf', transparent=True)
```