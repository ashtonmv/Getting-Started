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

**--- Cancel a large number of SLURM jobs at once ---**

Sometimes you goofed and submitted dozens or hundreds of jobs that you don't actually want to run. Instead of canceling them one at a time, you can cancel all of them at once. Copy the contents below to a file called `cancel_jobs` (or whatever you want to name it) located in `/home/your_username/bin/`.

```
#!/bin/bash

squeue -u your_username > jobs_qstat
grep 'your_username' jobs_qstat | cut -c 12-18,48-49 > job_ids

while read line; do
    ID=$(echo "$line" | cut -c 1-7)
    STAT=$(echo "$line" | cut -c 8-9)

    if [ "$ID" -ge "$1" ] && [ "$ID" -le "$2" ] && [ "$STAT" != "C" ]; then
        scancel $ID
    fi

done <job_ids

rm job_ids
rm jobs_qstat
```
If the following is not already in your ~/.bashrc, then add it:

```
if [ -d "$HOME/bin" ] ; then
    PATH="$HOME/bin:$PATH"
fi
```
Finally, `chmod 755 /home/your_username/bin/cancel_jobs`. To use this tool, type `cancel_jobs JOB_ID_1 JOB_ID_2`, where `JOB_ID_1` and `JOB_ID_2` should be the JOB_ID's of the first and last jobs defining the range you want to cancel.

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

**--- Save an iPython session to a file ---**

Have you ever had a really brilliant iPython session that went dozens of lines long and then later wished you had saved it to a `.py` file so you could run it again? You're not alone. Here's how to save an iPython session to a file:

```python
...
...
...
In [50]: %save filename 1-49
The following commands were written to file `filename`:
...