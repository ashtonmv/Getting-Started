# Part 3: Set up environments on your computers
The objectives in this part assume that you have an account on Hipergator and/or Stampede, so if those aren't set up for you yet, wait to do these things until they are.

It's up to you how much of the following you want to install on your own computer. You might find it convenient to
run software on your own machine, or you might want to keep your computer clean and just run everything on the supercomputers. Below are instructions for setting up programming environments on your supercomputer accounts, but the instructions would be more or less the same if you want to get set up on your own computer.

--------
## Create a virtual environment
The cleanest way to get a Python environment up and running is to use a Miniconda virtual environment. You can install Miniconda on your account by downloading the Linux 64-bit bash installer [here](https://conda.io/miniconda.html), uploading it to Hipergator, and then running it. To upload files to/from Hipergator, open your terminal (mac/linux users) or Bash on Ubuntu (windows users) and enter the command:

```bash
$ scp path/to/file/on/your/computer your_username@hpg2.rc.ufl.edu:/path/where/you/want/the/file/on/hipergator
```
for example,

```bash
$ scp Downloads/Miniconda3-latest-Linux-x86_64.sh mashton@hpg2.rc.ufl.edu:/home/mashton/
```
will put the file in my home directory. It will ask for your Hipergator password, so type that in.

To log into Hipergator, use this command, followed by your password:

```bash
$ ssh your_username@hpg2.rc.ufl.edu
```

That will put you in your home directory. For me, that's where my file is. You can see if the file is there by typing ``ls``. To execute the script you just uploaded, you will need to change its permissions. To do this, type:

```shell
$ chmod 755 Miniconda3-latest-Linux-x86_64.sh
$ ./Miniconda3-latest-Linux-x86_64.sh
```
That will execute the script, so just follow the instructions that it gives you after that. The default answers to every question should work just fine.

Now you should have a new directory, ``miniconda3``. To activate the miniconda environment, type

```
$ source miniconda3/bin/activate
```
and to deactivate it, simply type ``source deactivate``.

## Install Pymatgen and ipython with pip
If the environment is activated, you're ready to start installing the Python software packages we use in our group. To make sure your miniconda environment is working, type ``which python``, and it should say ``~/miniconda3/bin/python``. If it doesn't, deactivate and reactivate miniconda until it does. Then install Pymatgen, MPInterfaces and iPython:

```
$ pip install pymatgen
$ pip install ipython
```
`pip` is a convenient tool that came with your miniconda virtual environment and installs Python packages that are listed on the Python package index. **Pymatgen** is a large python package developed by a group called the Materials Project in California, and it has lots of very useful tools for computational materials research. **iPython** is an interactive python shell that you will be using to execute commands from both packages.


## Create an account with the Materials Project
Now's as good a time as any to create an account with the Materials Project, which is a huge database of materials' crystal structures that you will probably find helpful in your research. Head over to [their website](https://materialsproject.org/) and create an account. After you've made an account, go to [your dashboard](https://materialsproject.org/dashboard) and generate an API key. You will need this key in just a second. Feel free to browse around their website to see what kind of information they have in their database. All of their calculations are the same kind as the ones we perform in our group.

## Install MPInterfaces
**MPInterfaces** is a software package developed in our group that acts as a wrapper around Pymatgen to make certain tasks even easier. Visit the [Github page](https://github.com/henniggroup/mpinterfaces) to check it out. You're about to install that software on your Hipergator account. Navigate to your `software` directory under your home directory and issue the following commands:

```
$ module load git
$ git clone https://github.com/henniggroup/MPInterfaces.git
```
You'll notice that you have a new directory called `MPInterfaces`, and its contents are exactly the same as what's on the Github page. To let Python know that you've installed MPInterfaces, you need to add its location to your $PYTHONPATH environment variable ([learn about environment variables](https://www.digitalocean.com/community/tutorials/how-to-read-and-set-environmental-and-shell-variables-on-a-linux-vps)). To do this, open your `~/.bashrc` file with vim or emacs and add the following line:

```
export PYTHONPATH=$PYTHONPATH:/home/your_username/software/MPInterfaces
```
As a quick tip, your `~/.bashrc` file lists any commands you want to run every time you log into Hipergator. So you could add a command to source your miniconda environment to this file if you wanted to, or have it `echo` you a welcome message every time you log in. To activate the changes you just made, type `source ~/.bashrc`. If you type `echo $PYTHONPATH`, it should end with `/home/your_username/software/MPInterfaces`.

MPInterfaces is a tool that will set up and launch VASP jobs for you based on some reasonable default INCAR & KPOINTS parameters. The only thing it can't do is predict which structure you want to study, since that would require telepathy. Before it can launch jobs for you, though, it needs to know a few things about you. This is not a one-way relationship, after all.

You need to edit the mpint_config.yaml file that came with MPInterfaces to contain your settings. First, copy the config file to its permanent location:

```
$ cp /home/your_username/software/MPInterfaces/config_files/mpint_config.yaml /home/your_username/software/MPInterfaces/mpinterfaces/
```
Then we'll edit the version of the configuration file under `/home/your_username/software/MPInterfaces/mpinterfaces/` with the following (use vim or emacs to edit the file):

```
username: # your UF HPC username
mp_api: # your Materials Project API key. You can find your API key at https://materialsproject.org/open
normal_binary: /home/mashton/vasp.5.4.1/bin/vasp # /path/to/std/vasp/
executable. You are probably okay putting /home/mashton/vasp.5.4.1/bin/vasp
twod_binary: /home/mashton/vasp.5.4.1/bin/vasp_noz  # /path/to/2D/vasp/executable
vdw_kernel: /home/mashton/vasp.5.4.1/vdw_kernel.bindat  # /path/to/vdw_kernel.bindat file. Leave as null if the kernel is hard-coded into VASP. Ask the VASP manager in the group about this.
potentials: /home/mashton/POTCARS  # /path/to/your/POTCAR/files
queue_system: slurm
queue_template: config_files/
```

Save and quit that file, and now you're ready to use MPInterfaces.