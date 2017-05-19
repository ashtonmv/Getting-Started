# Part 3: Set up environments on your computers
The objectives in this part assume that you have an account on Hipergator and/or Stampede, so if those aren't set up for you yet, wait to do these things until they are.

It's up to you how much of the following you want to install on your own computer. You might find it convenient to
run software on your own machine, or you might want to keep your computer clean and just run everything on the supercomputers. Below are instructions for setting up programming environments on your supercomputer accounts, but the instructions would be more or less the same if you want to get set up on your own computer.

--------
## Creating a virtual environment
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

##Installing Software
If the environment is activated, you're ready to start installing the Python software packages we use in our group. To make sure your miniconda environment is working, type ``which python``, and it should say ``~/miniconda3/bin/python``. If it doesn't, deactivate and reactivate miniconda until it does. Then install Pymatgen, MPInterfaces and iPython:

```
$ pip install pymatgen
$ pip install mpinterfaces
$ pip install ipython
```
**Pymatgen** is a large python package developed by a group called the Materials Project in California, and it has lots of very useful tools for computational materials research. **MPInterfaces** is a python package developed in our group that acts as an extension, or wrapper, around Pymatgen to make certain actions easier to perform. Finally, **iPython** is an interactive python shell that you will be using to execute commands from both packages.

## Create an account with the Materials Project
Now's as good a time as any to create an account with the Materials Project, which is a huge database of materials' crystal structures that you will probably find helpful in your research. Head over to [their website](https://materialsproject.org/) and create an account. After you've made an account, go to [your dashboard](https://materialsproject.org/dashboard) and generate an API key. You will need this later. Feel free to browse around their website to see what kind of information they have in their database. All of their calculations are the same kind as the ones we perform in our group.