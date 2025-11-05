# Target Report Repository
Reproducible target discovery repository. 

## Downloading and Using the Singularity Image from Dockerhub and setting up the file structure

Create the following folder structure on your linux system:

```
-demo
    ├── scripts
    ├── envs
    ├── workdir
    └── data
```
          
The singularity container used in this work is available on dockerhub: [Link](https://hub.docker.com/repository/docker/boeings/r450.python310.ubuntu.22.04)

Navigate to the `envs` directory. 

The commands below contain system specific components which you will have to adapt to your system. 

```{bash}
## Go to envs folder ##
cd envs

## Activate the singularity module ##
ml Singularity/3.6.4

## Pull singularity container from dockerhub
singularity pull docker://boeings/r450.python310.ubuntu.22.04

## Navigating to the workdir ##
cd ../workdir

## Starting the singularity container ##
You need to adjust the file paths for the --bind flag as appropriate on your system.

General syntax:
singularity shell --cleanenv --bind </folder/on/host/system/data>:/data ../envs/r450.python310.ubuntu.22.04_latest.sif

Specific example:
singularity shell --cleanenv --bind /nemo:/nemo,/camp:/camp ../envs/r450.python310.ubuntu.22.04_latest.sif

## Or create venv environment in the current folder (change path if you'd like to store the venv environment files elsewhere) ##
This is only required once at the beginning of the project, so will be commented out for now. 
# python3.10 -m venv ../envs/demo_venv_310

## Activate venv environment ##
source ../envs/demo_venv_310/bin/activate

## Save venv.lock documentation file ##
## Save venv package environment     ##
pip freeze > ../envs/venv.lock
## Reload venv environment from venv.lock file ##
pip install -r ../envs/venv.lock

## Start python ##
cd ../workdir
python

OR
## Start R ##

## optional: set R package paths ##
cd ../workdir
R

```

