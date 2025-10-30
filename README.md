# Target Report Repository
Reproducible target discovery repository. 

# Downloading and Using the Singularity Image from Dockerhub and setting up the file structure

Create the following folder structure on your linux system:

```
-demo
    ├── scripts
    ├── envs
    ├── workdir
    └── data
```
          
This singularity container is available on dockerhub: [Link](https://hub.docker.com/repository/docker/boeings/r450.python310.ubuntu.22.04)

The image can be imported as follows ( `--cleanenv` is optional in many settings):

```{bash}
## For R 4.5.0 (downloading the docker image to your folder - see below the option to access the singularity environment on the institute's system)

## Go to envs folder
cd envs

ml Singularity/3.6.4
singularity pull docker://boeings/r450.python310.ubuntu.22.04

cd ../workdir

# Starting the singularity container

You need to adjust the file paths for the --bind flag as is appropriate for your system. 
singularity shell --cleanenv --bind /nemo:/nemo,/camp:/camp ../envs/r450.python310.ubuntu.22.04_latest.sif

# Or create venv environment in the current folder (change path if you'd like to store the venv environment files elsewhere)
python3.10 -m venv ../envs/demo_venv_310

# Activate venv environment
source ../envs/demo_venv_310/bin/activate
## Save venv.lock documentation file
### Save venv package environment
pip freeze > ../envs/venv.lock
### Reload
pip install -r ../envs/venv.lock


# Start python
cd ../workdir
python

# Start R
cd ../workdir
R

```

Setting up Renv environment
```{R}

```
