echo "copy conda env for FELINE project: the folder is ~108G"
sbatch FELINE_conda_env.sh

echo "a shell script records how I install each env"
envs.install.sh
