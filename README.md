# csMAHN

# Create conda environment

```shell
hostnamectl
#   ...
#   Operating System: CentOS Linux 7 (Core)
#        CPE OS Name: cpe:/o:centos:centos:7
#             Kernel: Linux 3.10.0-957.el7.x86_64
#       Architecture: x86-64
#   ...

conda --version
# conda 24.1.2

conda config --show channels
# channels:
#   - conda-forge
#   - bioconda
#   - defaults

# create the environment
conda activate base
# conda remove --name publish --all --yes
conda create --name publish --yes
conda activate publish
conda install --name publish --file require_publish --yes
pip install torch==2.0.1 torchvision==0.15.2 torchaudio==2.0.2 --index-url https://download.pytorch.org/whl/cpu
pip install came samap holoviews-samap
# create jupyter kernel
conda activate publish
python -m ipykernel install --user --name 'publish' --display-name 'publish'
##  referent https://irkernel.github.io/installation/#linux-panel
echo 'IRkernel::installspec(user = TRUE,name="R_publish",display="R_publish")' > temp_r_script.r
Rscript temp_r_script.r
rm temp_r_script.r
## check jupyter kernel
conda activate
jupyter kernelspec list
# Available kernels:
#   ...
#   publish           ....local/share/jupyter/kernels/publish
#   r_publish         ....local/share/jupyter/kernels/r_publish
#   ...

```
