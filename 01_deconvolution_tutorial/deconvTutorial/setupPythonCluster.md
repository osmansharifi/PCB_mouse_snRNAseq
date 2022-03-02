## Python setup on the cluster
In order to run deconvolution, you will need a python virtual environment with python3, tensorflow and tensorflow-probability. To setup this environment follow these steps:
```shell
## Create the virtual environment
virtualenv -p /usr/bin/python3 /share/<DIR>/<USERNAME>/pyvTf2

## Launch the virtual environment
source /share/<DIR>/<USERNAME>/pyvTf2/bin/activate ## Alias in your ~/.bash_profile to avoid retyping!

## Setup dependencies
pip3 install tensorflow
pip3 install tensorflow-probability
pip3 install sklearn

## Install <NAME>
pip3 install /share/quonlab/wkdir/njjohans/public/deconv_tutorial/deconv
```
Now anytime you want to run `<NAME` just launch the virtual environment on `pggb` before starting an R session!
