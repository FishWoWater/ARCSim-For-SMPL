# Homogenized Yarn-Level Cloth - ARCSim

This repository contains the macro-scale/cloth-simulation part of 'Homogenized Yarn-Level Cloth' (Siggraph 2020).

[Project Page: https://visualcomputing.ist.ac.at/publications/2020/HYLC/](https://visualcomputing.ist.ac.at/publications/2020/HYLC/)

For the micro-scale/pattern-optimization part, see the separate repository ["HYLC"](https://git.ist.ac.at/gsperl/HYLC/).

This is a fork / modified version of [ARCSim](http://graphics.berkeley.edu/resources/ARCSim/). The original readme file can be found [here](README_arcsim).

This repository includes source code within `v0.2.1/` (modified from the original ARCSim, especially the `hylc/` subfolder), the already fitted homogenized energies stored in `hylcmaterials/`, the simulation configuration files in `tests/`, and several additions to the `meshes/` folder.

## Usage

Check/call `exec.py` for how to compile and run a single simulation. The simulation config files are found in `tests/` subfolders.
E.g.: `python exec.py tests/2D/conf/stock_stretchX.json` (press space to start the simulation)

`python run_folder.py` can be used to (order and) simulate an entire folder of config files.

## Citation

Note that the [license](LICENSE) of ARCSim requires citing their work.
In addition, please consider citing our work if you find our modified version to be useful:
```bibtex
@article{sperl2020hylc,
  author    = {Sperl, Georg and Narain, Rahul and Wojtan, Chris},
  title     = {Homogenized Yarn-Level Cloth},
  journal   = {ACM Transactions on Graphics (TOG)},
  number    = {4},
  volume    = {39},
  year      = {2020},
  publisher = {ACM}
}
```
