**A fork version of [ARCSim-HYLC](https://git.ist.ac.at/gsperl/ARCSim-HYLC) to support SMPL Model**

This repo is based on ARCSim-HYLC for the compatibility of SMPL, which  is also a fork version of arcsim

addons

* SMPL forward kinematics (based on [SMPL_cpp](https://github.com/soulslicer/smpl_cpp)) 
* SMPL integration with ARCSim (forwarded SMPL model as the obstacle)
* SMPL shape interpolation and pose interpolation to simulate any SMPL pose
* a `play` option like original `replay`, for checkout the SMPL motion and don't care the cloth

flags:
* `is_smpl`
* `motions`->`'smpl_fps`  how fast smpl model act 
* `motions`->`smpl_motfile`
* `motions`->`smpl_initialization_steps` (for shape and pose interpolation)

Thanks for @[isantesteban](https://github.com/isantesteban) for valuable discussions
