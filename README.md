**A fork version of [ARCSim-HYLC](https://git.ist.ac.at/gsperl/ARCSim-HYLC) to support SMPL Model**

This repo is based on ARCSim-HYLC for the compatibility of SMPL, which  is also a fork version of arcsim

addons

* SMPL forward kinematics (based on [SMPL_cpp](https://github.com/soulslicer/smpl_cpp)) 
* SMPL integration with ARCSim (forwarded SMPL model as the obstacle)
* SMPL shape interpolation and pose interpolation to simulate any SMPL pose (for the details of interpolation, refer to [Discussion](https://github.com/isantesteban/vto-dataset/issues/1))
* a `play` option like original `replay`, for checkout the SMPL motion and don't care the cloth

fields:
* `is_smpl`
* `motions`->`smpl_fps`  how fast smpl model act 
* `motions`->`smpl_motfile`
* `motions`->`smpl_initialization_steps` (for shape and pose interpolation)



Usage

* Simulate a simple motion `hands raise up` with shape / pose interpolation 

  * `./build-Release/bin/arcsim_0.2.1 simulate smpl_demo/jump.json`

* Just checkout the motion and save each screen shot to some directory(optional)

  * `./build-Release/bin/arcsim_0.2.1 play smpl_demo/jump.json vis`

  * the visualization will be saved to `vis` directory, this is a useful tool for the visualization of MoCap sequences

    

Demo(simulate mode)

![](imgs/demo_jump.mp4)



TODOs

- [x] Support SMPL female model (now only male model is supported)
- [ ] Support CMU MoCap sequence reading from npz or json file (now read from txt motion file)



Thanks for @[isantesteban](https://github.com/isantesteban) for valuable discussions
