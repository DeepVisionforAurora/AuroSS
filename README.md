# Automatically Sketching AuroSS

This is a source code used in the paper "Automatically Sketching Auroral Skeleton Structure in All-sky Image for Measuring Aurora Arcs"<br>

Implementation for AuroSS detection needs:<br>
* psudoLGeneration_ASI2Ridge_wq.m (by calling tfridge.m or ridge_AreaOnASI.m). The scripts are used to generate skeleton pseudo-labels for training ASMs. Based on the probability map of the auroral regions output by the AAM, the ridges of each luminous region are detected.<br>
* U-net for skeleton extraction and aurora detection (pixel identification) is based on image segmentation [1].      U-net architecture is available at https://github.com/orobix/retina-unet.<br>
* Orientation_ASI2MLON_MLAT_wq.m This code is used to calculate the arc tilt to measure auroral orientation. The obtained AuroSS is converted into geomagnetic coordinate reference frame. We calculate the arc tilt to measure auroral orientation.<br>
* Cycle‐consistent Generative Adversarial Network (CycleGAN) [2] is used to generate ASI images. The method is used to test the effectiveness of our method on synthetic data. The code is available at https://github.com/junyanz/CycleGAN.<br>

## How To Use

Download the above code and packages and install the required toolbox and dependencies (detailed in Readme in each package).<br>
A demo to generate skeleton pseudo-labels is provided. Two optional maximum energy based ridge extraction methods are provided. Manually selected better ones are used as pseudo-labels to train the ASMs.<br>

For training U-net, the training data is organized as:<br>
ASI image dataset<br>
│<br>
└───test<br>
|    ├───Annotation<br>
|    └───images<br>
│<br>
└───training<br>
     ├───Annotation<br>
     └───images<br>
     
To train the AAM, the "Annotation" folder is a manually labeled image of the auroral area; to train ASM, the "Annotation" folder is the generated AuroSS pseudo-labels. Then run:<br>
   ``` python
   python prepare_datasets.py
```<br>

Specify the parameters in configuration.txt. Training and inference run by:<br>
   ``` python
   python run_training.py
   python run_testing.py
```
After obtaining the AuroSS which is saved under folder "skeImgs", run Orientation_ASI2MLON_MLAT_wq.m to calculate the arc tilt to measure auroral orientation.<br>
We used manually labeled skeleton images and real ASI images as training sets to train CycleGAN, run:<br>
      ```
      DATA_ROOT=./ASIdatast name=ASI_model th train.lua
      ```

## Data

The aurora images in this dataset are from the all-sky imagers at the Yellow River Station (YRS). And it is used in the paper Automatically sketching Auroral Skeleton Structure in All-sky Image for Measuring Aurora Arcs. The data source, acquisition instruments, data characteristics, preprocessing methods, and naming conventions of this dataset are detailed at https://doi.org/10.5281/zenodo.10464846.<br>

## Reference

[1] Ronneberger, Olaf, Philipp Fischer, and Thomas Brox. "U-net: Convolutional networks for biomedical image segmentation." Medical Image Computing and Computer-Assisted Intervention–MICCAI 2015: 18th International Conference, Munich, Germany, October 5-9, 2015, Proceedings, Part III 18. Springer International Publishing, 2015.<br>

[2]  Zhu, Jun-Yan, et al. "Unpaired image-to-image translation using cycle-consistent adversarial networks." Proceedings of the IEEE international conference on computer vision. 2017.
