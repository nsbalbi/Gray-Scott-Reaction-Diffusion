# Gray-Scott-Reaction-Diffusion

> MATLAB function used to generate and record Gray-Scott reaction-diffusion models

![Sample GIF 1](https://github.com/nsbalbi/Gray-Scott-Reaction-Diffusion/blob/master/Sample%20GIF%201.gif)

## Usage

Simply call the function with the desired system constants and the model will be generated and recorded! This function also allows for the customization of model grid size (performance limited), initial state, colormap, video name and file format, and many other paramters. Robert Munafo has a great breakdown <a href="https://mrob.com/pub/comp/xmorphia/" target="_blank">**here**</a> that explains how the system works and what each constant signifies. He also has an excellent compilation of example constants and classified outputs <a href="https://mrob.com/pub/comp/xmorphia/pearson-classes.html" target="_blank">**here**</a>.  See rxn_dfsn_gs.m for the full documentation and more examples.

## Examples 

The above sample was created with the following function call:
```MATLAB
rxn_dfsn_gs(0.04,0.0635,'FrameSpacing',40);
```

An exteremly wide variety of patterns and motion can be generated by just altering the initial state and the f and k constants. The below sample was created from the following call:
```MATLAB
rxn_dfsn_gs(0.01,0.041,'InitState','wavefront');
```

![Sample GIF 2](https://github.com/nsbalbi/Gray-Scott-Reaction-Diffusion/blob/master/Sample%20GIF%202.gif)

Experiment Away!

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)
