# hydraulic_wetro: python3 code

<!--
---
## Contents

* [Introduction](#introduction)
  * [Motivation](#motivation)
  * [Description](#A-brief-description-of-Wetropolis)
  * [Taster](#taster)
  * [References](#references)
* [Getting started](#getting-started)
* [Code overview](#files-overview)
  * [MATLAB](#matlab)
  * [Python](#python)
* [Preliminary simulations](#preliminary-simulations)
---
 -->

## Files overview

File/dir name                   |  Summary
:--------------------------:|:--------------------------:
```main_wetro_run.py```         | Main run script for initial test case
```main_wetro_fullsyst.py```    | Main run script for full system: 1st implementation
```init_cond.py```              | Initial condition functions
```flux_function.py```          | Numerical flux calculation for space discretisation
```cross_sections.py```         | Compute cross-sections h(A,s) and A(h,s)
```plot_wetro.py```             | Basic plotting routine for data generated in main scripts
```plot_wetro_anim.py```             | Basic plotting routine for data generated in main scripts
```plot_wetro_amp.py```             | Basic plotting routine for data generated in main scripts
```random_rainfall.py```     | Test script for generating and plotting random rainfall
```/configs```                  | Dir for config files:
```/configs/config#0.py``` etc. | Parameters for steady state test case
```/configs/config#1.py``` etc. | Parameters for flood wave test case
```/configs/config#2.py``` etc. | Parameters for full system


## Runninng the code

### Steady state test case

* run ```main_wetro_run.py``` with ```config#0.py``` specified as ```spec``` in line 29. This saves data in the dir  ```/configs/config#0```.
* run plotting routine ```plot_wetro_amp.py``` again with ```config#0.py``` specified as ```spec```.

### Flood wave test case

* run ```main_wetro_run.py``` with ```config#1.py``` specified as ```spec``` in line 29. This saves data in the dir  ```/configs/config#1```.
* run plotting routine ```plot_wetro_amp.py``` again with ```config#1.py``` specified as ```spec```.

### Full system

...
