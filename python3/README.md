# hydraulic_wetro: python3 code

![Wetro py3 dashboard](configs/config#2/t990.png)

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
```main_wetro_tests.py```         | Main run script for test cases
```main_wetro_fullsyst.py```    | Main run script for full system: 1st implementation
```init_cond.py```              | Initial condition functions
```flux_function.py```          | Numerical flux calculation for space discretisation
```cross_sections.py```         | Compute cross-sections h(A,s) and A(h,s)
```plot_wetro.py```             | Very basic plotting routine for data generated in ```main_wetro_tests.py```
```plot_wetro_anim.py```             | Plotting routine (animation) for data generated in ```main_wetro_tests.py```
```plot_wetro_amp.py```             | Plotting routine (interactive animation) for data generated in ```main_wetro_tests.py```
```random_rainfall.py```     | Test script for generating and plotting random rainfall
```/configs```                  | Dir for config files:
```/configs/config#0.py``` etc. | Parameters for steady state test case
```/configs/config#1.py``` etc. | Parameters for flood wave test case
```/configs/config#2.py``` etc. | Parameters for full system


## Runninng the code

### Steady state test case

* run ```main_wetro_tests.py``` with ```config#0.py``` specified as ```spec```. This saves data in the dir  ```/configs/config#0```.
* run plotting routine ```plot_wetro_amp.py```, again with ```config#0.py``` specified as ```spec```.

### Flood wave test case

* run ```main_wetro_tests.py``` with ```config#1.py``` specified as ```spec```. This saves data in the dir  ```/configs/config#1```.
* run plotting routine ```plot_wetro_amp.py``` again with ```config#1.py``` specified as ```spec```.

### Full system

...
