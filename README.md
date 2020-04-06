# hydraulic_wetro

## Hydraulic modelling of the flood and rainfall demonstrator Wetropolis

This repository contains the source code and documentation on the numerical modelling of Wetropolis, and builds upon the design model of Bokhove et al. (2020). See also Onno Bokhove's [Wetropolis Github page](https://github.com/obokhove/wetropolis20162020/) which tells the story of Wetropolis from its inception to the present day.

### Background
Urban flooding is a major hazard worldwide, brought about by intense rainfall and often exacerbated by the built environment. The tabletop flood-demonstrator Wetropolis illustrates in an idealised modelling environment how extreme hydroclimatic events can cause flooding of a city due to peaks in river levels and groundwater following intense rainfall. It aims to conceptualise the science of flooding in a way that is accessible to and directly engages the public and also provides a scientific testbed for flood modelling, control and mitigation, and data assimilation. As such, it is useful to the scientist, industrial practitioner, and general public. 

Physically, it comprises a winding river channel with parallel canal, a reservoir for water storage, a porous flow cell (analogous to a moor) with observable groundwater flow, and random rainfall, which may or may not lead to flooding in the idealised urban area of Wetropolis. The main river channel has three $\pi$--degree bends and one $(\pi/4)$--degree bend and is fed by water flowing into the domain at an upstream entry and leaving the domain at the downstream exit. The river bed is sloping down (uniformly with gradient 1 in 100); the river cross-sectional area is rectangular and uniform, and flanked on one side by a sloping flood plain outside of the urban area. Through the urban area, the rectangular channel is flanked on both sides by flat rectangular plains of higher elevation than the regular river channel, i.e., the cross-sectional area is T-shaped. Water enters the main channel in three places: (i) the upstream inflow, generally kept constant; (ii) overflow of a groundwater cell (or "moor") with porous material and fed by random daily rainfall; and (iii) overflow from a reservoir, also fed by random daily rainfall. The two overflows can be placed in three different spots along the river: upstream, midstream or downstream just before the city plain. The set-up is displayed in plan-view below.

INSERT FIG.

Rainfall is supplied randomly in space at four locations (reservoir, moor, reservoir and moor, or nowhere) and randomly in time at four rainfall amounts (1s, 2s, 4s, or 9s) during a 10s Wetropolis day (wd) via two skew-symmetric discrete probability distributions. The joint probabilities (rain amount times rain location) are determined daily as one of 16 possible outcomes from two asymmetric Galton boards, in which steel balls fall down every wd and according to the (imposed) discrete probability distributions. The most extreme daily rainfall event thus involves 9s rainfall on both moor and reservoir. Its design is based on simulations of a one-dimensional numerical model of the dynamics that uses a kinematic approximation to describe the river flow and a depth-averaged nonlinear diffusion equation for the groundwater cell; a stochastic rainfall generator determines the amount and location of rain per wd.  In order to create an extreme flood event in Wetropolis once every 5 to 10 minutes on average instead of, say, once every $100$ to $200$ years on average (as in reality), this preliminary modelling determined the length of the Wetropolis day to be 10s. Thus, Wetropolis is able to demonstrate random extreme rainfall and flood events in a physical model on reduced spatial and temporal scales (see Bokhove et al. (2020) for more details).

### Goal
In Bpkhove et al. (2020), a numerical model (based on the equations for open channel flow under the kinematic assumption) is used to determine the relevant time and length scales prior to its construction as a physical model. That is, it is a crude and inexpensive model suitable for design purposes. **This page tracks the further development of the hydrodynamic modelling, both mathematically and numerically, with a view to conducting Wetropolis-inspired experiments in flood mitigation and control.**




## References
* Bokhove, 0., Hicks, T., Zweers, W. and Kent, T. (2019): Wetropolis extreme rainfall and flood demonstrator: from mathematical design to outreach. *Accepted: HESS*. Preprint on [EarthArXiv](https://eartharxiv.org/59ymk/) and open review/discussion at [HESSD](https://www.hydrol-earth-syst-sci-discuss.net/hess-2019-191/).

* ...
----

## Getting started
### Versions
All of the source code is written in ...
