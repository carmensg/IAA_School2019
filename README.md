<img src="https://github.com/carmensg/IAA_School2019/blob/master/images/IAA-CSIC_School.png" width="900" align="left">




# Instituto de Astrofísica de Andalucía (CSIC) #
## 1st School on [Statistics, Data Mining and Machine Learning](https://www.granadacongresos.com/sostat) ##
## November 4th-7th, 2019 - Granada, Spain ##


It is an introductory school oriented at students and researchers of all levels. It is the first of its kind and will be followed by a School on Advanced Data Mining, Machine Learning and Artificial Intelligence in Fall 2020. The cycle of introductory and advanced schools will repeat in 2021 and 2022.
 
The school will have a strong hands-on aspect and the participants will repeat and apply the learned lessons in practical exercises. It is therefore indispensable that each participant bring their own laptop with the required software pre-installed.

This repository holds the tutorials and scripts used to illustrate the theoretical concepts during the School. 

### Tutors ###

* Matteo Bachetti (INAF, Osservatorio Astronomico di Cagliari, Italy)
* Gwendolyn Eadie (University of Toronto, Canada)
* Carmen Sánchez Gil (Universidad de Cádiz, Spain)
* Željko Ivezić (University of Washington, USA)
* Javier Pascual Granado (IAA-CSIC, Spain)

The tutorial is composed of theoretical and practical modules. 

The examples will be demonstrated in R and Python but _familiarity with the language is not a requirement_.

### SOFTWARE & INSTALLATION INSTRUCTIONS ###

The workshop and hands-on sessions will be based on **R** and **python**. 
In order to optimize the time spent in the examples and practical applications, we advise the participants to get the following software up and running in advance.

* [R](https://www.r-project.org)
* [Rstudio](https://www.rstudio.com)
* [python](https://www.python.org), with the following packages: [numpy,scipy, scikit-learn](https://scipy.org/install.html), [matplotlib](https://matplotlib.org/3.1.1/users/installing.html)
* [ipython notebook](https://ipython.org/notebook.html)
* [Jupyter Notebook](https://jupyter.org)
* [astroML](https://www.astroml.org)  (scroll down to "2. Installation of astroML")

We strongly recommend the all-in-one scientific Python/R installer [Anaconda](https://www.anaconda.com/distribution/). All the previous software and packages can be directly and efficiently installed through this platform.
It will be also useful the use of [Anaconda Navigator](http://docs.anaconda.com/anaconda/navigator/). 

More details on installation: [here](https://github.com/carmensg/IAA_School2019/blob/master/SOFTWARE_INSTALLATION.md)


**Getting the scripts**

In order to avoid problems with file paths it is advisable to clone this repository and work within it.  

To do so choose go to the command line and navigate to a location where you would like to work. Then type:

    git clone git@github.com:carmensg/IAA_School2019.git


This should be enough to get you ready for the examples we will be working on. 

If you do not have git installed you can get it typing:

    sudo apt-get install git


**Staying up to date**

This is a work in progress and it will be continuously updated so erros can be fixed and complementary material can be added. 

It is advisable to make sure you have the latest version before start working in you local directory. To do so, in the command line navigate to your copy of this repository and type:

    git pull 

Then you are certain to get all the bug fixes and improvements available.

Have fun!

### WORKSHOP MATERIAL ###

### Day 1  - Classical and Bayesian statistical inference (Gwen Eadie) ###


The slides used in this tutorial are available [here](https://github.com/carmensg/IAA_School2019/tree/master/lectures).

[Tutorial Examples](https://github.com/carmensg/IAA_School2019/tree/master/lectures)



### Day 2 - Data modeling and parameter estimation (Carmen Sánchez Gil) ###


The slides used in this tutorial are available [here](https://github.com/carmensg/IAA_School2019/tree/master/lectures).

[Tutorial Examples](https://github.com/carmensg/IAA_School2019/tree/master/lectures)



### Day 3 - Introduction to data mining: Searching for structure in data (Željko Ivezić) ###

For Day 3, you need to install [astroML](https://www.astroml.org)

To test your installation, please download and run this 
[testing notebook](https://github.com/carmensg/IAA_School2019/blob/master/lectures/Day3-ZeljkoIvezic/notebooks/astroMLtesting.ipynb) 

There are four jupyter python notebooks and two pdf lecture files that we will use in this section.

Please download them prior to class and run all the cells in notebooks to get some slow computations done in advance. 

The easiest way to clone these notebooks, and supporting files in subdirectory "figures", is 
to clone the entire IAA_School2019 repository (e.g. >git clone git@github.com:carmensg/IAA_School2019.git)

[1. Density Estimation](https://github.com/carmensg/IAA_School2019/tree/master/lectures/Day3-ZeljkoIvezic/notebooks/density_estimation.ipynb)

[2. Clustering](https://github.com/carmensg/IAA_School2019/tree/master/lectures/Day3-ZeljkoIvezic/notebooks/clustering.ipynb)

[3. Classification](https://github.com/carmensg/IAA_School2019/tree/master/lectures/Day3-ZeljkoIvezic/notebooks/classification.ipynb)

[4. Dimensionality Reduction](https://github.com/carmensg/IAA_School2019/tree/master/lectures/Day3-ZeljkoIvezic/notebooks/dimensionality_reduction.ipynb)

### Day 4 - Time series analysis (Matteo Bacchetti, Javier Pascual Granado) ###

The slides used in this tutorial are available [here](https://github.com/carmensg/IAA_School2019/tree/master/lectures).

[Tutorial Examples](https://github.com/carmensg/IAA_School2019/tree/master/lectures)



### References ###






