# LEAF-Self Organizing Map-Shiny
In this shiny app you can have a more detailed overview of the putative effectors predicted by LEAF at a feature-by-feature scale.
In the last tab of the app it is possible to retrieve information on proteins belonging to one or more cell on the map.

## Required packages
- shiny (used here v1.7.5)
- DT (used here v0.29)
- readr (used here v2.1.4)
- aweSOM (used here v1.3)
- bslib (used here v0.5.1)
  
## Usage
**Command line** 
```
R -e "shiny::runApp('./')"
```
The app will directly read, as inputs, the 2 files contained in ./app_LEAF/tables directory 

**Using the provided singularity container**
```
singularity exec -B binding/dirs ./LEAF1.0_SOM.simg R -e "shiny::runApp('./')"
```
**Copy and paste the URL written on terminal after the application opening**
