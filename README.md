# IndTB_dx
Portfolio modelling for TB epidemiology in India

## Project structure
- **/docs**: documentation for model equations and output figures/tables
- **/data**: data inputs and calibration targets
  - **/targets.csv**: all available calibration targets
  - **/targets_sel.csv**: selected calibration targets
- **/pars**: parametrised inputs, including demography inputs and care-cascade
- **/sim**: main dynamics model codes in python
- Scripts (initial number as the exec order):
  - **/scripts_py**: simulation/scripts for python-version model
  - **/scripts_r**: scripts for data processing/visualisation/output tables
- **/R_ver**: a standalone version for running the model in R

## Interventions included

### Treatment regimen

- Increasing treatment initiation
- Increasing cure rate
- Reducing relapse

### Improving diagnosis

- Reducing failing diagnosis due to sample collection
- Increasing sensitivity/specificity
- Strengthening the link between diagnosis and treatment

### Vaccination

- Preventing from infection
- Preventing from disease progression
- Preventing form relapse

### Mass-screening

- Annual screening at large scales

## Notes

- Model for data processing are constructed using rstan. Please use the [Docker](docker-compose.yml) to rebuild the environment.
- The main version of the simulation model is the python one. All analyses are generated with the version
- A R version [model in Odin](model/tb.R) replicates the python version in the same functionality but less performance. 


## License
[License](LICENSE)



