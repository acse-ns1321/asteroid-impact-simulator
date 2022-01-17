# Asteroid Impact Simulation Tool
![alt text](https://www.freethink.com/wp-content/uploads/2021/04/asteroid-impact-simulation_banner.jpg?resize=1800,675)
This software tool developed by Team Ceres. This had the ability to solve the system of differential equations describing meteoroid entry and compute the burst altitude, burst energy and horizontal path length from the entry point. This tool can also take these outputs and a location in the UK and determine the predicted extent of airblast damage on the ground and the postcodes and population affected.
The output is in the form of damage zone maps that plot the predicted risk sectors and the population affected by the simulated asteroid entry

## Installation

To install the module and any pre-requisites for the tool to function
```
pip install -r requirements.txt
```  

## Downloading postcode data 

To download the postcode data required for tool simulation
```
python download_data.py
```

## Automated testing

To run the pytest test suite
```
python -m pytest armageddon
```


## How it Works

Tool has different parts

- The solver - that calculates burst altitude, burst energy and horizontal path length from the entry point
- The locator - that inputs entry point and burst radii to caluclate postcodes and populations affected
- The damage radii calculator - caclulated the risk radii for the given postcodes and entry point
- The damage risk calculator - caclulated the risk probabilities  for the given postcodes and entry point
- The mappers - this maps all the damage risk on the postcodes, the radii of damage and the trajectory of the asteriod on the map


## Example of a Usage in Python

For example usage see `example.py` 
```
python example.py
```

## Example usage from the Scenario

For the scenario analysis see `scenario.py`

```
python scenario.py

```

## Authors

The contributors to this software tool include
        - Yin, Shiqi
        - Yao, Tianshun
        - NAKAMURA, Yuna
        - Sundararajan, Niranjana
        - YOU, ZHEXIN
        - Johnson, Eleda
        - Tang, Jieyi
        - Cheng, xiaoyuan
## Version History

[Version Files](Version.md)
## License

This project is licensed under the MIT License  - see the [LICENSE.md](license.md) file for details

