import armageddon

###################
# Airburst Solver #
###################

# Initialise the Planet class
earth = armageddon.Planet()
print("Here in A")
# Solve the atmospheric entry problem (for something similar to Chelyabinsk).
result = earth.solve_atmospheric_entry(radius=10, angle=20,
                                       strength=1e6, density=3000,
                                       velocity=19e3)
print("Here in B")
# Calculate the kinetic energy lost per unit altitude and add it
# as a column to the result dataframe
result = earth.calculate_energy(result)
print("Here in C")
# Determine the outcomes of the impact event
outcome = earth.analyse_outcome(result)
print("Here in D")
#################
# Damage Mapper #
#################

# Calculate the blast location and damage radius for several pressure levels
blast_lat, blast_lon, damage_rad = armageddon.damage_zones(outcome,
                                                           lat=51.2, lon=0.7,
                                                           bearing=-35.,
                                                           pressures=[1e3,
                                                                      3.5e3,
                                                                      27e3,
                                                                      43e3])
print("Here in E")
# Plot a circle to show the limit of the lowest damage level
damage_map = armageddon.plot_circle(blast_lat, blast_lon, damage_rad[0])
damage_map.save("damage_map.html")
print("Here in F")
# The PostcodeLocator tool
locator = armageddon.PostcodeLocator()
print("Here in G")
# Find the postcodes in the damage radii
postcodes = locator.get_postcodes_by_radius((blast_lat, blast_lon),
                                            radii=damage_rad)
print("Here in H")
# Find the population in each postcode
population = locator.get_population_of_postcode(postcodes)
print("Here in I")
# Alternatively find the postcode sectors in the damage radii,\
# and populations of the sectors
sectors = locator.get_postcodes_by_radius((blast_lat, blast_lon),
                                          radii=damage_rad, sector=True)
print("Here in J")
population_sector = locator.get_population_of_postcode(sectors, sector=True)
print("Here in K")
