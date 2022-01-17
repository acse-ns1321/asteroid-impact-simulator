# import armageddon

# ###################
# # Airburst Solver #
# ###################

# # Initialise the Planet class
# earth = armageddon.Planet()

# # Solve the atmospheric entry problem (for current UK breaking
#  news headline).
# result = earth.solve_atmospheric_entry(radius=65, angle=15,
#                                        strength=3e6, density=3200,
#                                        velocity=19e3)

# # Calculate the kinetic energy lost per unit altitude and add it
# # as a column to the result dataframe
# result = earth.calculate_energy(result)

# # Determine the outcomes of the impact event
# outcome = earth.analyse_outcome(result)

# # Calculate the blast location and damage radius for several pressure levels
# blast_lat, blast_lon, damage_rad = armageddon.damage_zones(outcome,
#                                                            lat=53.439,
# lon=-2.966,
#                                                            bearing=159.,
#                                                            pressures=[1e3,
#                                                                       3.5e3,
#                                                                       27e3,
#                                                                       43e3])
# print(f'blast center: ({blast_lat, blast_lon})')
# print(f'blast damage - for pressure of 1   km: ({damage_rad[0]})')
# print(f'blast damage - for pressure of 3.5 km: ({damage_rad[1]})')
# print(f'blast damage - for pressure of 27  km: ({damage_rad[2]})')
# print(f'blast damage - for pressure of 43  km: ({damage_rad[3]})')

# # Plots a circle defining the damage zone(s) around the blast center
# damage_map = armageddon.plot_circle(blast_lat, blast_lon, damage_rad)

# # add marker at blast center
# armageddon.plot_marker(blast_lat, blast_lon,
#                        popup="Event: {}, Altitude: {}, Energy: {}".format
# (outcome['outcome'],
#
# outcome['burst_altitude'],
#
# outcome['burst_energy'])

# # damage_map.save("damage_map.html")


# # Find postcodes within blast damage radii
# locator=armageddon.PostcodeLocator()

# postcodes=locator.get_postcodes_by_radius((blast_lat, blast_lon),
#                                             radii=damage_rad)


# # Visualize postcodes
# post_loc=armageddon.mapping.get_lat_long_of_postcodes(
#     [postcodes[0]], sector=True)
# post_map=armageddon.plot_circle(blast_lat, blast_lon, damage_rad)
# armageddon.plot_marker(blast_lat, blast_lon, map=post_map)
# armageddon.plot_multiple_markers(post_loc, popups=postcodes[0], map=post_map)

# # post_map.save("impacted_postcodes_map.html")

# # Calcuate population of effected postcodes
# population=locator.get_population_of_postcode(postcodes)


# # Visualize population centers with a HeatMap scaled by population density
# pop_sector_heatmap=armageddon.plot_circle(blast_lat, blast_lon, damage_rad)
# armageddon.plot_marker(blast_lat, blast_lon, map=pop_sector_heatmap)
# armageddon.mapping.heat_map_layer(
#     sector_loc, sector_population, map=pop_sector_heatmap, radius=50)

# # pop_sector_heatmap.save("impacted_population_map.html")

# # Calculate risk impact
# # example for possible damage zones for 10 scenarios to assess uncertainty
# earth=armageddon.Planet()
# result=armageddon.impact_risk(earth, pressure=1.e3, nsamples=10, sector=True)


# # Vizualize risk impact
# risk_loc=armageddon.mapping.get_lat_long_of_postcodes(
#     [result[0]["sector"].tolist()], sector=True)
# impact_scenario_map=armageddon.plot_marker(53.439, -2.966, map=None)
# armageddon.plot_multiple_markers(
#     risk_loc, popups=None, map=impact_scenario_map)

# for i in range(len(result)):

#     # Plots circle for 1e3 pressure damage radii for scenarios
# predicted to have impact
#     armageddon.mapping.plot_circle(
#         result[3][i], result[4][i], result[5][i][0], map=impact_scenario_map)
# armageddon.mapping.heat_map_layer(
#     risk_loc, [result[0]["risk"].tolist()], map=impact_scenario_map,
# radius=25)

# # impact_scenario_map .save("risk_impact_map.html")
