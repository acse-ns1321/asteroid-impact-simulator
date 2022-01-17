import numpy as np
import pandas as pd
from scipy.optimize import fsolve
import armageddon


def Calculate_r(x, *data):
    pressure = data
    return 3.14e+11 * x**(-1.3) + 1.8e+7 * x**(-0.565) - pressure


def damage_zones(outcome, lat, lon, bearing, pressures):
    """
    Calculate the latitude and longitude of the surface zero location and the
    list of airblast damage radii (m) for a given impact scenario.

    Parameters
    ----------

    outcome: Dict
        the outcome dictionary from an impact scenario
    lat: float
        latitude of the meteoroid entry point (degrees)
    lon: float
        longitude of the meteoroid entry point (degrees)
    bearing: float
        Bearing (azimuth) relative to north of meteoroid trajectory (degrees)
    pressures: float, arraylike
        List of threshold pressures to define airblast damage levels

    Returns
    -------

    blat: float
        latitude of the surface zero point (degrees)
    blon: float
        longitude of the surface zero point (degrees)
    damrad: arraylike, float
        List of distances specifying the blast radii
        for the input damage levels

    Examples
    --------

    >>> import armageddon
    >>> outcome = {'burst_altitude': 8e3, 'burst_energy': 7e3,
                   'burst_distance': 90e3, 'burst_peak_dedz': 1e3,
                   'outcome': 'Airburst'}
    >>> armageddon.damage_zones(outcome, 52.79, -2.95,
    135, pressures=[1e3, 3.5e3, 27e3, 43e3])
    """

    # Replace this code with your own. For demonstration we return lat, lon
    # and 1000 m
    R_p = 6371000
    damrad = []
    Ek = outcome['burst_energy']
    r = outcome['burst_distance']
    Zb = outcome['burst_altitude']
    lat = np.radians(lat)
    bearing = np.radians(bearing)
    sinphi_2 = np.sin(lat) * np.cos(r / R_p) + np.cos(lat) * \
        np.sin(r / R_p) * np.cos(bearing)
    blat = np.rad2deg(np.arcsin(sinphi_2))
    blon = np.rad2deg(np.arctan(((np.sin(bearing) * np.sin(r / R_p) *
                      np.cos(lat)) / (np.cos(r / R_p) - np.sin(lat)
                                      * sinphi_2)))) + lon

    for i in range(len(pressures)):
        data = pressures[i]
        a = fsolve(Calculate_r, [1], args=data)
        if a * Ek ** (2 / 3) - Zb ** 2 > 0:
            r_exposure = (a * Ek ** (2 / 3) - Zb ** 2) ** (1 / 2)
            damrad.append(float(r_exposure))
        else:
            damrad.append(0)

    blat = float(blat)
    blon = float(blon)
    return blat, blon, damrad


fiducial_means = {'radius': 10, 'angle': 20, 'strength': 1e6,
                  'density': 3000, 'velocity': 19e3,
                  'lat': 51.5, 'lon': 1.5, 'bearing': -45.}
fiducial_stdevs = {'radius': 1, 'angle': 1, 'strength': 5e5,
                   'density': 500, 'velocity': 1e3,
                   'lat': 0.025, 'lon': 0.025, 'bearing': 0.5}


def impact_risk(planet, means=fiducial_means, stdevs=fiducial_stdevs,
                pressure=27.e3, nsamples=100, sector=True):
    """
    Perform an uncertainty analysis to calculate the risk for each affected
    UK postcode or postcode sector

    Parameters
    ----------
    planet: armageddon.Planet instance
        The Planet instance from which to solve the atmospheric entry

    means: dict
        A dictionary of mean input values for the uncertainty analysis. This
        should include values for ``radius``, ``angle``, ``strength``,
        ``density``, ``velocity``, ``lat``, ``lon`` and ``bearing``

    stdevs: dict
        A dictionary of standard deviations for each input value. This
        should include values for ``radius``, ``angle``, ``strength``,
        ``density``, ``velocity``, ``lat``, ``lon`` and ``bearing``

    pressure: float
        The pressure at which to calculate the damage zone for each impact

    nsamples: int
        The number of iterations to perform in the uncertainty analysis

    sector: logical, optional
        If True (default) calculate the risk for postcode sectors, otherwise
        calculate the risk for postcodes

    Returns
    -------
    risk: DataFrame
        A pandas DataFrame with columns for postcode (or postcode sector) and
        the associated risk. These should be called ``postcode`` or ``sector``,
        and ``risk``.
    """
    gaussian_normalValue = {}
    for index in means:
        gaussian_normalValue[index] = np.random.normal(
            means[index], stdevs[index], nsamples)

    random_planets = []
    risk_dic = {}
    ori_lats = []
    ori_lons = []
    sur_lats = []
    sur_lons = []
    damards = []

    for i in range(nsamples):
        random_planets.append([gaussian_normalValue['radius'][i],
                               gaussian_normalValue['angle'][i],
                               gaussian_normalValue['strength'][i],
                               gaussian_normalValue['density'][i],
                               gaussian_normalValue['velocity'][i],
                               gaussian_normalValue['lat'][i],
                               gaussian_normalValue['lon'][i],
                               gaussian_normalValue['bearing'][i]])
        ori_lats.append(gaussian_normalValue['lat'][i])
        ori_lons.append(gaussian_normalValue['lon'][i])

        solvers = planet.solve_atmospheric_entry(radius=random_planets[i][0],
                                                 velocity=random_planets[i][4],
                                                 density=random_planets[i][3],
                                                 strength=random_planets[i][2],
                                                 angle=random_planets[i][1])
        solvers = planet.calculate_energy(solvers)
        outcome = planet.analyse_outcome(solvers)
        b_lat, b_lon, dam_rad = damage_zones(
            outcome, lat=random_planets[i][5], lon=random_planets[i][6],
            bearing=random_planets[i][7], pressures=[pressure])
        sur_lats.append(b_lat)
        sur_lons.append(b_lon)
        damards.append(dam_rad)
        locator = armageddon.locator.PostcodeLocator()
        if sector:
            sectors = locator.get_postcodes_by_radius(
                (b_lat, b_lon), radii=dam_rad, sector=True)
            pop_sector = locator.get_population_of_postcode(
                sectors, sector=True)
            try:
                for i in range(len(sectors[0])):
                    if sectors[0][i] not in risk_dic:
                        risk_dic[sectors[0][i]] = pop_sector[0][i]
                    else:
                        risk_dic[sectors[0][i]] += pop_sector[0][i]
            except BaseException:
                pass
        else:
            postcodes = locator.get_postcodes_by_radius(
                (b_lat, b_lon), radii=dam_rad)
            pop_postcodes = locator.get_population_of_postcode(postcodes)
            try:
                for i in range(len(postcodes[0])):
                    if postcodes[0][i] not in risk_dic:
                        risk_dic[postcodes[0][i]] = pop_postcodes[0][i]
                    else:
                        risk_dic[postcodes[0][i]] += pop_postcodes[0][i]
            except BaseException:
                pass

    key = []
    risk = []

    for i in risk_dic:
        risk_dic[i] = risk_dic[i] / nsamples
        key.append(i)
        risk.append(risk_dic[i])
    if sector is True:
        return pd.DataFrame({'sector': key, 'risk': risk}
                            ), ori_lats, ori_lons, sur_lats, sur_lons, damards
    else:
        return pd.DataFrame({'postcode': key, 'risk': risk}
                            ), ori_lats, ori_lons, sur_lats, sur_lons, damards
