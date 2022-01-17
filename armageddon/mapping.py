import folium
import numpy as np
from folium.plugins import HeatMap, MarkerCluster
import pandas as pd
from math import sin, cos, acos, asin, atan2, radians, degrees


def plot_circle(lat, lon, radius, map=None, **kwargs):
    """
    Plot a circle on a map (creating a new folium map instance if necessary).

    Parameters
    ----------
    lat: float
        latitude of circle to plot (degrees)
    lon: float
        longitude of circle to plot (degrees)
    radius: arraylike, float
        List of distances specifying the radius of circle(s) to plot (m)
    map: folium.Map
        existing map object

    Returns
    -------
    Folium map object

    Examples
    --------
    >>> import folium
    >>> armageddon.plot_circle(52.79, -2.95, 1e3, map=None)
    """
    # If the radius is int or float, change radius to list
    if isinstance(radius, (int, float)):
        radius = [radius]

    # If a map is not given, create a map
    if not map:
        map = folium.Map(location=[lat, lon], control_scale=True)

    # Decide colors which are used for showing damage zone circles.
    # zone1: purple, zone2: red, zone3: orange, zone4: yellow
    colors = ['#9370DB', '#DC143C', '#FF8000', '#FFFF00']

    # Plot color cicles with starting from zone1
    # To do so, sort the list of radius to fit zone number = color index+1.
    for i, rad in enumerate(sorted(radius, reverse=True)):
        folium.Circle([lat, lon], rad, fill=True,
                      fillOpacity=1., color=colors[i],
                      **kwargs).add_to(map)

    return map


def latlon_to_xyz(lat, lon):
    """Change lattitude and longitude into the rectangular coordinate system.
    The equatorial plane is the xy plane,
    and the axis of rotation is the z axis.

    Parameters
    ----------
    lat: float
        latitude(degree)
    lon: float
        longitude(degree)
    rlat: float
        latitude(rad)
    rlon: float
        longitude(rad)

    Returns
    ---------
    float
        Points on the rectangular coordinate system.
    """
    # Change degrees to radians
    rlat, rlon = radians(lat), radians(lon)

    return cos(rlat) * cos(rlon), cos(rlat) * sin(rlon), sin(rlat)


def xyz_to_latlon(x, y, z):
    """Change coodinate from xyz coordinate system to
    latitude & longitude coordinates.

    Parameter
    ----------
    x: float
       x coordinate of Equatorial plane
    y: float
        y coodinate of Equatorial plane
    z: float
        z coodinate of Arctic direction

    Returns
    ---------
    float
        Points on the earth surface(degree)
    """

    rlat = asin(z)
    coslat = cos(rlat)
    return degrees(rlat), degrees(atan2(y / coslat, x / coslat))


def halfway_on_sphere(lat, lon, elat, elon, z):
    """
    Calculate a point on the great circle rout of asteroid. If z= 0.5,
    the return shows lat & lon of harfway point.

    Parameter
    ---------
    lat: float
        latitude of zero point
    lon: float
        longitude of zero point
    elat: float
        latitude of entry point
    elon: float
        longitude of entry point
    z: float
        calculation point between entry and zero point.

    Return
    --------
    list
        latitude and longitude of interval point.
    """

    # Cange lattitude & longitude to xyz coodinate
    xyz0, xyz1 = latlon_to_xyz(lat, lon), latlon_to_xyz(elat, elon)

    # Calculate a distance between entry point and zero point.
    theta = acos(sum(x * y for x, y in zip(xyz0, xyz1)))
    v0 = sin(theta * (1 - z)) / sin(theta)
    v1 = sin(theta * z) / sin(theta)

    # Calculate latitude a& longitude of interval point.
    interval_lat, interval_lon = xyz_to_latlon(
        *(x * v0 + y * v1 for x, y in zip(xyz0, xyz1)))

    return [interval_lat, interval_lon]


def plot_line(lat, lon, elat, elon, map=None, n=100):
    """
    Plot a black lineconnecting entry point and zero point.
    Parameters
    ----------
    lat: float
        latitude of circle to plot (degrees)
    lon: float
        longitude of circle to plot (degrees)
    elat: float
        latitude of entry point (degrees)
    elon: float
        longitude of entry point(degrees)
    n: int
        number of rute divisions.Default value is 100.
    map: folium.Map
        existing map object

    Returns
    -------
    Folium map object

    Examples
    --------
    >>> plot_line(52.79, -2.95, 53.48, -2.24 , map=None)
    """

    # If a map is not given, create a map.
    if not map:
        map = folium.Map(location=[lat, lon], control_scale=True)

    Harf = []

    # Calculate points of intervals between entry and zero pints.
    for i in range(n):
        Intervals = halfway_on_sphere(lat, lon, elat, elon, z=i / (n + 1))
        Harf.append(Intervals)

    # Make a list of plotting points : [entry point, harf point, zero point]
    points = [[lat, lon], *Harf, [elat, elon]]

    # Plotting a line on map.
    folium.PolyLine(points,
                    color="black", weight=2.5, opacity=1).add_to(map)

    return map


def get_lat_long_of_postcodes(postcodes, sector=False):
    """
    Return location(latitude,longitude) of a list of postcode units or sectors.

    Parameters
    ----------
    postcodes : list of lists
        list of postcode units or postcode sectors
    sector : bool, optional
        if true return populations for postcode sectors,
        otherwise postcode units

    Returns
    -------
    list of lists
    Contains the latitude,longitude of input postcode units or sectors

    Examples
    --------
    >>> get_lat_log_of_postcode([['SW7 2AZ','SW7 2BT','SW7 2BU','SW7 2DD']])
    >>> get_lat_log_of_postcode([['SW7 2']], True)
    """

    # Get postcodes from csv
    postcodes_pd = pd.read_csv('./armageddon/resources/full_postcodes.csv')

    # Modify postcodes to no spaces for processing
    postcodes = [[x2.replace(" ", "_") for x2 in x1] for x1 in postcodes]
    postcodes_pd['Postcode'] = postcodes_pd['Postcode'].str.replace(" ", "_")

    # If sector flag is True―taking average of unit locations
    if sector:
        postcodes_pd = postcodes_pd.groupby(
            postcodes_pd['Postcode'].str.slice(stop=5),
            as_index=True).mean().reset_index()

    # Select postcodes
    select_postcodes = postcodes_pd[
        postcodes_pd['Postcode'].isin(postcodes[0])][['Latitude', 'Longitude']]

    return select_postcodes.values.tolist()


def heat_map_layer(locations, weights, map=None, radius=25):
    """
    Return heat map layer for follium map from
    a list of locations and a list of weights

    Parameters
    ----------
    locations : list of lists
        list of latitutde and longitude coordinates
        corresponding to postcode units or postcode sectors
    weights : list of lists, array-like
        list of weights to be plotted at locations

    Returns
    -------
    Follium map

    Examples
    --------
    >>> locations = get_lat_long_of_postcodes(postcodes, sector=False)
    >>> weights = [['10000', '20000', '30000', '40000']]
    >>> heat_map_layer(locations, weights, map = None, radius = 25)
    """

    # Calculate an average of latitude and longitude of given locations
    Avr_location = np.average(locations, axis=0)

    # If a map is not given, create a map
    if not map:
        map = folium.Map(location=Avr_location, control_scale=True)

    # Creating copy of locations
    combo = locations.copy()

    # Appending weight to the third column of combo.
    for i, a in enumerate(combo):
        a.append(float(weights[0][i]))

    # Initialize Follium HeatMap instance
    heat_map = HeatMap(combo, name=None, min_opacity=0.5,
                       max_zoom=18, radius=radius, blur=15, gradient=None,
                       overlay=False, control=True, show=True)
    heat_map.add_to(map)

    return map


def plot_marker(lat, lon, popup=None, map=None, **kwargs):
    """
    Plot a point on a map (creating a new folium map instance if necessary).

    Parameters
    ----------

    lat: float
        latitude of point to plot (degrees)
    lon: float
        longitude of point to plot (degrees)
    popup: str
        will plot a string label at point
    map: folium.Map
        existing map object

    Returns
    -------

    Folium map object

    Examples
    --------

    >>> import folium
    >>> armageddon.plot_point(52.79, -2.95, 1e3, map=None)
    """

    if popup is not None:
        if isinstance(popup, (str)) is False:
            popup = None

    if not map:
        map = folium.Map(location=[lat, lon], control_scale=True)

    folium.map.Marker(location=[lat, lon], popup=popup,
                      tooltip=None, icon=None,
                      draggable=False, **kwargs).add_to(map)

    return map


def plot_multiple_markers(locations, popups=None, map=None):
    """
    Return heat cluster of markers for follium map from
    a list of locations

    Parameters
    ----------
    locations : list of lists
        list of latitutde and longitude coordinates
        corresponding to postcode units or postcode sectors
    popup: list of str
        will plot a string label at points
    map: folium.Map
        existing map object

    Returns
    -------
    Follium map

    Examples
    --------

    >>> locations = get_lat_long_of_postcodes(postcodes, sector=False)
    >>> plot_multiple_markers(locations, popups= None, map = None)
    """

    Avr_location = np.average(locations, axis=0)

    if not map:
        map = folium.Map(location=Avr_location, control_scale=True)

    map = MarkerCluster(locations=locations, popups=popups,
                        icons=None, name='Location Markers',
                        overlay=True, control=True,
                        show=True, icon_create_function=None,
                        options=None).add_to(map)

    return map
