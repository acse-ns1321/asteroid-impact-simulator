"""Module dealing with postcode information."""
import numpy as np
import pandas as pd


def great_circle_distance(latlon1, latlon2):
    """
    Calculate the great circle distance (in metres) between pairs of
    points specified as latitude and longitude on a spherical Earth
    (with radius 6371 km).

    Parameters
    ----------

    latlon1: arraylike
        latitudes and longitudes of first point (as [n, 2] array for n points)
    latlon2: arraylike
        latitudes and longitudes of second point (as [m, 2] array for m points)

    Returns
    -------

    numpy.ndarray
        Distance in metres between each pair of points (as an n x m array)

    Examples
    --------

    >>> import numpy
    >>> fmt = lambda x: numpy.format_float_scientific(x, precision=3)}
    >>> with numpy.printoptions(formatter={'all', fmt}):
        print(great_circle_distance([[54.0, 0.0], [55, 0.0]], [55, 1.0]))
    [1.286e+05 6.378e+04]
    """
    # Radius of Earth in meters
    Rp = 6371000

    # Convert to np arrays and reshape to a 2D array
    point1 = np.array(latlon1) * np.pi/180
    point1 = np.reshape(point1, (-1, 2))

    point2 = np.array(latlon2) * np.pi/180
    point2 = np.reshape(point2, (-1, 2))

    # Initialize distance array
    distance = np.empty((len(point1), len(point2)), float)

    # Calculate distance using Vincenty Formula
    for i in range(0, len(point1)):
        for j in range(0, len(point2)):
            sum_1 = np.cos(point2[j][0]) \
                * np.sin(np.abs(point1[i][1] - point2[j][1]))
            sum_2 = np.cos(point1[i][0]) \
                * np.sin(point2[j][0]) - np.sin(point1[i][0]) \
                * np.cos(point2[j][0]) \
                * np.cos(np.abs(point1[i][1] - point2[j][1]))
            num = np.sqrt(sum_1 ** 2 + sum_2 ** 2)
            den = np.sin(point1[i][0]) \
                * np.sin(point2[j][0]) + np.cos(point1[i][0]) \
                * np.cos(point2[j][0]) \
                * np.cos(np.abs(point1[i][1] - point2[j][1]))
            ratio = num/den
            distance[i, j] = Rp * np.arctan(ratio)

    return distance


class PostcodeLocator(object):
    """Class to interact with a postcode database file."""

    def __init__(self, postcode_file='./armageddon/' +
                 'resources/full_postcodes.csv',
                 census_file='./armageddon/resources/' +
                 'population_by_postcode_sector.csv',
                 norm=great_circle_distance):
        """
        Parameters
        ----------

        postcode_file : str, optional
            Filename of a .csv file containing geographic
            location data for postcodes.

        census_file :  str, optional
            Filename of a .csv file containing census data by postcode sector.

        norm : function
            Python function defining the distance between
            points in latitude-longitude space.

        """
        self.norm = norm

        # Read data into pandas dataframe
        self.postcodes_data = pd.read_csv(postcode_file)
        self.population_data = pd.read_csv(census_file)

        # Convert into co-ordinates
        latutides = np.array(self.postcodes_data['Latitude'])\
            .reshape((self.postcodes_data['Longitude'].shape[0], 1))
        longitude = np.array(self.postcodes_data['Longitude'])\
            .reshape((self.postcodes_data['Longitude'].shape[0], 1))

        # Stack the co-ordinates
        self.latlon = np.hstack((latutides, longitude))

    def get_postcodes_by_radius(self, X, radii, sector=False):
        """
        Return (unit or sector) postcodes within specific distances of
        input location.

        Parameters
        ----------
        X : arraylike
            Latitude-longitude pair of centre location
        radii : arraylike
            array of radial distances from X
        sector : bool, optional
            if true return postcode sectors, otherwise postcode units

        Returns
        -------
        list of lists
            Contains the lists of postcodes closer than
            the elements of radii to the location X.


        Examples
        --------

        >>> locator = PostcodeLocator()
        >>> locator.get_postcodes_by_radius((51.4981, -0.1773), [0.13e3])
        >>> locator.get_postcodes_by_radius((51.4981, -0.1773),
                                            [0.4e3, 0.2e3], True)
        """
        # Initialise dataframe for postcode retrival
        latlon = pd.DataFrame(self.latlon)
        latlon['center_to_postcode'] = self.norm(self.latlon, X)
        latlon['Postcodes'] = self.postcodes_data['Postcode']
        result = []

        # Check distance from center for each postcode
        for r in radii:
            latlon['distance'] = latlon['center_to_postcode'] - r
            # Store impact postcode into new column
            impact_latlon = latlon.loc[latlon['distance'] < 0]
            impact_postcodes = impact_latlon['Postcodes']
            # If sector flag is true, calculate sector codes
            if sector:
                sector_codes = impact_postcodes.str[:-2]
                unique_sectors = sector_codes.unique().tolist()
                result.append(unique_sectors)
            # If sector flag is flase, calculate unit codes
            else:
                unique_postcodes = impact_postcodes.unique().tolist()
                result.append(unique_postcodes)

        return result

    def get_population_of_postcode(self, postcodes, sector=False):
        """
        Return populations of a list of postcode units or sectors.

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
            Contains the populations of input postcode units or sectors


        Examples
        --------

        >>> locator = PostcodeLocator()
        >>> locator.get_population_of_postcode([['SW7 2AZ',
                                                'SW7 2BT',
                                                'SW7 2BU',
                                                'SW7 2DD']])
        >>> locator.get_population_of_postcode([['SW7  2']], True)
        """

        # Delete the space of geography section in population csv
        # and assign them to a new column
        if postcodes == 0:
            return 0
        self.population_data["wts"] = self.population_data["geography"]\
            .apply(lambda x: x.replace(" ", ""))
        ppt_all = []
        # Sector postcode
        if sector:

            for dim_arr in postcodes:
                ppt_sec_list = []

                for sec in dim_arr:

                    sec_wts = sec.replace(" ", "")
                    # Slice out the data which has the same sector code
                    ppt_sec_DF = self.population_data.loc[
                        self.population_data.wts == sec_wts]
                    # Check whether is there a respective sector
                    if ppt_sec_DF.shape[0] == 0:
                        ppt_sec_list.append(0)
                    else:
                        ppt_sec_list.append(
                            ppt_sec_DF.iloc[0]
                            ["Variable: All usual residents; measures: Value"])

            # Append the population to the list of population of the postcodes
            ppt_all.append(ppt_sec_list)

        # Unit postcode
        else:
            self.postcodes_data["sec_wts"] = self.postcodes_data["Postcode"]\
                .apply(lambda x: x[0:-2].replace(" ", ""))

            for dim_arr in postcodes:

                dim_arr = [k[:-2].replace(" ", "") for k in dim_arr]
                ppt_unit_list = len(dim_arr) * [0]
                temp_set_unit = set(dim_arr)

                for key_sec_unit in temp_set_unit:
                    layout_unit_rows = self.postcodes_data.loc[
                        self.postcodes_data.sec_wts == key_sec_unit]
                    n_unit = layout_unit_rows.shape[0]
                    totalppt_sec_DF = self.population_data.loc[
                        self.population_data.wts == key_sec_unit]

                    for i in range(len(dim_arr)):

                        if key_sec_unit == dim_arr[i]:
                            # Check is there a respective unit or sector
                            if n_unit == 0:
                                ppt_unit_list[i] = 0
                            else:

                                if totalppt_sec_DF.shape[0] == 0:
                                    ppt_unit_list[i] = 0
                                else:
                                    totalppt_sec = totalppt_sec_DF.iloc[0][
                                        'Variable: All usual residents;' +
                                        'measures: Value']
                                    ppt_unit_list[i] = totalppt_sec / n_unit

                ppt_all.append(ppt_unit_list)

            self.postcodes_data.drop('sec_wts', axis=1, inplace=True)

        self.population_data.drop('wts', axis=1, inplace=True)
        return ppt_all
