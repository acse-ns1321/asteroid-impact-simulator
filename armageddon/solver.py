import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class Planet():
    """
    The class called Planet is initialised with constants appropriate
    for the given target planet, including the atmospheric density profile
    and other constants
    """

    def __init__(self, atmos_func='exponential',
                 atmos_filename='./armageddon/resources/' +
                 'AltitudeDensityTable.csv',
                 Cd=1., Ch=0.1, Q=1e7, Cl=1e-3, alpha=0.3, Rp=6371e3,
                 g=9.81, H=8000., rho0=1.2):
        """
        Set up the initial parameters and constants for the target planet
        Parameters
        ----------
        atmos_func : string, optional
            Function which computes atmospheric density, rho, at altitude, z.
            Default is the exponential function rho = rho0 exp(-z/H).
            Options are 'exponential', 'tabular' and 'constant'
        atmos_filename : string, optional
            Name of the filename to use with the tabular atmos_func option
        Cd : float, optional
            The drag coefficient
        Ch : float, optional
            The heat transfer coefficient
        Q : float, optional
            The heat of ablation (J/kg)
        Cl : float, optional
            Lift coefficient
        alpha : float, optional
            Dispersion coefficient
        Rp : float, optional
            Planet radius (m)
        rho0 : float, optional
            Air density at zero altitude (kg/m^3)
        g : float, optional
            Surface gravity (m/s^2)
        H : float, optional
            Atmospheric scale height (m)
        """

        # Input constants
        self.Cd = Cd
        self.Ch = Ch
        self.Q = Q
        self.Cl = Cl
        self.alpha = alpha
        self.Rp = Rp
        self.g = g
        self.H = H
        self.rho0 = rho0
        self.atmos_filename = atmos_filename
        self.tabular_dict = {}

        try:
            # set function to define atmoshperic density
            if atmos_func == 'exponential':
                # rhoa will change as z change
                self.rhoa = lambda z: self.rho0 * np.exp(-z / self.H)
            elif atmos_func == 'tabular':
                # between z_{i} and z_{i+1}, rho_a = rho_i * exp((z_i-z)/H_i)

                Density_Table = pd.read_csv(self.atmos_filename, sep=' ',
                                            skiprows=6,
                                            names=['Altitude',
                                                   'Atmospheric_density',
                                                   'H'])
                self.z_zero_rho = float(
                    Density_Table.iloc[0].Atmospheric_density)
                self.z_max = float(Density_Table.iloc[-1].Altitude)
                self.z_max_rho = float(
                    Density_Table.iloc[-1].Atmospheric_density)
                for i in range(len(Density_Table)):
                    self.tabular_dict[int(Density_Table.iloc[i].Altitude)] = (
                        float(Density_Table.iloc[i].Atmospheric_density),
                        float(Density_Table.iloc[i].H))
                self.rhoa = lambda x: self.z_zero_rho if int(x) <= 0 \
                    else self.tabular_dict[int(int(x/10)*10)][0] * np.exp(
                    (int(int(x/10)*10)-x) /
                    self.tabular_dict[int(int(x/10)*10)][1])\
                    if int(x) < self.z_max else self.z_max_rho

            elif atmos_func == 'constant':
                self.rhoa = lambda x: rho0
            else:
                raise NotImplementedError(
                    "atmos_func must be 'exponential', 'tabular' or\
                    'constant'")
        except NotImplementedError:
            print("atmos_func {} not implemented yet.".format(atmos_func))
            print("Falling back to constant density atmosphere for now")
            self.rhoa = lambda x: rho0

    def solve_atmospheric_entry(
            self, radius, velocity, density, strength, angle,
            init_altitude=100e3, dt=0.005, radians=False):
        """
        Solve the system of differential equations for a given impact scenario
        Parameters
        ----------
        radius : float
            The radius of the asteroid in meters
        velocity : float
            The entery speed of the asteroid in meters/second
        density : float
            The density of the asteroid in kg/m^3, which is a constant
        strength : float
            The strength of the asteroid (i.e. the maximum pressure it can
            take before fragmenting) in N/m^2
        angle : float
            The initial trajectory angle of the asteroid to the horizontal
            By default, input is in degrees. If 'radians' is set to True, the
            input should be in radians
        init_altitude : float, optional
            Initial altitude in m
        dt : float, optional
            The output time step, in s
        radians : logical, optional
            Whether angles should be given in degrees or radians. Default=False
            Angles returned in the dataframe will have the same units as the
            input
        Returns
        -------
        Result : DataFrame
            A pandas dataframe containing the solution to the system.
            Includes the following columns:
            'velocity', 'mass', 'angle', 'altitude',
            'distance', 'radius', 'time'
        """

        # Enter your code here to solve the differential equations
        # regulate the angle. the form of the degree is radians
        # if radians == False:
        #     angle = angle * np.pi / 180.
        # else:
        #     pass

        initial_data_asteroid = np.array(
            [velocity, density * 4 / 3 * np.pi * radius ** 3,
             angle * np.pi/180,
             init_altitude, 0, radius])
        data_asteroid, time_list = self.RK4(
            self.formula, initial_data_asteroid, dt, strength, density)

        # reset the form of degree.
        if not radians:
            data_asteroid[:, 2] = data_asteroid[:, 2] * 180 / np.pi
        else:
            pass

        # regulate the output
        # print(data_asteroid)
        velocity = data_asteroid[:, 0]
        mass = data_asteroid[:, 1]
        angle = data_asteroid[:, 2]
        altitude = data_asteroid[:, 3]
        distance = data_asteroid[:, 4]
        radius = data_asteroid[:, 5]

        return pd.DataFrame({'velocity': velocity,
                             'mass': mass,
                             'angle': angle,
                             'altitude': altitude,
                             'distance': distance,
                             'radius': radius,
                             'time': time_list},
                            index=range(data_asteroid.shape[0]))

    def formula(self, t, stepdata_asteroid, strength, density):
        """
        calculate the ode result directly by given parameters
           f[0]:dv/dt
           f[1]:dm/dt
           f[2]:d(theta)/dt
           f[3]:dz/dt
           f[4]:dx/dt
           f[5]:dr/dt
        parameters
        ----------
        stepdata_asteroid: list of parameters
            stepdata_asteroid = [velocity, mass, angle, altitude,
            distance, radius]
        density : float
            The density of the asteroid in kg/m^3, which is a constant
        strength : float
            The strength of the asteroid (i.e. the maximum pressure it can
            take before fragmenting) in N/m^2
        Returns
        -------
        Result :  array
            an array consisting of 6 numbers, as shown below.
            f[0, 1, 2, 3, 4, 5]
            f[0]:dv/dt
            f[1]:dm/dt
            f[2]:d(theta)/dt
            f[3]:dz/dt
            f[4]:dx/dt
            f[5]:dr/dt
        """
        f = np.zeros_like(stepdata_asteroid)
        rhoa = self.rhoa(stepdata_asteroid[3])
        A = np.pi*stepdata_asteroid[5]**2
        u = stepdata_asteroid.copy()
        f[0] = -(self.Cd*rhoa*A*(u[0])**2)/(2*(u[1])) + self.g*np.sin(u[2])
        f[1] = -(self.Ch*rhoa*A*(u[0])**3)/(2*self.Q)
        f[2] = self.g*np.cos(u[2])/(u[0]) - self.Cl*rhoa*A * \
            (u[0])/(2*(u[1])) - u[0]*np.cos(u[2])/(self.Rp+u[3])
        f[3] = -u[0]*np.sin(u[2])
        f[4] = u[0]*np.cos(u[2])/(1 + u[3]/self.Rp)
        f[5] = 0
        if rhoa*u[0]**2 >= strength:
            f[5] = (7/2*self.alpha*rhoa/density)**(1/2)*(u[0])
        return f

    def RK4(self, f, u0, dt, strength, density):
        """
        Solve the a coupled set of (total six) ordinary differential equations.
        Implement RK4 algorithm
        Parameters
        ----------
        f : list
            ODE function lists to solve
        u0 : list
            a list of initial parameters
            u0 = [velocity, mass, angle, altitude, distance, radius]
        dt : float
          The timestep
        strength : float
            The strength of the asteroid (i.e., the ram pressure above which
            fragmentation and spreading occurs) in N/m^2 (Pa)
        density : float
            The density of the asteroid in kg/m^3
        Returns
        -------
        Result : array
            u = [velocity, mass, angle, altitude, distance, radius]
            t = [time]
        """

        u = np.array(u0)
        t = np.array(0)
        u_all = [u0]
        t_all = [0]
        iteration = 0
        while iteration < 1e4:
            k1 = dt*f(t, u, strength, density)
            k2 = dt*f(t + 0.5*dt, u + 0.5*k1, strength, density)
            k3 = dt*f(t + 0.5*dt, u + 0.5*k2, strength, density)
            k4 = dt*f(t + dt, u + k3, strength, density)
            u = u + (1./6.)*(k1 + 2*k2 + 2*k3 + k4)
            if u[0] < 20 or u[3] < -1 or u[1] < 0:
                break
            u_all.append(u)
            t = t + dt
            iteration += 1
            t_all.append(t)
        return np.array(u_all), np.array(t_all)

    def calculate_energy(self, result):
        """
        Function to calculate the kinetic energy lost per unit altitude in
        kilotons TNT per km, for a given solution.
        Parameters
        ----------
        result : DataFrame
            A pandas dataframe with columns for the velocity, mass, angle,
            altitude, horizontal distance and radius as a function of time
        Returns : DataFrame
            Returns the dataframe with additional column ``dedz`` which is the
            kinetic energy lost per unit altitude
        """
        # the result DataFrame

        result = result.copy()

        # kinetic energy = 0.5*m*v^2
        result.insert(len(result.columns), 'dedz', np.array(np.nan))
        result.insert(len(result.columns), 'resi_en', np.array(np.nan))
        # central differential method
        result.resi_en = 1/2 * result.mass * result.velocity**2
        result.dedz = result.resi_en.diff(2)/result.altitude.diff(2)
        result.dedz = result.dedz.shift(-1)
        # convert to TNT
        result.dedz = result.dedz/(4.184e9)
        result.dedz[0] = result.dedz[1]
        del result['resi_en']
        return result.drop(len(result)-1)

        # Replace these lines with your code to add the dedz column to
        # the result DataFrame
        # result = result.copy()
        # result.insert(len(result.columns),
        #               'dedz', np.array(np.nan))
        #
        # return result

    def analyse_outcome(self, result):
        """
        Inspect a pre-found solution to calculate the impact and airburst stats
        Parameters
        ----------
        result : DataFrame
            pandas dataframe with velocity, mass, angle, altitude, horizontal
            distance, radius and dedz as a function of time
        Returns
        -------
        outcome : Dict
            dictionary with details of the impact event, which should contain
            the key ``outcome`` (which should contain one of the following :
            ``Airburst`` or ``Cratering``), as well as the following keys:
            ``burst_peak_dedz``, ``burst_altitude``, ``burst_distance``,
            ``burst_energy``
        """
        outcome = {}

        burst_peak_dedz = result.dedz.max()
        # find index
        dedz_max_index = result.index[(
            result['dedz'] == burst_peak_dedz)].tolist()
        current_index = dedz_max_index[0]
        burst_altitude = result.altitude[current_index]
        burst_distance = result.distance[current_index]
        initial_energy = 1 / 2 * \
            result['mass'][0] * result['velocity'][0] ** 2 / 4.184e12
        current_energy = 1 / 2 * \
            result['mass'][current_index] * \
            result['velocity'][current_index] ** 2 / 4.184e12
        burst_energy = np.abs(initial_energy - current_energy)
        if burst_altitude >= 0:
            outcome = 'Airburst'
        else:
            current = result[result.altitude < 0.].iloc[0]
            burst_peak_dedz = current["dedz"]
            burst_distance = current["distance"]
            residual_energy = 0.5 * \
                float(current['mass'] * (current['velocity'])**2)
            burst_energy = np.abs(initial_energy - residual_energy)
            burst_energy = max(residual_energy, burst_energy)
            outcome = 'Cratering'

        outcome = {'outcome': outcome,
                   'burst_peak_dedz': burst_peak_dedz,
                   'burst_altitude': burst_altitude,
                   'burst_distance': burst_distance,
                   'burst_energy': burst_energy}
        return outcome

    def chelyabinsk(self):
        """
        Determine asteroid parameters strength and radius that
        fits an observed energy deposition curve and plot a fig
        Parameters
        ----------
        Returns
        -------
        """
        # find best fitted parameter
        data = pd.read_csv(
            './armageddon/resources/ChelyabinskEnergyAltitude.csv',
            skiprows=1, names=[
                'altitude', 'dedz'])
        data_peak_dedz = data['dedz'].max()
        max_index = data['dedz'].idxmax()
        data_burst_altitude = data.loc[max_index, 'altitude'] * 1e3  # unit m
        print('exact burst altitude and peak dedz:')
        print([data_burst_altitude, data_peak_dedz])

        # guess a initial strength and raduis
        strength = 1e6  # N/m^2
        radius = 20
        # specified parameters
        velocity = 1.92e4  # m/s
        density = 3300  # kg/m 3
        angle = 18.3  # degrees

        raduis_all = []
        strength_all = []
        altitude_all = []
        dedz_all = []
        err1_all = []
        err2_all = []
        err_all = []
        raduis_all.append(radius)
        strength_all.append(strength)
        error2 = 50
        while error2 > 10:
            # reduce error2
            result = self.solve_atmospheric_entry(
                radius, velocity, density, strength, angle,
                init_altitude=43e3, dt=0.05, radians=False)
            result = self.calculate_energy(result)
            X = result['altitude'] / 1000  # km
            Y = result['dedz']
            altitude_all.append(X)
            dedz_all.append(Y)
            outcome = self.analyse_outcome(result)
            error1 = np.abs(outcome['burst_altitude'] -
                            data_burst_altitude) / 1000
            error2 = np.abs(outcome['burst_peak_dedz'] - data_peak_dedz)
            error = error1 + error2
            err1_all.append(error1)
            err2_all.append(error2)
            err_all.append(error)
            print('error1:')
            print(error1)
            print('error2:')
            print(error2)
            print('present fitted burst_altitude and burst peak dedz:')
            print([outcome['burst_altitude'], outcome['burst_peak_dedz']])

            # narrow error2 also will change error1
            if outcome['burst_peak_dedz'] > data_peak_dedz:
                if error2 > 3:
                    # when error2 is large, step of radius is large too
                    radius -= 0.5
                else:
                    radius -= 0.15
                    strength -= 200
            else:
                if error2 > 20:
                    # when error2 is large, step of radius is large too
                    radius += 0.5
                else:
                    radius += 0.15
                    strength += 200
            raduis_all.append(radius)
            strength_all.append(strength)
        print('min error:')
        min_error = np.min(err_all)
        print(min_error)
        min_error_index = err_all.index(min_error)

        final_radius = raduis_all[min_error_index]
        final_strength = strength_all[min_error_index]
        print('radius')  # final best fitted parameter
        print(final_radius)
        print('strength')  # final best fitted parameter
        print(final_strength)
        final_altitude = altitude_all[min_error_index]
        final_dedz = dedz_all[min_error_index]
        # plot exact energy deposition curve vis fitted curve

        plt.plot(final_altitude, final_dedz, '.',
                 label='fitted energy deposition curve')
        plt.plot(data['altitude'], data['dedz'], '*',
                 label='exact energy deposition curve')
        plt.xlabel('altitude(m)')
        plt.ylabel('Energy Per Unit Length (kt Km^-1)')
        plt.legend()
        plt.show()

    def solver_analytic(self, v0, z0, rho0, theta0, r0, mass, dt, H=8000):
        """
        Solve the system of analytical equations for a given impact scenario
        Parameters
        ----------
        r0 : float
            The radius of the asteroid in meters
        mass : float
            Mass of the asteroid in kg
        velocity : float
            The entery speed of the asteroid in meters/second
        rho0 : float
            The density of the asteroid in kg/m^3, which is a constant
        strength : float
            The strength of the asteroid (i.e. the maximum pressure it can
            take before fragmenting) in N/m^2
        theta0 : float
            The initial trajectory angle of the asteroid to the horizontal
            By default, input is in degrees. If 'radians' is set to True, the
            input should be in radians
        z0 : float, optional
            Initial altitude in m
        dt : float, optional
            The output time step, in s
        H : float
            Atmospheric scale height in m
        Returns
        -------
        Result : DataFrame
            A pandas dataframe containing the solution to the system.
            Includes the following column:
            'velocity', 'altitude'
            List
            A list which contains analytical solution of velocity
        """
        def RHS_analytic(t, state):
            f = np.zeros_like(state)

            velocity = state[0]
            # mass = state[1]
            angle = state[2]
            altitude = state[3]
            # distance = state[4]
            radius = state[5]

            f[0] = (-1.2 * np.exp(-altitude / H) * np.pi * radius**2 *
                    velocity**2) / (2 * mass)
            f[1] = 0
            f[2] = 0
            f[3] = -velocity * np.sin(angle)
            f[4] = velocity * np.cos(angle)
            f[5] = 0
            return f

        def RK4(f, state0, dt, t0=0.):
            # initialise array to store current state and time
            state = np.array(state0)
            t = np.array(t0)

            # initialise list to store state at each time step
            state_all = [state0]

            # initialise list to store all time steps
            t_all = [t0]

            # set limit to number ot iterations
            num_iteration = 0

            while state[3] > 0 and num_iteration < 1e4:
                k1 = dt * f(t, state)
                k2 = dt * f(t + 0.5 * dt, state + 0.5 * k1)
                k3 = dt * f(t + 0.5 * dt, state + 0.5 * k2)
                k4 = dt * f(t + dt, state + k3)
                state = state + (1. / 6.) * (k1 + 2 * k2 + 2 * k3 + k4)
                state_all.append(state)
                t = t + dt
                t_all.append(t)
                num_iteration += 1

            return np.array(state_all[:-1][:]), np.array(t_all[:-1][:])

        state0 = [v0, rho0 * 4. / 3. * np.pi * r0**3, theta0, z0, 0., r0]

        # use RK4 to solve system
        state_all, time_all = RK4(RHS_analytic, state0, dt)

        # obtain variables over time steps
        velocity = state_all[:, 0]
        # mass = state_all[:, 1]
        angle = state_all[:, 2]
        angle = angle / (2 * np.pi) * 360
        altitude = state_all[:, 3]
        # distance = state_all[:, 4]
        # radius = state_all[:, 5]
        a = (-rho0 * np.pi * r0**2 * H) / (np.sin(theta0) * 2 * mass)
        velocity_analytic = np.exp(
            a * np.exp(-altitude / H)) * (v0 / (np.exp(a * np.exp(-z0 / H))))
        return pd.DataFrame({
            'velocity': velocity,
            'altitude': altitude
        }), velocity_analytic
