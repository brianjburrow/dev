import numpy as np 

class single_dof_underdamped_harmonic_oscillator():
    def __init__(self, mass, stiffness, damping_coefficient, natural_frequency):
        self. mass = mass 
        self.stiffness = stiffness 
        self.damping_coefficient = damping_coefficient 
        self.natural_frequency = natural_frequency

        self._compute_derived_system_properties()
        pass 

    def _compute_derived_system_properties(self):
        self.critical_damping_coefficient = 2*self.mass*self.natural_frequency 
        
        self.damping_ratio = self.damping_coefficient / self.critical_damping_coefficient

        self.natural_frequency_of_damped_vibration = self.natural_frequency*np.sqrt(1 - self.damping_ratio**2)
        
        self.natural_period = 2*np.pi / self.natural_frequency 

        self.damped_natural_period = 2*np.pi / self.natural_frequency_of_damped_vibration
        pass

    def update_damping_ratio(self, new_damping_ratio):
        # allows for direct updating of the damping ratio, 
        # instead of updating the damping coefficient directly
        self.damping_ratio = new_damping_ratio 
        self.natural_frequency_of_damped_vibration = self.natural_frequency*np.sqrt(1 - self.damping_ratio**2)
        self.damped_natural_period = 2*np.pi / self.natural_frequency_of_damped_vibration
        self.damping_coefficient = self.damping_ratio * self.critical_damping_coefficient
        pass

    def position_func_free_vibration(self, times, init_positions, init_velocities):
        times_2d = times
        init_pos = init_positions
        init_vel = init_velocities
        u_t = np.exp(-self.damping_ratio*self.natural_frequency*times_2d)* \
            (init_pos*np.cos(self.natural_frequency_of_damped_vibration*times_2d) + \
                np.sin(self.natural_frequency_of_damped_vibration*times_2d)*(init_vel + \
                    self.damping_ratio*self.natural_frequency*init_pos ) / self.natural_frequency_of_damped_vibration )
        
        # returns an n_initial_condition x n_times array of positions
        return u_t 

    def position_func(self, forcing_amplitude, forcing_frequency, times, init_positions, init_velocities, forcing_type):
        transient_solution = self.position_func_free_vibration(times, init_positions, init_velocities)

        if forcing_type.upper() == "SINE":
            steady_state_solution = np.sin(forcing_frequency*times)*forcing_amplitude/(self.stiffness*(1 - (forcing_frequency/self.natural_frequency)**2) )
        else:
            print("COSINE Forcing not defined yet")
            osuhe
        return steady_state_solution, transient_solution

    def steady_state_solution(self, forcing_amplitude, forcing_frequency, times, init_position, init_velocities, forcing_type):
        z_m_interior_1 = (2*self.natural_frequency*self.damping_ratio)**
        z_m_interior_2 = (1 / forcing_frequency**2)*(self.natural_frequency**2 - forcing_frequency**2)**2
        z_m = np.sqrt(z_m_interior_1 + z_m_interior_2)

        phi = np.arctan((2*forcing_frequency*self.natural_frequency*self.damping_ratio)/(forcing_frequency**2 - self.natural_frequency**2))
        pass
    def position_func_zero_initials_damped_resonance_sine_forcing(self, forcing_amplitude, times, init_positions, init_velocities):
        u_st_0 = forcing_amplitude / self.stiffness # maximum static deflection
        coeff = u_st_0 / 2. / self.damping_ratio 

        interior_1 = np.cos(self.natural_frequency_of_damped_vibration * times) + (self.damping_ratio/np.sqrt(1 - self.damping_ratio**2))*np.sin(self.natural_frequency_of_damped_vibration*times)
        interior_2 = np.cos(self.natural_frequency*times)

        solution = coeff * (np.exp(-self.damping_ratio*self.natural_frequency*times)*interior_1 - interior_2)
        return solution

    def free_vibration_envelope_curve(self, times, init_positions, init_velocities):
        time_2d = times
        init_pos = init_positions
        init_vel = init_velocities

        envelope_coefficient = np.sqrt(init_pos**2 + ((init_vel + self.damping_ratio*init_pos)/self.natural_frequency_of_damped_vibration)**2)
        positive_envelope_curves = envelope_coefficient*np.exp(-self.damping_ratio*self.natural_frequency*time_2d)
        return positive_envelope_curves

    def demo_free_vibration_damping(self):
        from matplotlib import pyplot as plt
        damping_ratio_set = [0.02, 0.05, 0.1, 0.2]

        times = np.linspace(0, 20, 1000)*self.natural_period

        fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = [10, 10])
        for idx, ax in enumerate(axes.reshape(-1)):
            self.update_damping_ratio(damping_ratio_set[idx])

            print(self.damping_ratio)
            init_pos = np.array([1])
            init_vel = np.array([0])

            u_t = self.position_func_free_vibration(times, init_pos, init_vel)
            ax.plot(times, u_t, label = "$\zeta$ = %f Percent"%(100*self.damping_ratio))
            ax.legend()
            ax.set_xlabel("Time")
            ax.set_ylabel("Position")
        plt.show()
        pass

    def demo_free_vibration_envelope_curves(self):
        from matplotlib import pyplot as plt 
        
        times = np.linspace(0, 20, 1000)*self.natural_period
        
        u_t = self.position_func_free_vibration(times, 0.5, 1.)
        pos_envelope = self.free_vibration_envelope_curve(times, 0.5, 1.)

        self.update_damping_ratio(0)
        u_undamped = self.position_func_free_vibration(times, np.array([1]), np.array([1]))
        plt.figure(figsize = [10,10])
        plt.plot(times, u_t, c = 'black', label = "Damped System")
        plt.plot(times, pos_envelope, c = 'blue', label = None)
        plt.plot(times, -pos_envelope, c = 'blue', label = "Envelope Curve")
        plt.plot(times, u_undamped, c = 'red', label = "Undamped System", linewidth = 1, linestyle = ':')
        plt.xlabel("Time")
        plt.ylabel("Deformation")
        plt.legend()
        plt.show()
        pass

    def demo_damping_ratio_vs_frequency_ratio(self):
        from matplotlib import pyplot as plt
        damping_ratio_set = np.linspace(0, 1, 1000)

        frequency_ratio_set = np.zeros_like(damping_ratio_set)

        for index, damping_ratio in enumerate(damping_ratio_set):
            self.update_damping_ratio(damping_ratio)
            frequency_ratio_set[index] = self.natural_frequency_of_damped_vibration / self.natural_frequency

        plt.figure(figsize = [10,10])
        plt.plot(damping_ratio_set, frequency_ratio_set)
        plt.xlabel("Damping Ratio $\zeta$")
        plt.ylabel("$\omega_D$ / $\omega_n$ = $T_n$ / $T_D$")
        plt.show()
        pass

    def demo_logarithmic_decrement(self):
        from matplotlib import pyplot as plt 
        damping_ratio_set = np.linspace(0.001, 0.9, 1000)

        log_decrement_approx = 2*np.pi*damping_ratio_set
        log_decrement = 2*np.pi*damping_ratio_set/np.sqrt(1 - damping_ratio_set**2)
        plt.figure(figsize = [10,10])
        plt.plot(damping_ratio_set, log_decrement_approx, linestyle = "--", c = 'blue', label = "$\delta$ = 2$\pi$$\zeta$")
        plt.plot(damping_ratio_set, log_decrement, linestyle = '-', c = 'black', label = "$\delta$ = 2$\pi$$\zeta$/$\sqrt{1 - \zeta^2}$")

        plt.xlabel("Damping Ratio $\zeta$")
        plt.ylabel("Logarithmic Decrement $\delta$")
        plt.legend()
        plt.show()
        pass

    def demo_harmonic_vibration_of_undamped_systems(self):
        from matplotlib import pyplot as plt
        self.update_damping_ratio(0)

        times = np.linspace(0, 20, 1000)*self.natural_period # 20/(2*pi) periods
        forcing_frequency = self.natural_frequency / 8. 
        forcing_amplitude = 6

        steady_state, transient_state = self.position_func(forcing_amplitude, forcing_frequency, times, 1., 1., "SINE")

        fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = [10,10])

        axes[0].plot(times, forcing_amplitude*np.sin(forcing_frequency*times))
        axes[0].set_xlabel("Time")
        axes[0].set_ylabel("Applied Force")

        axes[1].plot(times, steady_state, label = "Steady State")
        axes[1].plot(times, steady_state + transient_state, label = "Total Response")
        axes[1].legend()
        axes[1].set_xlabel("Time")
        axes[1].set_ylabel("Position")

        plt.show()
        pass

    def demo_harmonic_vibration_of_damped_systems(self):
        from matplotlib import pyplot as plt
        self.update_damping_ratio(0.2)

        times = np.linspace(0, 20, 1000)*self.natural_period # 20/(2*pi) periods
        forcing_frequency = self.natural_frequency / 8. 
        forcing_amplitude = 6

        steady_state, transient_state = self.position_func(forcing_amplitude, forcing_frequency, times, 1., 1., "SINE")

        fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = [10,10])

        axes[0].plot(times, forcing_amplitude*np.sin(forcing_frequency*times))
        axes[0].set_xlabel("Time")
        axes[0].set_ylabel("Applied Force")

        axes[1].plot(times, steady_state, label = "Steady State")
        axes[1].plot(times, steady_state + transient_state, label = "Total Response")
        axes[1].legend()
        axes[1].set_xlabel("Time")
        axes[1].set_ylabel("Position")

        plt.show()
        pass

    def demo_resonant_frequency(self):
        from matplotlib import pyplot as plt 
        initial_velocity = 0 
        initial_position = 0

        forcing_frequency = self.natural_frequency
        forcing_amplitude = 1.0 
        times = np.linspace(0, 20*self.natural_period, 1000)
        u_steady, u_transient = self.position_func(forcing_amplitude, forcing_frequency, times, initial_position, initial_velocity, "SINE")
        plt.plot(times, u_steady + u_transient)
        plt.show()
        pass


        







if __name__ == "__main__":
    mass = 1
    stiffness = 10
    damping_coefficient = 0.02
    natural_frequency = 1
    oscillator = single_dof_underdamped_harmonic_oscillator(mass,
    stiffness,
    damping_coefficient,
    natural_frequency)

    #oscillator.demo_free_vibration_damping()

    oscillator = single_dof_underdamped_harmonic_oscillator(mass,
    stiffness,
    damping_coefficient,
    natural_frequency)
    #oscillator.demo_free_vibration_envelope_curves()
    #oscillator.demo_damping_ratio_vs_frequency_ratio()
    #oscillator.demo_logarithmic_decrement()
    #oscillator.demo_harmonic_vibration_of_undamped_systems()
    #oscillator.demo_harmonic_vibration_of_damped_systems()
    oscillator.demo_resonant_frequency()
    pass

    
