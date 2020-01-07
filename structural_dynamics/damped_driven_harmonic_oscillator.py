import numpy as np 
import matplotlib.pyplot as plt
class damped_driven_harmonic_oscillator:
    def __init__(self, mass, damping_coefficient, stiffness):
        self.mass = mass 
        self.stiffness = stiffness 
        self.damping_coefficient = damping_coefficient
        print(damping_coefficient)
        print(self.mass)
        print(self.stiffness)
        print(self.damping_coefficient)
        self._compute_derived_system_properties()

        pass

    def _compute_derived_system_properties(self):
        self.natural_frequency = np.sqrt(self.stiffness / self.mass)

        self.natural_frequency_of_damped_vibration = np.sqrt(self.stiffness / self.mass + self.damping_coefficient**2 / (4*self.mass**2))

        self.quality_factor = np.sqrt(self.mass*self.stiffness / self.damping_coefficient)

        self.natural_period = 2*np.pi / self.natural_frequency 

        self.damped_natural_period = 2*np.pi / self.natural_frequency_of_damped_vibration
        pass

    def _compute_system_coupled_force_properties(self, forcing_frequency, forcing_amplitude):
        self.forcing_frequency = forcing_frequency
        self.forcing_amplitude = forcing_amplitude
        value = self.stiffness - self.mass*forcing_frequency**2
        self.phase_angle = np.arctan(forcing_frequency*forcing_amplitude / value)

        self.amplitude = self.forcing_amplitude / np.sqrt( value + forcing_frequency**2 * self.damping_coefficient**2  )
        pass

    def free_overdamped_solution(self, initial_position, initial_velocity, times):
        C1 = initial_position 
        C2 = initial_velocity / self.natural_frequency
        
        radicand = self.damping_coefficient**2  - 4 * self.stiffness * self.mass 

        x_1 = C1*np.exp(times*( - self.damping_coefficient + np.sqrt(radicand)) / 2. / self.mass)

        x_2 = C2*np.exp(times*( - self.damping_coefficient - np.sqrt(radicand)) / 2. / self.mass)
        return x_1 + x_2 

    def free_critically_damped_solution(self, initial_position, initial_velocity, times):
        C1 = initial_position 
        C2 = initial_velocity 

        coeffs = np.eq(-self.damping_coefficient * times / 2. / self.mass)

        return C1 * coeffs + C2*times*coeffs 

    def free_underdamped_solution(self, initial_position, initial_velocity, times):
        C1 = initial_position 
        C2 = initial_velocity / self.natural_frequency

        coeffs = np.exp(-self.damping_coefficient * times / 2. / self.mass)

        x_1 = C1 * coeffs * np.cos(self.natural_frequency*times)

        x_2 = C2 * coeffs * np.sin(self.natural_frequency * times)

        return x_1 + x_2

    def cosine_driven_solution(self, driving_frequency, driving_amplitude, initial_position, initial_velocity, times):

        interior_1 = (self.stiffness - self.mass*driving_frequency**2)**2
        interior_2 = driving_frequency**2 * self.damping_coefficient**2
        rho = np.sqrt(interior_1 + interior_2)

        phase_angle = np.arctan(driving_frequency * self.damping_coefficient / (self.stiffness - self.mass*driving_frequency**2))

        steady_state_response = driving_amplitude * np.cos(driving_frequency*times - phase_angle)/rho 
        
        if self.quality_factor < 0.5:
            transient_response = self.free_overdamped_solution(initial_position, initial_velocity, times)
        elif self.quality_factor == 0.5:
            transient_response = self.free_critically_damped_solution(initial_position, initial_velocity, times)
        elif self.quality_factor > 0.5:
            transient_response = self.free_underdamped_solution(initial_position, initial_velocity, times)

        return steady_state_response + transient_response

if __name__ == "__main__":
    mass = 4. # KG

    damping_coefficient = 1. # N s / m

    stiffness = 1. # N / m

    force_amplitude = 1 # N

    force_frequency = .75 # rad / s

    initial_position = 10 # m

    initial_velocity = 0.01 # m / s

    oscillator = damped_driven_harmonic_oscillator(mass, damping_coefficient, stiffness)

    times = np.linspace(0, 5 * oscillator.damped_natural_period, 1000)

    u_t = oscillator.cosine_driven_solution(force_frequency, force_amplitude, initial_position, initial_velocity, times)

    plt.plot(times/oscillator.damped_natural_period, u_t)
    plt.show()
