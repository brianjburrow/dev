import numpy as np
class free_vibration:
    def __init__(self, mass, stiffness, damping_coefficient, initial_pos = 0, initial_velocity = 0, damping_type = "VISCOUS"):

        self.stiffness = stiffness
        self.mass = mass 
        self.damping_coefficient = damping_coefficient

        self.initial_pos = initial_pos 
        self.initial_velocity = initial_velocity

        self.damping_type = damping_type.upper()

        damping_options = ["VISCOUS", "UNDAMPED", "COULOMB"] 

        if damping_type not in damping_options:
            print("Supported damping types:", damping_options)
            print("Damping Type %s not supported."%damping_type)
            assert(damping_type in damping_options)

        self._compute_system_properties()

        self._compute_integration_constants()
        
        self._set_position_function_reset()

    def _compute_system_properties(self):
        self.natural_freq = np.sqrt(self.stiffness / self.mass)

        self.damping_ratio = self.damping_coefficient / (2*np.sqrt(self.stiffness * self.mass))

        if self.damping_type != "UNDAMPED":
            self.damped_natural_freq = self.natural_freq*np.sqrt(1 - self.damping_ratio)

        if self.damping_type == "VISCOUS":
            self.damped_period = 2*np.pi / self.damped_natural_freq

        self.logarithmic_decrement = 2*np.pi * self.damping_ratio / np.sqrt(1 - self.damping_ratio**2)
        pass

    def _compute_integration_constants(self):
        if self.damping_type == "VISCOUS":
            if self.damping_ratio == 0:
                self.damping_type = "UNDAMPED"

            if self.damping_ratio < 1:
                # Under Damped
                self.amplitude = np.sqrt(self.initial_pos**2 + ((self.initial_velocity + self.damping_ratio*self.natural_freq*self.initial_pos)/self.damped_natural_freq)**2)
                self.phase_angle = np.arctan(self.damped_natural_freq*self.initial_pos / (self.initial_velocity + self.damping_ratio*self.natural_freq*self.initial_pos))
                self.compute_position_response = self._position_response_viscously_underdamped

            if self.damping_ratio == 1:
                self.const_1 = self.initial_velocity + self.natural_freq*self.initial_pos
                self.compute_position_response = self._position_response_viscously_critically_damped

            if self.damping_ratio > 1:
                self.const_1 = 1./(2.*np.sqrt(self.damping_ratio**2 - 1))
                self.const_2 = self.initial_velocity / self.natural_freq + self.initial_pos*(self.damping_ratio + np.sqrt(self.damping_ratio**2 - 1))
                self.const_3 = -self.initial_velocity / self.natural_freq + self.initial_pos*(-self.damping_ratio + np.sqrt(self.damping_ratio**2 - 1))
                self.exp_const_2 = self.natural_freq*np.sqrt(self.damping_ratio**2 - 1)
                self.exp_const_3 = -self.natural_freq*np.sqrt(self.damping_ratio**2 - 1)
                self.compute_position_response = self._position_response_viscously_overdamped

        if self.damping_type == "UNDAMPED":
            self.amplitude = np.sqrt(self.initial_pos**2 + (self.initial_velocity / self.natural_freq)**2)
            self.phase_angle = np.arctan(self.natural_freq*self.initial_pos / self.initial_velocity)
            self.compute_position_response = self._position_response_undamped
        pass

        
    def _position_response_undamped(self, times):
        position_t = self.amplitude*np.sin(self.natural_freq*times + self.phase_angle)
        return position_t

    def _position_response_viscously_underdamped(self, times):
        position_t = self.amplitude*np.exp(-self.damping_ratio*self.natural_freq*times)
        position_t *= np.sin(self.damped_natural_freq*times + self.phase_angle)
        return position_t

    def _position_response_viscously_overdamped(self, times):
        position_t = self.const_1*np.exp(-self.damping_ratio*self.natural_freq*times)
        position_t *= (self.const_2*np.exp(self.exp_const_2*times) + self.const_3*np.exp(self.exp_const_3*times))
        return position_t

    def _position_response_viscously_critically_damped(self, times):
        position_t = np.exp(-self.natural_freq*times)*(self.initial_pos + self.const_1*times)
        return position_t 

    def _set_position_function_reset(self):
        # we will overwrite self.compute_position_response for speed
        # but if the user changes the damping ratio/type, we may need to
        # reset it to default.  So, we store it in self.reset_function
        # see self.compute_position_response for details
        self.reset_function = self.compute_position_response
        pass
        

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    mass = 1 # kg
    stiffness = 400 # N / m
    damping_coeff = 200
    damping_coeff = 2*np.sqrt(stiffness*mass)
    a = free_vibration(mass, stiffness, damping_coeff, initial_pos = 1, initial_velocity = 0)

    times = np.linspace(0, 10, 1000)

    x_t = a.compute_position_response(times)

    plt.plot(times, x_t)
    plt.show()
