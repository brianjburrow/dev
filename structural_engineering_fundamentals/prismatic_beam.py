import numpy as np

class prismatic_element:
    def __init__(self, A, I, L, E, orientation, nnAN = 0, nnAP = 0):
        # A, I, L are geometric properties of the prismatic_element
        # each is assumed uniform
        self.k_11 = A*E/L
        self.k_22 = 12*E*I/(L**3)
        self.k_23 = -6*E*I/(L**2)
        self.k_33 = 4*E*I/L 
        # E is the material property (uniform) of the beam
        self.klBB = None 
        self.klAA = None 
        self.klAB = None 
        self.klBA = None
        self.rotation_matrix = None

        self.node_number_at_negative_end = nnAN 
        self.node_number_at_positive_end = nnAP

        self.orientation = orientation 
        self.set_up_local_system()


        pass

    def compute_generic_positive_stiffness_matrix(self):
        mat = np.zeros([3, 3])
        mat[0,0] = self.k_11 
        mat[1,1] = self.k_22 
        mat[1,2] = self.k_23 
        mat[2,1] = self.k_23 
        mat[2,2] = self.k_33
        return mat

    def compute_klBB(self):
        mat = self.compute_generic_positive_stiffness_matrix()
        mat[1,2] = -mat[1,2]
        mat[2,1] = -mat[2,1]
        self.klBB = mat 
        pass 

    def compute_klAA(self):
        self.klAA = self.compute_generic_positive_stiffness_matrix()
        pass 

    def compute_klBA(self):
        mat = -self.compute_generic_positive_stiffness_matrix()
        mat[2,1] = -mat[2,1]
        mat[2,2] = -mat[2,2]/2 
        self.klBA = mat 
        pass 

    def compute_klAB(self):
        mat = -self.compute_generic_positive_stiffness_matrix()
        mat[1,2]  = -mat[2,1]
        mat[2,2]   = -mat[2,2]/2 
        self.klBA = mat 
        pass 

    def compute_rotation_matrix(self):
        rotation_matrix = np.zeros([3,3])
        rotation_matrix[0,0] = np.cos(self.orientation)
        rotation_matrix[1,1] = np.cos(self.orientation)
        rotation_matrix[0,1] = np.sin(self.orientation)
        rotation_matrix[1,0] = -np.sin(self.orientation)
        rotation_matrix[2,2] = 1
        self.rotation_matrix = rotation_matrix 
        pass

    def set_up_local_system(self):
        self.compute_klBB()
        self.compute_klAA()
        self.compute_klAB()
        self.compute_klBA()
        self.compute_rotation_matrix()
        pass


