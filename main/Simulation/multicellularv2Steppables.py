
import numpy as np

from cc3d.core.PySteppables import *


def distance(coord_1, coord_2, p=None):

    if p:
        print(coord_1 - coord_2)
    x = np.linalg.norm(coord_1 - coord_2)
    return x
    
    

class ConstraintInitializerSteppable(SteppableBasePy):
            
    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)

    def start(self):
        # any code in the start function runs before MCS=0
        # Set Parameter values here and create persistent variables
        self.field = CompuCell.getConcentrationField(self.simulator, "glucose")
        self.x_lattice_size=500
        self.y_lattice_size=500
        # lambda chemo
        self.lambda_chemo=1.0 
        self.lambda_persistence = 3.0
        self.tau_p = 50
        self.tau_s = 5000 ## I should include this in a dictionary so that I can pass it to the mitosis steppable
        
        
        ################################################################################################
        ################################################################################################
        
        ## CM: I am going to try to initialize a population of 200 cells here:
        
        size = 5.0 # Size of the cells
        num_of_cells = 200 # How many cells to initialize
        
        # The initialization border gives a breadth of 10 units from any of the borders
        x_min = 10 
        x_max = 490
        y_min = 10
        y_max = 490
        
        
        for i in range (num_of_cells): 

                x = np.random.randint(x_min,x_max)
                y = np.random.randint(y_min,y_max)
                z = 0
                
                self.cell_field[x:x + size - 1, y:y + size - 1, 0] = self.new_cell(self.CELL1)
                
                
        
        ##############################################################################################
        ##############################################################################################


        for cell in self.cellList:
            
           
            if cell.type == 1:
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "glucose") # Cell and field it responds to
   
                cd.setLambda(self.lambda_chemo) # Strength of attraction
                    # cd.assignChemotactTowardsVectorTypes([self.MEDIUM, self.BACTERIUM])
                cell.lambdaVolume= 5.0 # In this simulation I am scanning lambda volume between min_value and max_value
                cell.targetVolume= 50.0
                angle = np.random.uniform(0,np.pi*2)
                # Make sure ExternalPotential plugin is loaded
                cell.lambdaVecX = -(self.lambda_persistence*np.cos(angle))  # force component pointing along X axis - towards positive X's
                cell.lambdaVecY = -(self.lambda_persistence*np.sin(angle))  # force component pointing along Y axis - towards negative Y's       
                
                # we will give the cell a way to keep track of its velocity by tracking the last time and x position we changed direction
                # cell.dict['lastDirectionChangeTime'] = 100
                cell.dict['lastDirectionChangeTime'] = 0
                cell.dict['lastDirectionChangePosition'] = [cell.xCOM, cell.yCOM]
                

                

        #patchy gradient
        secrConst = 1 #starting constant
        
        self.dict = {'starting_location': None}

        self.initialize_season()



    def step(self,mcs):
                # # any code in the start function runs before MCS=0
        # # Set Parameter values here and create persistent variables
        # arguments are (name of the data series, x, y)
        # iterating over all cells in simulation

        ## Checking to see if it is time to update the season.
        if not (mcs%self.tau_s) and (mcs > 0):
            self.initialize_season()

        for cell in self.cellList:
            
            if cell.type == 1:

            # This line below is to fix the fact that sometimes the cell's center of mass does not move far
            # enough away after the lambda sign switches, and therefore the cell flips multiple times in a row.
                if mcs - cell.dict['lastDirectionChangeTime'] >= self.tau_p:
                    diffx = cell.xCOM - cell.dict['lastDirectionChangePosition'][0]
                    diffy = cell.yCOM - cell.dict['lastDirectionChangePosition'][1]
                    diffs = [diffx, diffy]
                    norm = np.linalg.norm(diffs)
                    if not (mcs%10):
                        cell.lambdaVecX = -(self.lambda_persistence*(diffx/norm))
                        cell.lambdaVecY = -(self.lambda_persistence*(diffy/norm))
                    cell.dict['lastDirectionChangeTime'] = mcs
                    cell.dict['lastDirectionChangePosition'] = [cell.xCOM, cell.yCOM]


            
    def finish(self):
        """
        Finish Function is called after the last MCS
        """

    def on_stop(self):
        # this gets called each time user stops simulation
        return
        
    def initialize_season(self, baseline = .50, baseline_add = 0):

            #if first starting locations, initialize empty list and pick random location
            secrConst = baseline+baseline_add
            if self.dict['starting_location'] is None:
                self.dict['starting_location'] = np.zeros(shape=(2)) # Initialize empty coordinates
                self.dict['starting_location'][0] = np.random.randint(500) # Setting the x-coordinate
                self.dict['starting_location'][1] = np.random.randint(500) # Setting the y-coordinate 
                
            else:
                self.dict['starting_location'][0] = np.random.randint(500) # Setting the x-coordinate
                self.dict['starting_location'][1] = np.random.randint(500) # Setting the y-coordinate

            bL_corner = np.array([0,0])
            uR_corner = np.array([500,500])
            uL_corner = np.array([0,500])
            bR_corner = np.array([500, 0])
            distances = [distance(self.dict['starting_location'], bL_corner), distance(self.dict['starting_location'], uR_corner),
                         distance(self.dict['starting_location'], uL_corner), distance(self.dict['starting_location'], bR_corner)]
            max_distance = max(distances)

            #iterate through all pixels, pixel is either secrConstant*distance or zero
            for x, y, z in self.everyPixel(1, 1, 1):
                tracker = self.cellField[x,y,z] 
                d = distance(np.array([x,y]), self.dict['starting_location'])
                dist = abs(max_distance-d)/max_distance
                outcome = np.random.binomial(1, dist)
                if tracker:
                    if outcome:
                        self.field[x,y,z] = secrConst*dist+baseline
                    else:
                        self.field[x,y,z] = 0
                else:
                    if outcome:
                        self.field[x,y,z] = secrConst*dist+baseline
                    else:
                        self.field[x,y,z] = 0
            

        
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)
    
    
    def step(self, mcs):
        
        self.tau_s = 5000 ## Needs to match with the variable in the main steppable class

        if not (mcs%self.tau_s) and (mcs > 0):
            cells_to_divide=[]
            for cell in self.cellList:
                    cells_to_divide.append(cell)
            for cell in cells_to_divide:
                # Other valid options
                # self.divide_cell_orientation_vector_based(cell,1,1,0)
                self.divide_cell_along_major_axis(cell)
                # self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        # Target volume should stay the same.
        self.parent_cell.targetVolume /= 1.0

        self.clone_parent_2_child()
     
  
class AdhesionMoleculesSteppables(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
                

    def step(self, mcs):
        
        self.tau_s = 5000 ## Needs to match with the variable in the main steppable class
        
        if mcs == 0:   
            for cell in self.cell_list:
                adhesion_molecule_vector = self.adhesionFlexPlugin.getAdhesionMoleculeDensityVector(cell)
                # for i, elt in enumerate(adhesion_molecule_vector):
                    # adhesion_molecule_vector[i] = np.random.uniform()
                new_molecules = np.random.uniform(0,1,len(adhesion_molecule_vector))
                self.adhesionFlexPlugin.setAdhesionMoleculeDensityVector(cell, new_molecules)
                
        #if not(mcs%self.tau_s) and (mcs>0):
            

        
class DeathSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
        
        self.tau_s = 5000 ## Needs to match with the variable in the main steppable class

        if not (mcs%self.tau_s) and (mcs > 0):
            if len(self.cellList) > 200:
                num_to_die = len(self.cellList) - 200
                while num_to_die > 0:
                    for cell in self.cellList:
                        if (cell.type == 1) and (cell.targetVolume > 0):
                            coin = np.random.uniform(0,1)
                            if coin > 0.5:
                                cell.targetVolume = 0
                                cell.lambdaVolume = 100000
                                num_to_die -= 1
  