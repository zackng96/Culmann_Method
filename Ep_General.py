import numpy as np
from utilities import *
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Introduce parameter classes
params = Display_Parameters()
s_params = Sand_Parameters()
f_params = Failure_Parameters()


###############################################################################
# Geometry and Display Units
###############################################################################
# Define default ERSS position and parameters (mm)
SoilLen1 = 5000  # length of which 1st slope extends
SoilLen2 = 2000  # length of which 2nd slope extends
SoilLen_list = [SoilLen1, SoilLen2,
                1000,
                5000]
SoilHt1 = 3500  # elevation of 1st slope peak above wall top
SoilHt2 = 2500  # elevation of 2nd slope peak above 1st peak
SoilHt_list = [SoilHt1, SoilHt2,
               -2000,
               3000]

ERSSHt = 6000 # wall height
alpha = 90 # deg
ERSSTopRight = 3000 - ERSSHt * np.tan((alpha - 90) / 180 * np.pi)
ERSSTopLeft = 1000
ERSSBotRight = 3000
ERSSBotLeft = 0
Wall_Geometry = [ERSSBotLeft, ERSSBotRight, ERSSTopRight, ERSSTopLeft]

# Units
params.unit = 'm'
params.force_unit = 'kN/m'

# Storage of geometric properties in dict
coord_dict = {}
coord_dict['SoilLen_list'] = SoilLen_list
coord_dict['SoilHt_list'] = SoilHt_list
coord_dict['Wall_Geometry'] = Wall_Geometry
coord_dict['Wall_Ht'] = ERSSHt
coord_dict['Wall_alpha'] = alpha

###############################################################################
# Soil Parameters
###############################################################################
s_params.phi = 33 # Culmann analysis assumes cohesionless
s_params.delta = 22
s_params.UnitWeight = 18  # soil unit weight in kN/m3

###############################################################################
# Active Earth Pressure
###############################################################################
# Boundary Control
f_params.xBound = 20000
f_params.yBound = 18000

# Adopt suitable scale by trial and error and convert WeightList to 
# equivalent length

# Trial failure slip control parameters
f_params.FendX = f_params.xBound
f_params.FstartX = 1000
f_params.FstepX = 2000
f_params.Weight2Len = 10

# Introducing coordinates
Wall_Dim(coord_dict, f_params)
Soil_Dim(coord_dict, f_params)
xE, yE = Wall(coord_dict)
xS, yS = Soil(coord_dict)

# Plot main features (wall + soil)
fig, ax = Plotter(xE, yE, xS, yS, params, f_params)

# Plot failure planes and compile failure surfaces for shoelace method later
Failure_Planes(coord_dict, f_params, params, show_FP = True, ax = ax)

# Determine area of each polygon 
WeightFn(coord_dict, s_params, f_params)

Calc_Pa(coord_dict, f_params, s_params, params, 
        show_Wline = True, show_Pline = True, 
        show_Wpoint = False, show_Rpoint = True, show_Rline = True, show_max_Pa = True,
        ax = ax)

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())

ax.set_title('Culmann Graphical Approach (Active Ratio: {}mm/kN)'.format(f_params.Weight2Len))
plt.show()

###############################################################################
# Passive Earth Pressure
###############################################################################

# Boundary Control
f_params.xBound = 50000
f_params.yBound = 20000

# Adopt suitable scale by trial and error and convert WeightList to 
# equivalent length
# Trial failure slip control parameters
f_params.FendX = f_params.xBound
f_params.FstartX = 20000
f_params.FstepX = 2000
f_params.Weight2Len = 1

# Introducing coordinates
Wall_Dim(coord_dict, f_params)
Soil_Dim(coord_dict, f_params)
xE, yE = Wall(coord_dict)
xS, yS = Soil(coord_dict)

# Plot main features (wall + soil)
fig, ax = Plotter(xE, yE, xS, yS, params, f_params)

# Plot failure planes and compile failure surfaces for shoelace method later
Failure_Planes(coord_dict, f_params, params, show_FP = True, ax = ax)

# Determine area of each polygon
WeightFn(coord_dict, s_params, f_params)

Calc_Pp(coord_dict, f_params, s_params, params, 
        show_Wline = True, show_Pline = True, 
        show_Wpoint = True, show_Rpoint = True, show_Rline = True, show_min_Pp = True,
        ax = ax)

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())

ax.set_title('Culmann Graphical Approach (Passive Ratio: {}mm/kN)'.format(f_params.Weight2Len))
plt.show()


