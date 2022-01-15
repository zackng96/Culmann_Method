import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def Scale(unit):
    if unit == 'm':
        scale = .001
    
    elif unit == 'cm':
        scale = .01
        
    elif unit == 'mm':
        scale = 1
    
    elif unit == 'ft':
        scale = 0.00328084
        
    elif unit == 'kN/m':
        scale = 1
        
    elif unit == 'kips/ft':
        scale = 1/14.5939029372063648294
        
    return scale

def CoordGenerator(ListOf2DCoords):
    # unpack 3DCoords using loop and slicing
    x = []
    y = []
    for coord in ListOf2DCoords:
        x.append(coord[0])
        y.append(coord[1])
    return x, y

def Wall_Dim(coord_dict, f_params):
    WG = coord_dict['Wall_Geometry']
    WH = coord_dict['Wall_Ht']
    E1 = [WG[0], 0, 0]
    E2 = [WG[1], 0, 0]
    E3 = [WG[2], WH, 0]
    E4 = [WG[3], WH, 0]
    xBound, yBound = f_params.xBound, f_params.yBound
    coord_dict['Wall_Coordinates'] = [E1,E2,E3,E4]
    return E1, E2, E3, E4

def Wall(coord_dict):
    # Defines ERSS polygon (x and y coordinates)
    xE, yE = CoordGenerator(coord_dict['Wall_Coordinates'])
    coord_dict['Wall_X'] = xE
    coord_dict['Wall_Y'] = yE
    return xE, yE

def Soil_Dim(coord_dict, f_params):
    SoilLen_list = coord_dict['SoilLen_list']
    SoilHt_list = coord_dict['SoilHt_list']
    WG = coord_dict['Wall_Geometry']
    WH = coord_dict['Wall_Ht']
    xBound, yBound = f_params.xBound, f_params.yBound
    
    if len(SoilLen_list) != len(SoilHt_list):
        return 'Check input parameters'
    
    ERSSBotRight = WG[1]
    ERSSTopRight = WG[2]
    S1 = [ERSSBotRight, 0, 0]
    S2 = [ERSSTopRight, WH, 0]
    SoilSurfaces = [S1, S2]
    cum_SoilLen = 0
    cum_SoilHt = 0
    for SoilLen, SoilHt in zip(SoilLen_list,SoilHt_list):
        cum_SoilLen += SoilLen
        cum_SoilHt += SoilHt
        S = [ERSSBotRight + cum_SoilLen, WH + cum_SoilHt, 0]
        SoilSurfaces.append(S)
#    S3 = [ERSSBotRight + SoilLen1, ERSSHt + SoilHt1, 0]
#    S4 = [ERSSBotRight + SoilLen1 + SoilLen2, ERSSHt + SoilHt1 + SoilHt2, 0]
    # force close polygon
    Sm2 = [xBound, WH + cum_SoilHt, 0]
    Sm1 = [xBound, 0, 0]
    SoilSurfaces.append(Sm2)
    SoilSurfaces.append(Sm1)
    coord_dict['Soil_Coordinates'] = SoilSurfaces
    return SoilSurfaces

def Soil(coord_dict):
    # Defines soil polygon (x and y coordinates)
    xS, yS = CoordGenerator(coord_dict['Soil_Coordinates'])  
    coord_dict['Soil_X'] = xS
    coord_dict['Soil_Y'] = yS
    return xS, yS

def Plotter(xE, yE, xS, yS, params, f_params, line_load_pos = 0, line_load_mag = 0):
    scale = Scale(params.unit)
        
    xE = np.asarray(xE)
    yE = np.asarray(yE)
    xS = np.asarray(xS)
    yS = np.asarray(yS)
    
    fig, ax = plt.subplots(dpi = 100)
    ax.set_aspect('equal')
    ax.set_xlabel(params.unit)
    ax.set_ylabel(params.unit)
    ax.set_ylim([-1*f_params.yBound*.5*scale, f_params.yBound*scale])
    ax.set_xlim([0, f_params.xBound*scale])
    ax.fill(xE*scale, yE*scale, color=params.WallColor)  # wall
    ax.fill(xS*scale, yS*scale, color=params.SoilColor)  # soil
    


    return fig, ax
    
def Failure_Planes(coord_dict, f_params, params, show_FP = True, ax = None):
    scale = Scale(params.unit)
    # Plot failure planes and compile failure surfaces for shoelace method later
    Fdict = {}
    FArea_list = []
    SlipLen_list = []
    WG = coord_dict['Wall_Geometry']
    WH = coord_dict['Wall_Ht']
    ERSSBotRight = WG[1]
    ERSSTopRight = WG[2]
    SoilLen_list = coord_dict['SoilLen_list']
    SoilHt_list = coord_dict['SoilHt_list']
    SoilSurfaces = coord_dict['Soil_Coordinates']
    xBound = f_params.xBound

    SoilLen_array = np.insert(np.flip(np.cumsum(np.asarray(SoilLen_list))), len(SoilLen_list), 0) 
    SoilHt_array = np.insert(np.flip(np.cumsum(np.asarray(SoilHt_list))), len(SoilHt_list), 0)
    
    count = 1
    cum_SoilLen = SoilLen_array[0]
    cum_SoilHt = SoilHt_array[0]
    eqn_update = True
    for SlipLen in range(f_params.FendX, f_params.FstartX, -1*f_params.FstepX):
        SlipLen += ERSSBotRight
        pos_check = np.argwhere(SoilLen_array + ERSSBotRight <= SlipLen)
        
        if bool(list(pos_check)):
            i = pos_check[0][0]
        else:
            i = 0
        
        if cum_SoilLen != SoilLen_array[i]:
            cum_SoilLen = SoilLen_array[i] # update cum_SoilLen and call for eqn update
            cum_SoilHt = SoilHt_array[i] # update cum_SoilLen and call for eqn update
            eqn_update = True

        if eqn_update:
            # construct polyline params using boundary pts within range
            if i == 0:
                X = [SlipLen, ERSSBotRight + SoilLen_array[i+1]]
                Y = [WH + cum_SoilHt, WH + cum_SoilHt]
                
            elif i != (len(SoilLen_array)-1):
                X = [ERSSBotRight + cum_SoilLen, ERSSBotRight + SoilLen_array[i-1]]
                Y = [WH + cum_SoilHt, WH + SoilHt_array[i-1]]

            else:
                X = [ERSSTopRight, ERSSBotRight + SoilLen_array[i-1]]
                Y = [WH, WH + SoilHt_array[i-1]]
                
            coefficients = np.polyfit(X, Y, 1)
            m, c = coefficients[0], coefficients[1]
#            print('m = {0}, c = {1}, i = {2}'.format(m,c,i))
            eqn_update = False

        if SlipLen < xBound:
            Fdict["F{0}".format(count)] = [SlipLen, m*SlipLen + c, 0]
            SS = SoilSurfaces[0:(len(SoilSurfaces)-i-1)]
            FArea = SS + [Fdict["F{0}".format(count)]]
            FArea_list.append(FArea)
            SlipLen_list.append(SlipLen)
            
            if show_FP:
                x, y = CoordGenerator([SoilSurfaces[0], Fdict["F{0}".format(count)]])
                FP, = ax.plot(scale*np.asarray(x), scale*np.asarray(y), color=params.SlipPlaneColor, label='Failure Plane')
                
            count += 1
            
    coord_dict['Fdict'] = Fdict
    coord_dict['FArea_list'] = FArea_list
    coord_dict['SlipLen_list'] = SlipLen_list
    print("Total of {} failure slips considered".format(count-1))
    return FArea_list

def WeightFn(coord_dict, s_params, f_params, line_load_pos = 0, line_load_mag = 0):
    # line_load_mag is in kN
    Weight_list = []
    Wline_list = []
    gamma = s_params.UnitWeight
    Weight2Len = f_params.Weight2Len
    FArea_list = coord_dict['FArea_list']
    SlipLen_array = np.asarray(coord_dict['SlipLen_list'])
    line_load_pos = np.asarray(line_load_pos)
    line_load_mag = np.asarray(line_load_mag)
    
    for coords,SlipLen in zip(FArea_list,SlipLen_array):
        edges = len(coords)
        sum1 = 0
        sum2 = 0
    
        for i in range(0, edges - 1):
            sum1 = sum1 + coords[i][0] * coords[i + 1][1]
            sum2 = sum2 + coords[i][1] * coords[i + 1][0]
    
        # Add xn.y1
        sum1 = sum1 + coords[edges - 1][0] * coords[0][1]
        # Add x1.yn
        sum2 = sum2 + coords[0][0] * coords[edges - 1][1]
    
        area = abs(sum1 - sum2) / 2
        weight = area * gamma/1000/1000
        ll_weight = np.sum(line_load_mag*(line_load_pos < SlipLen))
        weight += ll_weight
        print('@ SlipLen = {0}mm, ll_weight = {1}, weight = {2}'.format(SlipLen, ll_weight, weight))
        Weight_list.append(weight)
        Wline_list.append((weight * Weight2Len))
    coord_dict['Wline_list'] = Wline_list
    return Weight_list  # in mm2

def WPline(coord_dict, s_params, params, p_type = 'Pa', show_Wline = True, show_Pline = True, ax = None):
    SoilSurfaces = coord_dict['Soil_Coordinates']
    WG = coord_dict['Wall_Geometry']
    WH = coord_dict['Wall_Ht']
    alpha = coord_dict['Wall_alpha']
    SoilHt_list = coord_dict['SoilHt_list']
    scale = Scale(params.unit)
    
    S1 = SoilSurfaces[0]
    ERSSBotRight = WG[1]
    phi = s_params.phi
    delta = s_params.delta
    
    if p_type == 'Pa':
        Wgradient = np.tan(phi / 180 * np.pi)
        WeightX = (WH + sum(abs(np.asarray(SoilHt_list)))) / Wgradient
        theta = 180 - alpha - delta - phi
        
        xW = [S1[0], ERSSBotRight + WeightX]
        yW = [S1[1], WH + sum(abs(np.asarray(SoilHt_list)))]
        
        PY = WeightX * np.tan(theta / 180 * np.pi)
        PX = WeightX
        
        PGradient = PY / PX
        xP = [S1[0], ERSSBotRight + PX]
        yP = [S1[1], -1 * PY]
    
    if p_type == 'Pp':
        Wgradient = np.tan(phi / 180 * np.pi)
        WeightX = (WH + sum(abs(np.asarray(SoilHt_list)))) / Wgradient
        theta = 180 - alpha - delta - phi
        
        xW = [S1[0], ERSSBotRight + WeightX]
        yW = [S1[1], -1*(WH + sum(abs(np.asarray(SoilHt_list))))]
        
        PY = WeightX * np.tan(theta / 180 * np.pi)
        PX = WeightX
        
        PGradient = PY / PX
        xP = [S1[0], ERSSBotRight + PX]
        yP = [S1[1], PY]
    if show_Wline:
        Wline, = ax.plot(scale*np.asarray(xW), scale*np.asarray(yW), color=params.WeightLineColor, label='Weight Line')
        
    if show_Pline:
        Pline, = ax.plot(scale*np.asarray(xP), scale*np.asarray(yP), color=params.PLineColor, label='P Line')
    return (xW, yW), (xP, yP)

def Calc_Pa(coord_dict, f_params, s_params, params,
            show_Wline = True, show_Pline = True,
            show_Wpoint = True, show_Rpoint = True, show_Rline = True, show_max_Pa = True,
            ax = None):
    SoilSurfaces = coord_dict['Soil_Coordinates']
    WG = coord_dict['Wall_Geometry']
    WH = coord_dict['Wall_Ht']
    alpha = coord_dict['Wall_alpha']
    SoilHt_list = coord_dict['SoilHt_list']
    FArea_list = coord_dict['FArea_list']
    Wline_list = coord_dict['Wline_list']
    Fdict = coord_dict['Fdict']
    scale = Scale(params.unit)
    fscale = Scale(params.force_unit)
    
    S1 = SoilSurfaces[0]
    ERSSBotRight = WG[1]
    phi = s_params.phi
    delta = s_params.delta

    WPline(coord_dict, s_params, params, p_type = 'Pa', show_Wline = show_Wline, show_Pline = show_Pline, ax = ax)
    
    Wcoords_list = []
    for Wline in Wline_list:
        WdeltaX = Wline * np.cos(phi / 180 * np.pi)
        WdeltaY = Wline * np.sin(phi / 180 * np.pi)
        Wcoords_list.append([S1[0] + WdeltaX, S1[1] + WdeltaY, S1[2]])
        if show_Wpoint:
            Wpoint = ax.scatter(scale*np.asarray([S1[0] + WdeltaX]), scale*np.asarray([S1[1] + WdeltaY]), s=params.WcoordsSize, zorder=3, color=params.WcoordsColor, label='W Point')
        
    # Using sine rule to determine length of R vector and subsequently, Rcoordslist
    count = 0
    Rline_list = []
    Rcoords_list = []
    for coord in FArea_list:
        kappa = np.arctan((coord[-1][1] - S1[1]) / (coord[-1][0] - S1[0])) * 180 / np.pi
        eta = kappa - phi
        theta = 180 - eta - (180 - alpha - delta)
        Rline = np.sin((180 - alpha - delta) / 180 * np.pi) * Wline_list[count] / np.sin(theta / 180 * np.pi)
        Rline_list.append(Rline)
        RdeltaX = Rline * np.cos(kappa / 180 * np.pi)
        RdeltaY = Rline * np.sin(kappa / 180 * np.pi)
        Rcoords_list.append([S1[0] + RdeltaX, S1[1] + RdeltaY, S1[2]])
        if show_Rpoint:
            Rpoint = ax.scatter(scale*np.asarray([S1[0] + RdeltaX]), scale*np.asarray([S1[1] + RdeltaY]), s=params.RcoordsSize, zorder=3, color=params.RcoordsColor, label='R Point')
        count += 1

    # Create smooth curve using interpolation fn
    xR, yR = CoordGenerator(Rcoords_list)
    interp_model = interp1d(xR, yR, kind = f_params.interp_kind)
    xR = np.linspace(np.asarray(xR).min(), np.asarray(xR).max(), f_params.interp_interval)
    yR = interp_model(xR)
    if show_Rline:
        Rline, = ax.plot(scale*np.asarray(xR), scale*np.asarray(yR), color=params.RcoordsColor, label='Resultant Curve')
    
    Distance_list = []
    for n in range(len(Wcoords_list)):
        Distance_list.append(np.sqrt(np.sum((Rcoords_list[n][0] - Wcoords_list[n][0]) ** 2
                                           + (Rcoords_list[n][1] - Wcoords_list[n][1]) ** 2)))
    if params.force_unit == 'kN/m':
        print("Pa is {0} kN/m".format(max(Distance_list) / f_params.Weight2Len))
    elif params.force_unit == 'kips/ft':
        print("Pa is {0} kips/ft".format(max(Distance_list)/f_params.Weight2Len * fscale))

    max_Pa = max(Distance_list)
    max_Pa_index = Distance_list.index(max_Pa)
    xD, yD = CoordGenerator([Rcoords_list[max_Pa_index], Wcoords_list[max_Pa_index]])
    if show_max_Pa:
        max_Pa_line, = ax.plot(scale*np.asarray(xD), scale*np.asarray(yD), color=params.max_Pa_Color, label='Maximum Pa = {}kN/m'.format(np.around(max_Pa/f_params.Weight2Len,2)))
        x, y = CoordGenerator([SoilSurfaces[0], Fdict["F{0}".format(max_Pa_index+1)]])
        max_FP, = ax.plot(scale*np.asarray(x), scale*np.asarray(y), color=params.final_FP_Color, label='Rupture Plane')
        
def Calc_Pp(coord_dict, f_params, s_params, params,
            show_Wline = True, show_Pline = True,
            show_Wpoint = True, show_Rpoint = True, show_Rline = True, show_min_Pp = True,
            ax = None):
    SoilSurfaces = coord_dict['Soil_Coordinates']
    WG = coord_dict['Wall_Geometry']
    WH = coord_dict['Wall_Ht']
    alpha = coord_dict['Wall_alpha']
    SoilHt_list = coord_dict['SoilHt_list']
    FArea_list = coord_dict['FArea_list']
    Wline_list = coord_dict['Wline_list']
    Fdict = coord_dict['Fdict']
    scale = Scale(params.unit)
    fscale = Scale(params.force_unit)
    
    S1 = SoilSurfaces[0]
    ERSSBotRight = WG[1]
    phi = s_params.phi
    delta = s_params.delta

    WPline(coord_dict, s_params, params, p_type = 'Pp', show_Wline = show_Wline, show_Pline = show_Pline, ax = ax)
    
    Wcoords_list = []
    for Wline in Wline_list:
        WdeltaX = Wline * np.cos(phi / 180 * np.pi)
        WdeltaY = Wline * np.sin(phi / 180 * np.pi)
        Wcoords_list.append([S1[0] + WdeltaX, -1*(S1[1] + WdeltaY), S1[2]])
        if show_Wpoint:
            Wpoint = ax.scatter(scale*np.asarray([S1[0] + WdeltaX]), -1*scale*np.asarray([S1[1] + WdeltaY]), s=params.WcoordsSize, zorder=3, color=params.WcoordsColor, label='W Point')
        
    # Using sine rule to determine length of R vector and subsequently, Rcoordslist
    count = 0
    Rline_list = []
    Rcoords_list = []
    sorting_list = []
    for coord in FArea_list:
        
        kappa = np.arctan((coord[-1][1] - S1[1]) / (coord[-1][0] - S1[0])) * 180 / np.pi
        eta = kappa + phi
        theta = 180 - eta - (180 - alpha + delta)
        Rline = np.sin((180 - alpha + delta) / 180 * np.pi) * Wline_list[count] / np.sin(theta / 180 * np.pi)
        Rline_list.append(Rline)
        RdeltaX = Rline * np.cos(kappa / 180 * np.pi)
        RdeltaY = Rline * np.sin(kappa / 180 * np.pi)
        Rcoords_list.append([S1[0] + RdeltaX, S1[1] + RdeltaY, S1[2]])
        if show_Rpoint:
            Rpoint = ax.scatter(scale*np.asarray([S1[0] + RdeltaX]), scale*np.asarray([S1[1] + RdeltaY]), s=params.RcoordsSize, zorder=3, color=params.RcoordsColor, label='R Point')
        count += 1
    
    # Create smooth curve using interpolation fn
    xR, yR = CoordGenerator(Rcoords_list)
    interp_model = interp1d(xR, yR, kind = f_params.interp_kind)
    # xR = np.linspace(np.asarray(xR).min(), np.asarray(xR).max(), f_params.interp_interval)
    # yR = interp_model(xR)
    if show_Rline:
        Rline, = ax.plot(scale*np.asarray(xR), scale*np.asarray(yR), color=params.RcoordsColor, label='Resultant Curve')
    
    Distance_list = []
    for n in range(len(Wcoords_list)):
        Distance_list.append(np.sqrt(np.sum((Rcoords_list[n][0] - Wcoords_list[n][0]) ** 2
                                           + (Rcoords_list[n][1] - Wcoords_list[n][1]) ** 2)))

    if params.force_unit == 'kN/m':
        print("Pp is {0} kN/m".format(min(Distance_list) / f_params.Weight2Len))
    elif params.force_unit == 'kips/ft':
        print("Pp is {0} kips/ft".format(min(Distance_list) / f_params.Weight2Len * fscale))

    min_Pp = min(Distance_list)
    min_Pp_index = Distance_list.index(min_Pp)
    xD, yD = CoordGenerator([Rcoords_list[min_Pp_index], Wcoords_list[min_Pp_index]])
    if show_min_Pp:
        min_Pp_line, = ax.plot(scale*np.asarray(xD), scale*np.asarray(yD), color=params.min_Pp_Color, label='Minimum Pp = {}kN/m'.format(np.around(min_Pp/f_params.Weight2Len,2)))
        x, y = CoordGenerator([SoilSurfaces[0], Fdict["F{0}".format(min_Pp_index+1)]])
        min_FP, = ax.plot(scale*np.asarray(x), scale*np.asarray(y), color=params.final_FP_Color, label='Rupture Plane')
        
# Provides all colour settings and configurations
class Display_Parameters:
    def __init__(self):
        self.WallColor = 'xkcd:dark grey'
        self.SoilColor = 'xkcd:light yellow'
        self.SlipPlaneColor = 'xkcd:grey'
        self.WeightLineColor = 'xkcd:teal'
        self.PLineColor = 'k'
        self.WcoordsColor = 'xkcd:lime green'
        self.WcoordsSize = 10
        self.RcoordsColor = 'xkcd:royal blue'
        self.RcoordsSize = 10
        self.max_Pa_Color = 'xkcd:pink'
        self.min_Pp_Color = 'xkcd:pink'
        self.final_FP_Color = 'xkcd:red'
        self.unit = 'mm'
        self.force_unit = 'kN/m'


class Failure_Parameters:
    def __init__(self):
        self.xBound = 10000
        self.yBound = 10000
        self.FstartX = 0
        self.FendX = 10000  # Furthest daylight point for failure slip
        self.FstepX = 2000  # Step size. Decrease for granularity
        self.Weight2Len = 10
        self.interp_kind = 'linear'
        self.interp_interval = 500


class Sand_Parameters:
    def __init__(self):
        self.phi = 33 # friction angle
        self.delta = 22  # wall friction
        self.UnitWeight = 18  # kN/m3


class ERSS:
    def __init__(self):
        self.ERSSHt = 6000
        self.alpha = 90
        self.ERSScoords = [0, 1000, 3000 - self.ERSSHt * np.tan((self.alpha - 90) / 180 * np.pi), 3000]  # from bottom left hand corner in clockwise direction
