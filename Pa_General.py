import numpy as np
from utilities import *
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Define ERSS position and parameters
SoilHt1 = 3500  # elevation of 1st slope peak above wall top
SoilLen1 = 5000  # length of which 1st slope occurs
SoilHt2 = 2500  # elevation of 2nd slope peak above 1st peak
SoilLen2 = 10000  # length of which 2nd slope occurs
ERSSHt = 6000
alpha = 90
ERSSTopRight = 3000 - ERSSHt * math.tan((alpha - 90) / 180 * math.pi)
ERSSTopLeft = 1000
ERSSBotRight = 3000
ERSSBotLeft = 0
xBound = 30000
yBound = 10000  # bottom bound set to 10% buffer by default

phi = 33
c = 0  # assumed cohesionless
delta = 22
gamma = 18  # soil unit weight in kN/m3

# Trial failure slip control parameters
FstartX = 0
FendX = 20000  # Furthest daylight point for failure slip
FstepX = 500  # Step size. Decrease for granularity

# Colors
WallColor = 'xkcd:grey'
SoilColor = 'xkcd:dull yellow'
SlipPlaneColor = 'r'
WeightLineColor = 'xkcd:teal'
PLineColor = 'xkcd:royal blue'
WcoordsColor = 'xkcd:lime green'
RcoordsColor = 'k'
max_Pa_Color = 'xkcd:turquoise'

# Point Sizes
WcoordsSize = 10
RcoordsSize = 10

# Conversion parameter
Weight2Len = 15  # mm/kN
interp_kind = 'cubic'
interp_interval = 500

# Introducing coordinates
E1 = [ERSSBotLeft, 0, 0]
E2 = [ERSSBotRight, 0, 0]
E3 = [ERSSTopRight, ERSSHt, 0]
E4 = [ERSSTopLeft, ERSSHt, 0]
xE, yE = xyGenerator([E1, E2, E3, E4])

S1 = E2
S2 = E3
S3 = [ERSSBotRight + SoilLen1, ERSSHt + SoilHt1, 0]
S4 = [ERSSBotRight + SoilLen1 + SoilLen2, ERSSHt + SoilHt1 + SoilHt2, 0]
S5 = [xBound, ERSSHt + SoilHt1 + SoilHt2, 0]
S6 = [xBound, 0, 0]
xS, yS = xyGenerator([S1, S2, S3, S4, S5, S6])

# Plot main features (wall + soil)
plt.figure(figsize=(8, 8))
plt.axis('equal')
plt.xlabel('mm')
plt.ylabel('mm')
plt.fill(xE, yE, color=WallColor)  # wall
plt.fill(xS, yS, color=SoilColor)  # soil

# Plot failure planes and compile failure surfaces for shoelace method later
Fdict = {}
FArea = []
count = 0
for SlipLen in range(FstartX, FendX, FstepX):
    if SlipLen > ERSSBotRight + SoilLen1 + SoilLen2:
        Fdict["F{0}".format(count)] = [SlipLen, ERSSHt + SoilHt1 + SoilHt2, 0]
        x, y = xyGenerator([S1, Fdict["F{0}".format(count)]])
        FS, = plt.plot(x, y, color=SlipPlaneColor, label='Failure Plane')
        FArea.append([S1, Fdict["F{0}".format(count)], S4, S3, S2])

    elif SlipLen > ERSSBotRight + SoilLen1:
        delX = SlipLen - ERSSBotRight - SoilLen1
        delY = delX / SoilLen2 * SoilHt2
        Fdict["F{0}".format(count)] = [SlipLen, ERSSHt + SoilHt1 + delY, 0]
        x, y = xyGenerator([S1, Fdict["F{0}".format(count)]])
        FS, = plt.plot(x, y, color=SlipPlaneColor, label='Failure Plane')
        FArea.append([S1, Fdict["F{0}".format(count)], S3, S2])

    elif SlipLen > ERSSBotRight:
        delX = SlipLen - ERSSBotRight
        delY = delX / SoilLen1 * SoilHt1
        Fdict["F{0}".format(count)] = [SlipLen, ERSSHt + delY, 0]
        x, y = xyGenerator([S1, Fdict["F{0}".format(count)]])
        FS, = plt.plot(x, y, color=SlipPlaneColor, label='Failure Plane')
        FArea.append([S1, Fdict["F{0}".format(count)], S2])

    count += 1
print("Total of {} failure slips considered".format(count + 1))

# Plot weight line (angled at phi to horizontal)
Wgradient = math.tan(phi / 180 * math.pi)
WeightX = (ERSSHt + SoilHt1 + SoilHt2) / Wgradient
xW = [S1[0], ERSSBotRight + WeightX]
yW = [S1[1], ERSSHt + SoilHt1 + SoilHt2]
Wline, = plt.plot(xW, yW, color=WeightLineColor, label='Weight Line')

# Plot P line (angled at pi-alpha-delta away from horizontal)
theta = 180 - alpha - delta - phi
PY = WeightX * math.tan(theta / 180 * math.pi)
PX = WeightX
PGradient = PY / PX
xP = [S1[0], ERSSBotRight + PX]
yP = [S1[1], -1 * PY]
Pline, = plt.plot(xP, yP, color=PLineColor, label='P Line')

# Determine area of each polygon using shoelace method to find W length
# Firstly, define area function that takes in list of coords in FArea to compute wedge weight: Use WeightFn
WeightList = []
for coords in FArea:
    WeightList.append(WeightFn(coords) * (gamma / 1000 / 1000))

# Adopt suitable scale by trial and error and convert WeightList to equivalent length
WlineList = [i * Weight2Len for i in WeightList]

# W1 to Wn using WlineList
WcoordsList = []
for Wline in WlineList:
    WdeltaX = Wline * math.cos(phi / 180 * math.pi)
    WdeltaY = Wline * math.sin(phi / 180 * math.pi)
    WcoordsList.append([S1[0] + WdeltaX, S1[1] + WdeltaY, S1[2]])
    Wpoint = plt.scatter([S1[0] + WdeltaX], [S1[1] + WdeltaY], s=WcoordsSize, zorder=3, color=WcoordsColor, label='W Point')

# Using sine rule to determine length of R vector and subsequently, Rcoordslist
count = 0
RlineList = []
RcoordsList = []
for coord in FArea:
    kappa = math.atan((coord[1][1] - S1[1]) / (coord[1][0] - S1[0])) * 180 / math.pi
    eta = kappa - phi
    theta = 180 - eta - (180 - alpha - delta)
    Rline = math.sin((180 - alpha - delta) / 180 * math.pi) * WlineList[count] / math.sin(theta / 180 * math.pi)
    RlineList.append(Rline)
    RdeltaX = Rline * math.cos(kappa / 180 * math.pi)
    RdeltaY = Rline * math.sin(kappa / 180 * math.pi)
    RcoordsList.append([S1[0] + RdeltaX, S1[1] + RdeltaY, S1[2]])
    Rpoint = plt.scatter([S1[0] + RdeltaX], [S1[1] + RdeltaY], s=RcoordsSize, zorder=3, color=RcoordsColor, label='R Point')
    count += 1

# Create smooth curve using interpolation fn
xR, yR = xyGenerator(RcoordsList)
interp_model = interp1d(xR, yR, kind=interp_kind)
xR = np.linspace(np.asarray(xR).min(), np.asarray(xR).max(), interp_interval)
yR = interp_model(xR)
Rline, = plt.plot(xR, yR, color=RcoordsColor, label='Resultant Curve')

Distancelist = []
for n in range(len(WcoordsList)):
    Distancelist.append(np.sqrt(np.sum((RcoordsList[n][0] - WcoordsList[n][0]) ** 2
                                       + (RcoordsList[n][1] - WcoordsList[n][1]) ** 2)))
print("Pa is {0} kN/m".format(max(Distancelist) / Weight2Len))

max_Pa = max(Distancelist)
max_Pa_index = Distancelist.index(max_Pa)
xD, yD = xyGenerator([RcoordsList[max_Pa_index], WcoordsList[max_Pa_index]])
max_Pa_line, = plt.plot(xD, yD, color=max_Pa_Color, label='Maximum Pa = {}kN/m'.format(np.around(max_Pa/Weight2Len,2)))

#remove duplicates in legend
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())

plt.title('Culmann Graphical Approach (Active, Ratio: 15mm/kN)')
plt.xlim((-.1 * xBound, 1.1 * xBound))
plt.ylim((-.1 * yBound, yBound))
plt.show()
