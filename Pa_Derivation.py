import rhinoscriptsyntax as rs
import math
from System.Drawing import Color
rs.ClearCommandHistory()
rs.DeleteObjects(rs.AllObjects(select = True))

#Define ERSS position and parameters
rs.AddLayer("Layer 1")
rs.CurrentLayer(layer = "Layer 1")
rs.LayerColor("Layer 1",Color.Black)
SoilHt1 = 3500 #elevation of 1st slope peak above wall top
SoilLen1 = 5000 #length of which 1st slope occurs
SoilHt2 = 2500 #elevation of 2nd slope peak above 1st peak 
SoilLen2 = 10000 #length of which 2nd slope occurs
ERSSHt = 6000
alpha = 90
ERSSTopRight = 3000 - ERSSHt*math.tan((alpha-90)/180*math.pi)
ERSSTopLeft = 1000
ERSSBotRight =3000
ERSSBotLeft = 0
Bound = 50000

phi = 33
c = 0
delta = 22
gamma = 18
E1 = [ERSSBotLeft,0,0]
E2 = [ERSSBotRight,0,0]
E3 = [ERSSTopRight,ERSSHt,0]
E4 = [ERSSTopLeft,ERSSHt,0]
rs.AddLine(E1,E2)
rs.AddLine(E2,E3)
rs.AddLine(E3,E4)
rs.AddLine(E4,E1)

S1 = E2
S2 = E3
S3 = [ERSSBotRight+SoilLen1,ERSSHt+SoilHt1,0]
S4 = [ERSSBotRight+SoilLen1+SoilLen2,ERSSHt+SoilHt1+SoilHt2,0]
S5 = [Bound,ERSSHt+SoilHt1+SoilHt2,0]
S6 = [Bound,0,0]
rs.AddLine(S1,S2)
rs.AddLine(S2,S3)
rs.AddLine(S3,S4)
rs.AddLine(S4,S5)
rs.AddLine(S5,S6)
rs.AddLine(S6,S1)

#Plot failure planes and compile failure surfaces for shoelace method later
rs.AddLayer("Layer 2")
rs.CurrentLayer(layer = "Layer 2")
rs.LayerColor("Layer 2", Color.Blue)
Fdict = {}
FArea = []
FstartX = 0  
FendX = 20000  #iterate
FstepX = 1000  #iterate
count = 0
for SlipLen in range(FstartX, FendX, FstepX):
    if SlipLen > ERSSBotRight+SoilLen1+SoilLen2:
        Fdict["F{0}".format(count)] = [SlipLen,ERSSHt+SoilHt1+SoilHt2,0]
        rs.AddLine(S1,Fdict["F{0}".format(count)])
        FArea.append([S1, Fdict["F{0}".format(count)], S4, S3, S2])
        
    elif SlipLen > ERSSBotRight+SoilLen1:
        delX = SlipLen - ERSSBotRight - SoilLen1
        delY = delX/SoilLen2 * SoilHt2
        Fdict["F{0}".format(count)] = [SlipLen,ERSSHt+SoilHt1+delY,0]
        rs.AddLine(S1,Fdict["F{0}".format(count)])
        FArea.append([S1, Fdict["F{0}".format(count)], S3, S2])
        
    elif SlipLen > ERSSBotRight:
        delX = SlipLen - ERSSBotRight
        delY = delX/SoilLen1 * SoilHt1
        Fdict["F{0}".format(count)] = [SlipLen,ERSSHt+delY,0]
        rs.AddLine(S1,Fdict["F{0}".format(count)])
        FArea.append([S1, Fdict["F{0}".format(count)], S2])
    
    #Fdict["F{0}".format(count)] = [ERSSBotRight+(SoilHt-ERSSHt)+i,SoilHt,0]
    #rs.AddLine(S1,Fdict["F{0}".format(count)])
    #if Fdict["F{0}".format(count)][0] <= ERSSBotRight + (SoilHt-ERSSHt):
    #    FArea.append([S1, Fdict["F{0}".format(count)], S2])
    #else:
    #    FArea.append([S1, Fdict["F{0}".format(count)], S3, S2])
    count += 1
print("Total of {} failure slips considered".format(count+1))

#Plot weight line (angled at phi to horizontal)
rs.AddLayer("Layer 3")
rs.CurrentLayer(layer = "Layer 3")
rs.LayerColor("Layer 3",Color.Red)
Wgradient = math.tan(phi/180*math.pi)
WeightX =(ERSSHt+SoilHt1+SoilHt2)/Wgradient
rs.AddLine(S1,[ERSSBotRight+WeightX,ERSSHt+SoilHt1+SoilHt2,0])

#Plot P line (angled at pi-alpha-delta away from horizontal)
rs.AddLayer("Layer 4")
rs.CurrentLayer(layer = "Layer 4")
rs.LayerColor("Layer 4",Color.Purple)
theta = 180 - alpha - delta - phi
PY = WeightX*math.tan(theta/180*math.pi)
PX = WeightX
rs.AddLine(S1,[ERSSBotRight+PX,-1*PY,0])
PGradient = PY/PX

#Determine area of each polygon using shoelace method to find W length
#Firstly, define area function that takes in list of coords in FArea to compute wedge weight
def WeightFn(coords):
    edges = len(coords)
    sum1 = 0
    sum2 = 0
  
    for i in range(0,edges-1):
        sum1 = sum1 + coords[i][0] *  coords[i+1][1]
        sum2 = sum2 + coords[i][1] *  coords[i+1][0]
  
    #Add xn.y1
    sum1 = sum1 + coords[edges-1][0]*coords[0][1]   
    #Add x1.yn
    sum2 = sum2 + coords[0][0]*coords[edges-1][1]   
  
    area = abs(sum1 - sum2) / 2
    return area #in mm2

WeightList = []
for coords in FArea:
    WeightList.append(WeightFn(coords)*(gamma/1000/1000))


#Adopt suitable scale by trial and error and convert WeightList to equivalent length
Weight2Len = 15 #mm
WlineList = [i*Weight2Len for i in WeightList]

#W1 to Wn using WlineList
WcoordsList = []
for Wline in WlineList:
    WdeltaX = Wline*math.cos(phi/180*math.pi)
    WdeltaY = Wline*math.sin(phi/180*math.pi)
    WcoordsList.append([S1[0]+WdeltaX,S1[1]+WdeltaY,S1[2]])
    rs.AddPoint([S1[0]+WdeltaX,S1[1]+WdeltaY,S1[2]])

#Using sine rule to determine length of R vector and subsequently, Rcoordslist
rs.AddLayer("Layer 5")
rs.CurrentLayer(layer = "Layer 5")
rs.LayerColor("Layer 5",Color.Yellow)
count = 0
RlineList = []
RcoordsList = []
for coord in FArea:
    kappa = math.atan((coord[1][1] - S1[1])/(coord[1][0]-S1[0]))*180/math.pi
    eta = kappa - phi
    theta = 180 - eta - (180 - alpha - delta)
    Rline = math.sin((180-alpha-delta)/180*math.pi)*WlineList[count]/math.sin(theta/180*math.pi)
    RlineList.append(Rline)
    RdeltaX = Rline*math.cos(kappa/180*math.pi)
    RdeltaY = Rline*math.sin(kappa/180*math.pi)
    RcoordsList.append([S1[0]+RdeltaX,S1[1]+RdeltaY,S1[2]])
    rs.AddPoint([S1[0]+RdeltaX,S1[1]+RdeltaY,S1[2]])
    count += 1
rs.AddCurve(RcoordsList, degree = 3)

rs.CurrentLayer(layer = "Layer 4")
Distancelist = []
for n in range(len(WcoordsList)):
    rs.AddLine(RcoordsList[n],WcoordsList[n])
    Distancelist.append(rs.Distance(RcoordsList[n],WcoordsList[n]))
print("Pa is {0} kN/m".format(max(Distancelist)/Weight2Len))

rs.CurrentLayer(layer = "Layer 4")
