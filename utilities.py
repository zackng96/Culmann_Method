import math
import matplotlib.pyplot as plt


def WeightFn(coords):
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
    return area  # in mm2

def xyGenerator(ListOf3DCoords):
    # unpack 3DCoords using loop and slicing
    x = []
    y = []
    for coord in ListOf3DCoords:
        x.append(coord[0])
        y.append(coord[1])

    return x,y
