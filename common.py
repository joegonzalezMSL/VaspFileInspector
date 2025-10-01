# -*- coding: utf-8 -*-

import math
import sys

D2R = 3.141592653589793238462643/180.0
R2D = 180.0/3.141592653589793238462643

class Point:
    def __init__(self):
        x = 0
        y = 0
        z = 0
        xs = 0
        ys = 0
        zs = 0
        globID = 0
        typeID = 0 

        rnn = []
        naborImg = []
        naborIds = []
        naborTypes = []