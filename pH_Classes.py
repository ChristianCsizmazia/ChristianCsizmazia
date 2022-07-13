# import pandas as pd
import numpy as np
# import os
# import logging
# import threading
# import time
import math

from Buffer_Input_Class import buffer_input
# from matplotlib import pyplot as plt
# from sympy import symbols
# from numba import jit
# from datetime import datetime

class one_buffer(object):
    def __init__(self, H, OH, Kw, HA, A, Ka):
        # General ions
        self.H = H
        self.OH = OH
        self.Kw = Kw

        # Buffer A
        self.HA = HA
        self.A = A
        self.Ka = Ka

    def find_y(self, c, x, d):
        return (c - x * self.Ka * self.HA) / (x * self.Ka + d)

    def find_xy(self):
        H = self.H
        OH = self.OH
        Kw = self.Kw

        HA = self.HA
        A = self.A
        Ka = self.Ka

        c = Ka * (OH * HA) - Kw * A
        d = -1 * Kw - Ka * OH

        x3 = Ka
        x2 = -1 * H * Ka - OH * Ka + d - HA * Ka
        x1 = -1 * OH * d + H * d + Ka * H * OH + Ka * Kw + OH * Ka * HA + c
        x0 = H * OH * d - Kw * d - OH * c
        x_123 = np.roots([x3, x2, x1, x0])

        y_123 = []
        for i in range(len(x_123)):
            y_123.append(self.find_y(c, x_123[i], d))
        return [x_123, y_123]

    def h_calculation1(self, x):
        return self.Kw / (self.OH - x)

    def h_calculation2(self, y):
        return (self.Ka * (self.HA + y)) / (self.A - y)

    def check_positive(self, a):
        if a >= 0:
            return True
        else:
            return False

    def check_and_append(self, list1, list2, hc1, hc2):
        if self.check_positive(hc1) != True or isinstance(hc1, complex) != False:
            return
        list1.append(hc1)
        list2.append(abs(hc1 - hc2))

    def find_pH(self):
        local_array = self.find_xy()
        calc_h = []
        reasonable_h = []
        list_check = []

        for i in range(len(local_array[0])):
            calc_h.append(self.h_calculation1(local_array[0][i]))

        for i in range(len(calc_h)):
            self.check_and_append(reasonable_h, list_check, calc_h[i], self.h_calculation2(local_array[1][i]))

        if len(reasonable_h) != 0:
            correct_h = [x for _, x in sorted(zip(list_check, reasonable_h))]
            return -1 * math.log10(correct_h[0])
        else:
            return np.nan

class two_buffers(one_buffer):
    def __init__(self, H, OH, Kw, HA, A, Ka, HB, B, Kb):
        one_buffer.__init__(self, H, OH, Kw, HA, A, Ka)

        # Buffer B
        self.HB = HB
        self.B = B
        self.Kb = Kb

    def h_calculation3(self, z):  # unused
        return (self.Kb * (self.HB + z)) / (self.B - z)

    def find_xy(self):
        H = self.H
        OH = self.OH
        Kw = self.Kw

        HA = self.HA
        A = self.A
        Ka = self.Ka

        HB = self.HB
        B = self.B
        Kb = self.Kb

        c = Ka * (OH * HA) - Kw * A
        d = -1 * Kw - Ka * OH
        e = Kb * (OH * HB) - Kw * B
        f = -1 * Kw - Kb * B
        g = H * OH - Kw

        x4 = Ka * Kb
        x3 = (f * Ka + d * Kb - (Ka * Kb) * (H + OH + HA + HB))
        x2 = (d * f - Ka * (f * H + f * OH + f * HA - e) - Kb * (d * H + d * OH + d * HB - c) + (Ka * Kb) * (g + HA * OH + HB * OH))
        x1 = (c * f + d * e - d * f * (H + OH)) + Ka * (f * g + f * HA * OH - e * OH - e * OH) + Kb * (d * g + d * HB * OH - c * OH)
        x0 = d * f * g - c * f * OH - d * e * OH
        x_1234 = np.roots([x4, x3, x2, x1, x0])

        y_1234 = []
        for i in range(len(x_1234)):
            y_1234.append(self.find_y(c, x_1234[i], d))
        return [x_1234, y_1234]

"""
class three_buffers(two_buffers):  # unused, untested!
    def __init__(self, H, OH, Kw, HA, A, Ka, HB, B, Kb, HC, C, Kc):
        two_buffers.__init__(self, H, OH, Kw, HA, A, Ka, HB, B, Kb)

        # Buffer C
        self.HC = HC
        self.C = C
        self.Kc = Kc

    def h_calculation4(self, u):
        return (self.Kc * (self.HC + u)) / (self.C - u)
    
    def find_xy(self):
        H = self.H
        OH = self.OH
        Kw = self.Kw

        HA = self.HA
        A = self.A
        Ka = self.Ka

        HB = self.HB
        B = self.B
        Kb = self.Kb

        HC = self.HC
        C = self.C
        Kc = self.Kc

        c = Ka * (OH * HA) - Kw * A
        d = -1 * Kw - Ka * OH
        e = Kb * (OH * HB) - Kw * B
        f = -1 * Kw - Kb * B
        g = Kc * OH * HC - Kw * C
        i = - 1 * Kw * -1 * Kb * OH

        x5 = Ka * Kb * Kc
        x4 = -1 * H * Ka * Kc - 1 * OH * Ka * Kb * Kc + f * Ka * Kc + d * Kb * Kc + i * Ka * Kb - 1 * HA * Ka * Kb * Kc - 1 * HB * Ka * Kb * Kc - 1 * HC * Ka * Kb * Kc
        x3 = H * OH * Ka * Kb * Kc - 1 * H * f * Ka * Kc - 1 * H * d * Kb * Kc - 1 * H * i * Ka * Kb -1 * OH * f * Ka * Kc -1 * OH * d * Kb * Kc -1 * OH * i * Ka * Kb + d * f * Kc + f * i * Ka + d * i * Kb - 1 * Kw * Ka * Kb * Kc + OH * HA * Ka * Kb * Kc + c * Kb * Kc - 1* HA * i * Ka * Kb -1 * HA * f * Ka * Kc + OH * HB * Ka * Kb* Kc + e * Ka * Kc - 1 * HB * i * Ka * Kb - 1 * HB * d * Kb * Kc + OH * HC * Ka * Kb * Kc + g * Ka * Kb - 1 * HC * d * Kb * Kc - 1 * HC * f * Ka * Kc
        x2 = H * OH * f * Ka * Kc + H * OH * d * Kb * Kc + H * OH * i * Ka * Kb - 1 * H * d * f * Kc -1 * H * f * Ka -1 * H * d * i * Kb -1 * OH * d * f * Kc -1 * OH * f * i * Ka -1 * OH * d * i * Kb -1 * f * Kw * Ka * Kc -1 * d * Kw * Kb * Kc -1 * i * Kw * Ka * Kb + d * f * i -1 * OH * c * Kb * Kc + OH * HA * i * Ka * Kb + c * i * Kb + OH * HA * f * Ka * Kc + c * f * Kc - 1 * HA * f * i * Ka -1 * OH * e * Ka * Kc + OH * HB * i * Ka * Kb + e * i * Ka + OH * HB * d * Kb * Kc + d * e * Kc - HB * d * i * Kb - OH * g * Ka * Kb + OH * HC * d * Kb * Kc + d * g * Kb + OH * HC * f * Ka * Kc + f * g * Ka - 1 * HC * d * f * Kc
        x1 = H * OH * d * f * Kc + H * OH * f * i * Ka + H * OH * d * i * Kb - H * d * f * i -1 * OH * d * f * i -1 * d * f * Kw * Kc - f * i * Kw * Ka -1 * d * i * Kw * Kb -1 * OH * c * i * Kb -1 * OH * c * f * Kc + OH * HA  * f * i * Ka + c * f * i -1 * OH * e * i * Ka -1 * OH * e * d * Kc + OH * HB * d * i * Kb + d * i * e -1 * OH * d * g * Kb - OH * g * f * Ka+ OH * HC * d * f * Kc + d * f * g
        x0 = -1 * d * f * i * Kw - 1 * OH * c * f * i -1 * OH * d * e * i - OH * d * f * g + H * OH * d * f * i
        x_12345 = np.roots([x5, x4, x3, x2, x1, x0])

        y_12345 = []
        for i in range(len(x_12345)):
            y_12345.append(self.find_y(c, x_12345[i], d))
        return [x_12345, y_12345]
"""

