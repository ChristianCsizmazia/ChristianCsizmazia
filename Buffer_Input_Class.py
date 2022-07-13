# import pandas as pd
import numpy as np
# import os
# import logging
# import threading
# import time
import math

# from matplotlib import pyplot as plt
# from sympy import symbols
# from numba import jit
# from datetime import datetime

class buffer_input(object):
    def __init__(self, c, pH, pKa):
        # Buffer
        self.c = c
        self.pH = pH
        self.pKa = pKa

    def find_ratio(self, a, b):  #
        return 10 ** (a - b)

    def negative_to_power_ten(self, a):
        return 10 ** (-a)

    def henderson_hasselbalch(self, a, b, c):
        return a - math.log10(b / c)

    def calculate(self):
        Kw = 1e-14
        Ka = self.negative_to_power_ten(self.pKa)
        x = self.find_ratio(self.pH, self.pKa)

        A = x * self.c / (1 + x)
        HA = self.c - A

        H = self.negative_to_power_ten(self.pH)
        OH = Kw / H

        return [H, OH, Kw, HA, A, Ka]
        

