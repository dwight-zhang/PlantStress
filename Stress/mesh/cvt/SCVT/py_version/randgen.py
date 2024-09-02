"""randgen.py
define:
    - function:
        :set_randome_seed(nrank): set seed for random 
        :uniform_random(): generate one point on unit sphere randomly
        :random_generator(dens_max): genrate one point on unit sphere randomly
        accoding to the density function
author:
    ZHANG Donghang
    zdh@lsec.cc.ac.cn
"""

import numpy as np
import time
import dens_f

def set_random_seed(nrank):
    '''set the seed for the random number generator
    
    :param nrank: int
    '''
    np.random.seed() 
    k = np.random.randint(1, 100) 
    seed = np.zeros(k, dtype=int)

    date_time = time.localtime()
    date_time_values = [
        date_time.tm_year, date_time.tm_mon, date_time.tm_mday,
        date_time.tm_hour, date_time.tm_min, date_time.tm_sec,
        date_time.tm_wday, date_time.tm_yday, date_time.tm_isdst
    ]

    for i in range(k):
        seed[i] = date_time_values[7 - (i % 8)] + (i + 1) * (nrank + 1) * 100

    np.random.seed(seed[0]) 

def uniform_random():
    '''generate a point on the unit sphere randomely

    :return: numpy.ndarray
    '''
    while True:
        crd = np.random.rand(3) * 2.0 - 1.0
        dd = np.sum(crd**2)

        if dd <= 1.0:
            dd = np.sqrt(dd)
            pt = crd / dd
            return pt


def random_generator(dens_max):
    '''sample a point on the unit sphere randomely accoding to the density
    function dens_f using the Selection Method

    :param dens_max: float
    :return: numpy.ndarray
    '''
    pt = np.zeros(3)

    while True:
        pt = uniform_random()
        u = np.random.random()
        if u <= (dens_f(pt) / dens_max):
            break

    return pt
