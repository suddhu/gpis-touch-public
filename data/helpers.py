#!/usr/bin/env python

# Sudharshan Suresh <suddhu@cmu.edu>, Nov 2019
# Helper functions for get_contacts.py

# Copyright: This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License v3 as published by
# the Free Software Foundation. This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of any FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License v3 for more details. You should have received a copy of the GNU General Public License v3
# along with this program; if not, you can access it online at http://www.gnu.org/licenses/gpl-3.0.html.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pdb
import math
import time 
import scipy.io as sio
matplotlib.use('TkAgg')

# generate the ground truth square object ([4,2])
def gen_rect1():
    data = sio.loadmat('shapes/rect/rect1.mat')
    vertices = data['rect1_shape']
    return vertices

def gen_circle():
    data = sio.loadmat('shapes/circle/circle.mat')
    vertices = data['circle_shape']
    return vertices

def gen_ellip2():
    data = sio.loadmat('shapes/ellip/ellip2.mat')
    vertices = data['ellip2_shape']
    return vertices

def gen_butter():
    data = sio.loadmat('shapes/butter/butter.mat')
    vertices = data['butter_shape']
    return vertices

# https://stackoverflow.com/questions/23411688/drawing-polygon-with-n-number-of-sides-in-python-3-2
def genInit(sides, radius=1, rotation=0, translation=None):
    one_segment = math.pi * 2 / sides

    points = [
        (math.sin(one_segment * i + rotation) * radius,
         math.cos(one_segment * i + rotation) * radius)
        for i in range(sides)]

    if translation:
        points = [[sum(pair) for pair in zip(point, translation)]
                  for point in points]

    return np.asarray(points)

# take user input ordered vertices (convex shape)
def getContacts(): 
    contactsAndNormals = np.asarray(plt.ginput(0,0))

    contacts =  contactsAndNormals[0::2,:]
    normals = contactsAndNormals[1::2,:] - contactsAndNormals[0::2,:]

    return contacts, normals

def getDateTime():
    return time.strftime("%Y%m%d-%H%M")
