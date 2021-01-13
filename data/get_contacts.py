#!/usr/bin/env python

# Sudharshan Suresh <suddhu@cmu.edu>, Nov 2019
# Takes a series of input contacts and normal directions and saves them to txt file

# Copyright: This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License v3 as published by
# the Free Software Foundation. This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of any FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License v3 for more details. You should have received a copy of the GNU General Public License v3
# along with this program; if not, you can access it online at http://www.gnu.org/licenses/gpl-3.0.html.

from helpers import *
matplotlib.use('TkAgg')
import os
import argparse

# RECT1_SIDE = 0.09/2.0
# CIRCLE_RADIUS = 0.0525
# ELLIP2_R1 = 0.105/2.0 
# ELLIP2_R2 = 0.13089/2.0

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--shape', type=str, help='Shape variable', required=True)
    args = parser.parse_args()
    shape =args.shape
    print('Shape: ' + shape)

    # generate the ground truth object
    if shape == 'rect1':
        v_gt = gen_rect1()
    elif shape == 'circle':
        v_gt = gen_circle()
    elif shape == 'ellip2':
        v_gt = gen_ellip2()
    elif shape == 'butter':
        v_gt = gen_butter()
    else:
        print("ERROR: invalid shape!")
        exit()
    
    # generate initial object
    n_init = 10
    v = genInit(n_init, 0.1)

    # display
    plt.axis('equal')
    plt.title(shape)
    print("Select contact points")

    # plot ground truth objects
    # plt.scatter(v_gt[:,0], v_gt[:, 1], c = 'k')
    for i in range(0, len(v_gt)- 1):
        plt.plot([v_gt[i,0], v_gt[(i+1), 0]], [v_gt[i, 1], v_gt[(i+1), 1]], c = 'k', linestyle = '--')
    
    plt.plot([v_gt[0,0], v_gt[len(v_gt) - 1, 0]], [v_gt[0, 1], v_gt[len(v_gt) -1, 1]], c = 'k', linestyle = '--')
    plt.plot([0],[0], 'ro')

    # init plot (not needed)
    # plt.scatter(v[:,0], v[:, 1], c = 'b')
    # for i in range(0, n_init):
    #     plt.plot([v[i,0], v[(i+1) % n_init, 0]], [v[i, 1], v[(i+1) % n_init, 1]], c = 'b', linestyle = '--')
        
    # choose noisy contact points
    c, n = getContacts()

    for i in range(0, len(n)):
        plt.arrow(c[i,0], c[i,1], n[i,0], n[i,1], head_width=0.02, head_length=0.02)

    folder = "contacts/"

    if not os.path.exists(folder):
        os.makedirs(folder)
    cn = np.hstack((c,n))
    np.savetxt(folder + 'contacts-' + shape + '-' + getDateTime() + '.txt', cn, delimiter=',')
    # plt.close()
    plt.show(block = True)
