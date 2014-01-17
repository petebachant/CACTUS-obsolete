# -*- coding: utf-8 -*-
"""
Created on Sun Aug 11 09:50:18 2013

@author: pete
"""

from numpy import zeros, ones
import numpy as np

def quatrot(v, theta, nR, origin):
    """Perform rotation of vector(s) v around normal vector nR using the
    quaternion machinery. \
    v: input vector(s) (can be m x 3 for m vectors to rotate) \
    Theta: rotation angle (rad)
    nR: normal vector around which to rotate
    Origin: origin point of rotation
    
    vR: Rotated vector(s) (size m x 3 for m input vectors)"""

    # Force normalize nR
    nR = nR/np.sqrt(np.sum(nR**2))
    
    # Quaternion form of v
    O = origin[ones(size(v,1),1), :] 
    # I am 12 years old and who uses "O" as a variable name?
    vO = v - O
    p = [zeros(np.size(v,1),1),vO]
    
    # Rotation quaternion and conjugate
    q = [cos(theta/2), nR*sin(theta/2)]
    qbar=[q[1], -q[2:4]]
    
    # These should be matrices
    QL = [q[1] -q[2] -q[3] -q[4]
        q(2)  q(1) -q(4)  q(3)
        q(3)  q(4)  q(1) -q(2)
        q(4) -q(3)  q(2)  q(1)]
       
    QbarR = [qbar(1) -qbar(2) -qbar(3) -qbar(4)
           qbar(2)  qbar(1)  qbar(4) -qbar(3)
           qbar(3) -qbar(4)  qbar(1)  qbar(2)
           qbar(4)  qbar(3) -qbar(2)  qbar(1)]
    
    # Rotate p
    pR = p*(QbarR*QL).transpose()
    vR = pR[:,2:4] + O
    
    return vR

