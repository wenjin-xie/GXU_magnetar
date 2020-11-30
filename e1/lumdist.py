#!/usr/bin/env python
# -*- coding: utf-8 -*-

' change redshift to distance(Mpc)'
__author__='Wenjin Xie'




from numpy import array, ndarray, sqrt, sin, sinh, maximum
from scipy.integrate import quad
from math import sqrt as msqrt


def cosmo_param(omega_m=None, omega_lambda=None, omega_k=None, q0=None):
    nk = omega_k is not None
    nlambda = omega_lambda is not None 
    nomega = omega_m is not None
    nq0 = q0 is not None
    # Check which two parameters are defined, and then determine the other two

    if nomega and nlambda:
        if not nk:   
            omega_k = 1 - omega_m - omega_lambda
        if not nq0:   
            q0 = omega_m / 2. - omega_lambda

    if nomega and nk:
        if not nlambda:   
            omega_lambda = 1. - omega_m - omega_k
        if not nq0:   
            q0 = -1 + omega_k + 3 * omega_m / 2

    if nlambda and nk:   
        if not nomega:   
            omega_m = 1. - omega_lambda - omega_k
        if not nq0:   
            q0 = (1 - omega_k - 3. * omega_lambda) / 2.

    if nomega and nq0:   
        if not nk:   
            omega_k = 1 + q0 - 3 * omega_m / 2.
        if not nlambda:   
            omega_lambda = 1. - omega_m - omega_k

    if nlambda and nq0:   
        if not nk:   
            omega_k = 1 - 2 * q0 - 3 * omega_lambda
        if not nomega:   
            omega_m = 1. - omega_lambda - omega_k

    if nk and nq0:   
        if not nomega:   
            omega_m = (1 + q0 - omega_k) * 2 / 3.
        if not nlambda:   
            omega_lambda = 1. - omega_m - omega_k

    #Set default values
    if omega_k is None:
        omega_k = 0       #Default is flat space
    if omega_lambda is None:
        omega_lambda = 0.7
    if omega_m is None:
        omega_m = 1 - omega_lambda
    if q0 is None:
        q0 = (1 - omega_k - 3 * omega_lambda) / 2.

    return omega_m, omega_lambda, omega_k, q0

def ldist(z, q0=None, lambda0=None):

    term1 = (1. + z) ** 2
    term2 = 1. + 2. * (q0 + lambda0) * z
    term3 = z * (2. + z) * lambda0
    denom = (term1 * term2 - term3)
    if denom>0:
        out = 1. / msqrt(denom) # since the function is used with scalar arguments
                                                                                        # I use math.sqrt instead of numpy.sqrt for
                                                # performance reasons
    else:
        out = 0.
    return out


def lumdist(z, h0=None, k=None, lambda0=None, omega_m=None, q0=None, silent=None):
    '''Syntax: result = lumdist(z, H0 = ,k=, Lambda0 = ])
    Returns luminosity distance in Mpc'''

    scal=False
    scalret = lambda x : x[0] if scal else x

    if isinstance(z, list):
        z = array(z)
    elif isinstance(z, ndarray):
        pass
    else:
        scal = True
        z = array([z])
    n = len(z)

    omega_m, lambda0, k, q0 = cosmo_param(omega_m, lambda0, k, q0)

    # Check keywords
    c = 2.99792458e5                  #  speed of light in km/s
    if h0 is None:
        h0 = 70
   # if not silent:
   #     print 'LUMDIST: H0:', h0, ' Omega_m:', omega_m, ' Lambda0', lambda0, ' q0: ', q0, ' k: ', k#, format='(A,I3,A,f5.2,A,f5.2,A,f5.2,A,F5.2)'

    # For the case of Lambda = 0, we use the closed form from equation 5.238 of
    # Astrophysical Formulae (Lang 1998).   This avoids terms that almost cancel
    # at small q0*z better than the more familiar Mattig formula.
    #
    if lambda0 == 0:   
        denom = sqrt(1 + 2 * q0 * z) + 1 + q0 * z
        dlum = (c * z / h0) * (1 + z * (1 - q0) / denom)
        return scalret(dlum)*3.086*10**24

        # For non-zero lambda
    else:   
        dlum = z * 0.0
        for i in range(n):
            if z[i] <= 0.0:   
                dlum[i] = 0.0
            else:   
                lz = quad(ldist, 0, z[i], args=(q0, lambda0))
                dlum[i] = lz[0]

        if k > 0:   
            dlum = sinh(sqrt(k) * dlum) / sqrt(k)
        else:   
            if k < 0:   
                dlum = maximum(sin(sqrt(-k) * dlum) / sqrt(-k), 0)
        return scalret(c * (1 + z) * dlum / h0)

#-----------------------return Mpc

if __name__=='__main__':
    lumdist()*20000