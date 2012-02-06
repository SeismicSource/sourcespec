# -*- coding: utf-8 -*-
# Parameter file fof source_spec.py
# This file has to be valid Python code and can contain scripting


# GENERAL PARAMETERS
#DEBUG=True # True: print debug information
DEBUG=False # True: print debug information
DOPLOTS=True
#DOPLOTS=False
# -------------------


# INITIAL PARAMETERS
#vp=5.5
vs=3.055
rho=2700 #density
rps=0.62 #radiation pattern coefficient
#qs=12000.
# Brune stress drop
#bsd=20e6 #MPa
bsd=5e6 #Mpa
# -------------------


# SPECTRUM PARAMETERS
# S-wave window
pre_s_time   = 1 #sec
s_win_length = 5 #sec

# Band-pass frequencies for accelerometers and velocimeters
# TODO: calculate from sampling rate?
bp_freqmin_acc = 1.0
bp_freqmax_acc = 50.0
bp_freqmin_vel = 0.25
bp_freqmax_vel = 50.0

# Filter cut frequencies for accelerometers and velocimeters
freq1_acc = 1.0
freq2_acc = 30.0
freq1_vel = 0.5
freq2_vel = 50.0
# -------------------


# INVERSION PARAMETERS
# Spectral weighting:
#   weight for f<=f_weight
#   1      for f> f_weight
f_weight = 7. #Hz
weight = 5.
# Initial values for fc (optional)
fc_0 = 8.
# initial value for t_star
t_star_0 = 0.045
# -------------------


# POST-PROCESSING PARAMETERS
# Min and max acceptable corner frequencies
min_corner_freq = 2
max_corner_freq = 50
