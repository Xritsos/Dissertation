"""
Joule Heating Statistical Error Code
Christos Psychalas, Democritus University of Thrace, Xanthi, Greece
Daedalus Science Study
"""

import datetime
import warnings
import time
import pyglow
import PySimpleGUI as Sg
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from netCDF4 import Dataset
from matplotlib import ticker
from cmcrameri import cm
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ######################## CONSTANTS ########################
qe = 1.602176565 * 10 ** (-19)  # electron charge in coulomb
me = 9.10938356 * 10 ** (-31)   # electron mass in kg
mO = 15.9994   # Oxygen atomic mass in g/mol
mO2 = 2 * mO   # Oxygen molecular mass in g/mol
mN = 14.0067   # Nitrogen atomic mass in g/mol
mN2 = 2 * mN   # Nitrogen molecular mass in g/mol
mNO = mN + mO  # Nitric oxide mass in g/mol

boltzmann = 1.380645852 * 10 ** (-16)  # Boltzmann's constant in cm^2 * g * s^(-2) * K^(-1)
nAvogadro = 6.02214086 * 10 ** 23      # Avogadro's constant in mol^(-1)
fb = 1.5  # Burnside Factor, the factor that connects theoretical and practical estimations of O+

# Masses in kg
mkO = mO / (nAvogadro * 1000)
mkO2 = mO2 / (nAvogadro * 1000)
mkN = mN / (nAvogadro * 1000)
mkN2 = mN2 / (nAvogadro * 1000)
mkNO = mNO / (nAvogadro * 1000)

ccm = 10 ** 6  # constant to convert cm^(-3) to m^(-3)

# Altitudes
heights = np.zeros((72, 144, 57), order='F')

# Altitude used in Lat - Alt plot
heights_la = np.zeros(57, order='F')

# Magnetic field
Bx = np.zeros((72, 144, 57), order='F')
By = np.zeros((72, 144, 57), order='F')
Bz = np.zeros((72, 144, 57), order='F')

Bx_noisy = np.zeros((72, 144, 57), order='F')
By_noisy = np.zeros((72, 144, 57), order='F')
Bz_noisy = np.zeros((72, 144, 57), order='F')

# Electric field
Ex = np.zeros((72, 144, 57), order='F')
Ey = np.zeros((72, 144, 57), order='F')
Ez = np.zeros((72, 144, 57), order='F')

Ex_noisy = np.zeros((72, 144, 57), order='F')
Ey_noisy = np.zeros((72, 144, 57), order='F')
Ez_noisy = np.zeros((72, 144, 57), order='F')

# Neutral wind
Unx = np.zeros((72, 144, 57), order='F')
Uny = np.zeros((72, 144, 57), order='F')
Unz = np.zeros((72, 144, 57), order='F')

Unx_noisy = np.zeros((72, 144, 57), order='F')
Uny_noisy = np.zeros((72, 144, 57), order='F')
Unz_noisy = np.zeros((72, 144, 57), order='F')

# Densities
NO = np.zeros((72, 144, 57), order='F')
NO2 = np.zeros((72, 144, 57), order='F')
NN2 = np.zeros((72, 144, 57), order='F')
NOp = np.zeros((72, 144, 57), order='F')
NO2p = np.zeros((72, 144, 57), order='F')
NNOp = np.zeros((72, 144, 57), order='F')
Ne = np.zeros((72, 144, 57), order='F')

NO_noisy = np.zeros((72, 144, 57), order='F')
NO2_noisy = np.zeros((72, 144, 57), order='F')
NN2_noisy = np.zeros((72, 144, 57), order='F')
NOp_noisy = np.zeros((72, 144, 57), order='F')
NO2p_noisy = np.zeros((72, 144, 57), order='F')
NNOp_noisy = np.zeros((72, 144, 57), order='F')
Ne_noisy = np.zeros((72, 144, 57), order='F')

# Temperatures
Te = np.zeros((72, 144, 57), order='F')
Ti = np.zeros((72, 144, 57), order='F')
Tn = np.zeros((72, 144, 57), order='F')

Te_noisy = np.zeros((72, 144, 57), order='F')
Ti_noisy = np.zeros((72, 144, 57), order='F')
Tn_noisy = np.zeros((72, 144, 57), order='F')

# Collision frequencies
nu_Op_sum = np.zeros((72, 144, 57), order='F')
nu_O2p_sum = np.zeros((72, 144, 57), order='F')
nu_NOp_sum = np.zeros((72, 144, 57), order='F')
nu_e_sum = np.zeros((72, 144, 57), order='F')

nu_Op_sum_noisy = np.zeros((72, 144, 57), order='F')
nu_O2p_sum_noisy = np.zeros((72, 144, 57), order='F')
nu_NOp_sum_noisy = np.zeros((72, 144, 57), order='F')
nu_e_sum_noisy = np.zeros((72, 144, 57), order='F')

# Conductivities
pedersen_con = np.zeros((72, 144, 57), order='F')
hall_con = np.zeros((72, 144, 57), order='F')
parallel_con = np.zeros((72, 144, 57), order='F')

pedersen_con_noisy = np.zeros((72, 144, 57), order='F')
hall_con_noisy = np.zeros((72, 144, 57), order='F')
parallel_con_noisy = np.zeros((72, 144, 57), order='F')

# Ion velocities perpendicular to magnetic field
Vi_vertx = np.zeros((72, 144, 57), order='F')
Vi_verty = np.zeros((72, 144, 57), order='F')
Vi_vertz = np.zeros((72, 144, 57), order='F')

Vi_vertx_noisy = np.zeros((72, 144, 57), order='F')
Vi_verty_noisy = np.zeros((72, 144, 57), order='F')
Vi_vertz_noisy = np.zeros((72, 144, 57), order='F')

# Heating rates
Joule_Heating = np.zeros((72, 144, 57), order='F')
Frictional_Heating = np.zeros((72, 144, 57), order='F')
Ohmic_Heating = np.zeros((72, 144, 57), order='F')

Joule_Heating_noisy = np.zeros((72, 144, 57), order='F')
Frictional_Heating_noisy = np.zeros((72, 144, 57), order='F')
Ohmic_Heating_noisy = np.zeros((72, 144, 57), order='F')

# Cross sections
C_Op = np.zeros((72, 144, 57), order='F')
C_O2p = np.zeros((72, 144, 57), order='F')
C_NOp = np.zeros((72, 144, 57), order='F')
C_ion = np.zeros((72, 144, 57), order='F')

C_Op_noisy = np.zeros((72, 144, 57), order='F')
C_O2p_noisy = np.zeros((72, 144, 57), order='F')
C_NOp_noisy = np.zeros((72, 144, 57), order='F')
C_ion_noisy = np.zeros((72, 144, 57), order='F')

# Perpendicular currents
J_pedersen = np.zeros((72, 144, 57), order='F')
J_hall = np.zeros((72, 144, 57), order='F')
J_ohmic = np.zeros((72, 144, 57), order='F')
J_dens = np.zeros((72, 144, 57), order='F')

J_pedersen_noisy = np.zeros((72, 144, 57), order='F')
J_hall_noisy = np.zeros((72, 144, 57), order='F')
J_ohmic_noisy = np.zeros((72, 144, 57), order='F')
J_dens_noisy = np.zeros((72, 144, 57), order='F')


# Convert east north up to earth centered earth fixed
# lat,lon from TIE-GCM correspond to perfect sphere(spherical latitude and longitude)
def enu_ecef(lat_phi, lon_lmd, Fe, Fn, Fup):
    fac = np.pi / 180
    lat_phi = lat_phi * fac
    lon_lmd = lon_lmd * fac

    north_temp_unit = [-np.cos(lon_lmd) * np.sin(lat_phi), - np.sin(lon_lmd) * np.sin(lat_phi), np.cos(lat_phi)]
    east_temp_unit = [-np.sin(lon_lmd), np.cos(lon_lmd), 0]
    up_temp_unit = [np.cos(lon_lmd) * np.cos(lat_phi), np.sin(lon_lmd) * np.cos(lat_phi), np.sin(lat_phi)]

    Fnorth_vector = [Fn * north_temp_unit[0], Fn * north_temp_unit[1], Fn * north_temp_unit[2]]
    Feast_vector = [Fe * east_temp_unit[0], Fe * east_temp_unit[1], Fe * east_temp_unit[2]]
    Fup_vector = [Fup * up_temp_unit[0], Fup * up_temp_unit[1], Fup * up_temp_unit[2]]

    Fx = Fnorth_vector[0] + Feast_vector[0] + Fup_vector[0]
    Fy = Fnorth_vector[1] + Feast_vector[1] + Fup_vector[1]
    Fz = Fnorth_vector[2] + Feast_vector[2] + Fup_vector[2]

    return Fx, Fy, Fz


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ GET INPUT FROM TIE-GCM AND I-GRF $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def models_input(file_name, timer, lat_value=-1, lon_value=-1, pressure_level=-1):
    warnings.filterwarnings("ignore")
    start_time = time.time()

    print('Getting data from models........')

    global glat_in
    global glon_in
    global glev_in
    global title
    global map_time

    # Get TIE-GCM file from gui
    tiegcm_file = "TIEGCM_FILES_LIFETIME_2015/" + file_name

    # Open TIE-GCM file
    tiegcm = Dataset(tiegcm_file)

    # Input parameters from TIE-GCM
    zgmid_in = tiegcm.variables['ZGMID'][:]  # geometric height at midpoints
    glat_in = tiegcm.variables['lat'][:]     # geographic latitude in deg(earth as perfect sphere)
    glon_in = tiegcm.variables['lon'][:]     # geographic longitude in deg(earth as perfect sphere)
    glev_in = tiegcm.variables['lev'][:]     # midpoint levels
    time_in = tiegcm.variables['time'][:]    # minutes since 2015-1-1 0:0:0
    O_in = tiegcm.variables['O_CM3'][:]      # atomic oxygen density (neutral) in cm^(-3)
    O2_in = tiegcm.variables['O2_CM3'][:]    # molecular oxygen density (neutral) in cm^(-3)
    N2_in = tiegcm.variables['N2_CM3'][:]    # molecular nitrogen density (neutral) in cm^(-3)
    Op_in = tiegcm.variables['OP'][:]        # atomic oxygen density (ion) in cm^(-3)
    O2p_in = tiegcm.variables['O2P'][:]      # molecular oxygen density (ion) in cm^(-3)
    NOp_in = tiegcm.variables['NOP_LAM'][:]  # nitric oxide density (ion) in cm^(-3)
    Te_in = tiegcm.variables['TE'][:]        # electron temperature in kelvin
    Tn_in = tiegcm.variables['TN'][:]        # neutral temperature in kelvin
    Ti_in = tiegcm.variables['TI'][:]        # ion temperature in kelvin
    Un_north_in = tiegcm.variables['VN_si'][:]  # neutral meridional wind (+north) in m/sec
    Un_east_in = tiegcm.variables['UN_si'][:]   # neutral zonal wind (+east) in m/sec
    Un_up_in = tiegcm.variables['WN_si'][:]     # neutral vertical wind (+up) in m/sec

    Ee_in = tiegcm.variables['EEX_si'][:]  # electric field (+east) in V/m
    En_in = tiegcm.variables['EEY_si'][:]  # electric field (+north) in V/m
    Eu_in = tiegcm.variables['EEZ_si'][:]  # electric field (+up) in V/m

    # Close TIE-GCM file
    tiegcm.close()

    timeg = time_in[timer]
    time_IGRF = datetime.datetime(2015, 1, 1, 0, 0, 0)         # first time step of the TIE-GCM run
    real_time = time_IGRF + datetime.timedelta(minutes=timeg)  # time user defined with time-step

    # Distinguish Map from Vertical profile
    lev_range = 0
    lat_range = 0
    lon_range = 0

    lev_start = 0
    lat_start = 0
    lon_start = 0

    # Lat - Alt map profile
    if lat_value == -1 and pressure_level == -1:
        lev_range = len(glev_in) - 1
        lat_range = len(glat_in)
        lon_start = lon_value
        lon_range = lon_start + 1
        title = ' Lon: ' + str(glon_in[lon_value]) + ' ,' + str(real_time)
        map_time = real_time

    # Lat - Lon map profile
    if lat_value == -1 and lon_value == -1:
        lev_start = pressure_level
        lev_range = lev_start + 1
        lat_range = len(glat_in)
        lon_range = len(glon_in)
        title = ' Pressure level:' + str(glev_in[pressure_level]) + ' , ' + str(real_time)
        map_time = real_time

    # Vertical profile
    if lat_value != -1 and lon_value != -1:
        lat_start = lat_value
        lon_start = lon_value
        lat_range = lat_start + 1
        lon_range = lon_start + 1
        lev_range = len(glev_in) - 1
        title = ' Lat: ' + str(glat_in[lat_value]) + ' ,' + ' Lon: ' + str(glon_in[lon_value]) + ' , ' + str(real_time)

    for lev in range(lev_start, lev_range):
        temp_height = 0
        for lat in range(lat_start, lat_range):
            for lon in range(lon_start, lon_range):

                heights[lat, lon, lev] = zgmid_in[timer, lev, lat, lon] / 1e5  # altitude in km

                # Average heights for Lat - Alt map
                temp_height = temp_height + heights[lat, lon, lev]
                heights_la[lev] = round(temp_height / lat_range)

                # Run I-GRF model using pyglow to get magnetic field
                # Create pyglow point
                pt = pyglow.Point(real_time, glat_in[lat], glon_in[lon], heights[lat, lon, lev], user_ind=False)

                # Run I-GRF
                pt.run_igrf()  # default version is I-GRF-12

                # Magnetic field in ENU
                Be = pt.Bx  # in tesla
                Bn = pt.By  # in tesla
                Bu = pt.Bz  # in tesla

                # Magnetic field from ENU to ECEF (in tesla)
                Bx[lat, lon, lev], By[lat, lon, lev], Bz[lat, lon, lev] = enu_ecef(glat_in[lat], glon_in[lon], Be, Bn, Bu)

                # Electric field from ENU to ECEF (in V/m)
                Ee_temp = Ee_in[timer, lev, lat, lon]
                En_temp = En_in[timer, lev, lat, lon]
                Eu_temp = Eu_in[timer, lev, lat, lon]
                Ex[lat, lon, lev], Ey[lat, lon, lev], Ez[lat, lon, lev] = enu_ecef(glat_in[lat], glon_in[lon], Ee_temp, En_temp, Eu_temp)

                # Neutral wind from ENU to ECEF (in m/sec)
                Un_e_temp = Un_east_in[timer, lev, lat, lon]
                Un_n_temp = Un_north_in[timer, lev, lat, lon]
                Un_u_temp = Un_up_in[timer, lev, lat, lon]
                Unx[lat, lon, lev], Uny[lat, lon, lev], Unz[lat, lon, lev] = enu_ecef(glat_in[lat], glon_in[lon], Un_e_temp, Un_n_temp, Un_u_temp)

                # Assign densities (in cm^(-3))
                NO[lat, lon, lev] = O_in[timer, lev, lat, lon]
                NO2[lat, lon, lev] = O2_in[timer, lev, lat, lon]
                NN2[lat, lon, lev] = N2_in[timer, lev, lat, lon]
                NOp[lat, lon, lev] = Op_in[timer, lev, lat, lon]
                NO2p[lat, lon, lev] = O2p_in[timer, lev, lat, lon]
                NNOp[lat, lon, lev] = NOp_in[timer, lev, lat, lon]

                # Force charge neutrality ignoring other minor ion densities
                Ne[lat, lon, lev] = NOp[lat, lon, lev] + NO2p[lat, lon, lev] + NNOp[lat, lon, lev]

                # Assign temperatures (in kelvin)
                Te[lat, lon, lev] = Te_in[timer, lev, lat, lon]
                Ti[lat, lon, lev] = Ti_in[timer, lev, lat, lon]
                Tn[lat, lon, lev] = Tn_in[timer, lev, lat, lon]

                # Measurements for each ion cannot be made separately so we take the average
                # no coordinate system assigned to satellite, assumed as point measurement
                # perpendicular ion velocity corresponds to measured ion velocity
                # in contrast to neutral wind, which we get as model input

                nu_Op_N2 = 6.82 * NN2[lat, lon, lev] * 10 ** (-10)
                nu_Op_O2 = 6.64 * NO2[lat, lon, lev] * 10 ** (-10)

                Tr = (Ti[lat, lon, lev] + Tn[lat, lon, lev]) / 2  # in kelvin
                nu_Op_O = fb * (3.67 * NO[lat, lon, lev] * 10 ** (-11) * Tr ** (1 / 2) * (1 - 0.064 * np.log10(Tr)) ** 2)

                nu_Op_sum_temp = nu_Op_N2 + nu_Op_O2 + nu_Op_O

                # nu-O2p = nu(O2p-N2) + nu(O2p-O) + nu(O2p-O2) (in Hz)
                # densities in cm^(-3)
                nu_O2p_N2 = 4.13 * NN2[lat, lon, lev] * 10 ** (-10)
                nu_O2p_O = 2.31 * NO[lat, lon, lev] * 10 ** (-10)
                nu_O2p_O2 = 2.59 * NO2[lat, lon, lev] * 10 ** (-11) * Tr ** (1 / 2) * (1 - 0.073 * np.log10(Tr)) ** 2

                nu_O2p_sum_temp = nu_O2p_N2 + nu_O2p_O + nu_O2p_O2

                # nu-NOp = nu(NOp-N2) + nu(NOp-O) + nu(NOp-O2) (in Hz)
                # densities in cm^(-3)
                nu_NOp_N2 = 4.34 * NN2[lat, lon, lev] * 10 ** (-10)
                nu_NOp_O = 2.44 * NO[lat, lon, lev] * 10 ** (-10)
                nu_NOp_O2 = 4.27 * NO2[lat, lon, lev] * 10 ** (-10)

                nu_NOp_sum_temp = nu_NOp_N2 + nu_NOp_O + nu_NOp_O2

                # ################ GYRO-FREQUENCIES(OMEGAS) ################
                # Magnetic field vector (in tesla)
                B = [Bx[lat, lon, lev], By[lat, lon, lev], Bz[lat, lon, lev]]
                # Electric field vector(in volt/meter)
                E = [Ex[lat, lon, lev], Ey[lat, lon, lev], Ez[lat, lon, lev]]
                # Magnetic field magnitude (in tesla)
                Bnorm = np.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)
                # Magnetic field unit vector
                bunit = [B[0] / Bnorm, B[1] / Bnorm, B[2] / Bnorm]

                # qe(in coulomb), mk(masses in kg), omegas(in Hz)
                omega_Op = (qe * Bnorm) / mkO
                omega_O2p = (qe * Bnorm) / mkO2
                omega_NOp = (qe * Bnorm) / mkNO

                # Electric field perpendicular to magnetic field
                # Evert = E cross bunit
                Evertx = E[1] * bunit[2] - E[2] * bunit[1]
                Everty = E[2] * bunit[0] - E[0] * bunit[2]
                Evertz = E[0] * bunit[1] - E[1] * bunit[0]

                # E vertical vector
                Evert = [Evertx, Everty, Evertz]

                # Neutral wind vector(in meter/sec)
                Un = [Unx[lat, lon, lev], Uny[lat, lon, lev], Unz[lat, lon, lev]]

                # Neutral wind perpendicular to magnetic field
                # Unvert = Un cross bunit
                Un_vertx = Un[1] * bunit[2] - Un[2] * bunit[1]
                Un_verty = Un[2] * bunit[0] - Un[0] * bunit[2]
                Un_vertz = Un[0] * bunit[1] - Un[1] * bunit[0]

                # Un perpendicular to magnetic field vector
                Un_vert = [Un_vertx, Un_verty, Un_vertz]

                # Unvert cross B vector
                UnvertXBx = Un_vert[1] * B[2] - Un_vert[2] * B[1]
                UnvertXBy = Un_vert[2] * B[0] - Un_vert[0] * B[2]
                UnvertXBz = Un_vert[0] * B[1] - Un_vert[1] * B[0]
                UnvXB = [UnvertXBx, UnvertXBy, UnvertXBz]

                # Estar: perpendicular electric field in the neutral frame Estar = Evert + Unvert cross B
                # vector addition
                Estar_x = Evert[0] + UnvXB[0]
                Estar_y = Evert[1] + UnvXB[1]
                Estar_z = Evert[2] + UnvXB[2]

                Estar = [Estar_x, Estar_y, Estar_z]

                # Estar cross bunit
                x = Estar[1] * bunit[2] - Estar[2] * bunit[1]
                y = Estar[2] * bunit[0] - Estar[0] * bunit[2]
                z = Estar[0] * bunit[1] - Estar[1] * bunit[0]

                EstarXbunit = [x, y, z]

                # Ion velocities vectors(in meter/sec) (in neutral frame, star) perpendicular to magnetic field
                # extracted from ion momentum equation
                # ############# O+ #############
                Vi_Op_star_x = (nu_Op_sum_temp * omega_Op * Estar[0] + omega_Op ** 2 * EstarXbunit[0]) / \
                               (Bnorm * (nu_Op_sum_temp ** 2 + omega_Op ** 2))
                Vi_Op_star_y = (nu_Op_sum_temp * omega_Op * Estar[1] + omega_Op ** 2 * EstarXbunit[1]) / \
                               (Bnorm * (nu_Op_sum_temp ** 2 + omega_Op ** 2))
                Vi_Op_star_z = (nu_Op_sum_temp * omega_Op * Estar[2] + omega_Op ** 2 * EstarXbunit[2]) / \
                               (Bnorm * (nu_Op_sum_temp ** 2 + omega_Op ** 2))
                # ############# O2+ #############
                Vi_O2p_star_x = (nu_O2p_sum_temp * omega_O2p * Estar[0] + omega_O2p ** 2 * EstarXbunit[0]) / \
                                (Bnorm * (nu_O2p_sum_temp ** 2 + omega_O2p ** 2))
                Vi_O2p_star_y = (nu_O2p_sum_temp * omega_O2p * Estar[1] + omega_O2p ** 2 * EstarXbunit[1]) / \
                                (Bnorm * (nu_O2p_sum_temp ** 2 + omega_O2p ** 2))
                Vi_O2p_star_z = (nu_O2p_sum_temp * omega_O2p * Estar[2] + omega_O2p ** 2 * EstarXbunit[2]) / \
                                (Bnorm * (nu_O2p_sum_temp ** 2 + omega_O2p ** 2))
                # ############ NO+ ############
                Vi_NOp_star_x = (nu_NOp_sum_temp * omega_NOp * Estar[0] + omega_NOp ** 2 * EstarXbunit[0]) / \
                                (Bnorm * (nu_NOp_sum_temp ** 2 + omega_NOp ** 2))
                Vi_NOp_star_y = (nu_NOp_sum_temp * omega_NOp * Estar[1] + omega_NOp ** 2 * EstarXbunit[1]) / \
                                (Bnorm * (nu_NOp_sum_temp ** 2 + omega_NOp ** 2))
                Vi_NOp_star_z = (nu_NOp_sum_temp * omega_NOp * Estar[2] + omega_NOp ** 2 * EstarXbunit[2]) / \
                                (Bnorm * (nu_NOp_sum_temp ** 2 + omega_NOp ** 2))

                # Changing frame from neutral το ECEF frame(not star)
                Vi_Op_x = Vi_Op_star_x + Un_vert[0]
                Vi_O2p_x = Vi_O2p_star_x + Un_vert[0]
                Vi_NOp_x = Vi_NOp_star_x + Un_vert[0]

                Vi_Op_y = Vi_Op_star_y + Un_vert[1]
                Vi_O2p_y = Vi_O2p_star_y + Un_vert[1]
                Vi_NOp_y = Vi_NOp_star_y + Un_vert[1]

                Vi_Op_z = Vi_Op_star_z + Un_vert[2]
                Vi_O2p_z = Vi_O2p_star_z + Un_vert[2]
                Vi_NOp_z = Vi_NOp_star_z + Un_vert[2]

                # Ion velocity perpendicular to magnetic field
                Vi_vertx[lat, lon, lev] = (Vi_Op_x + Vi_O2p_x + Vi_NOp_x) / 3
                Vi_verty[lat, lon, lev] = (Vi_Op_y + Vi_O2p_y + Vi_NOp_y) / 3
                Vi_vertz[lat, lon, lev] = (Vi_Op_z + Vi_O2p_z + Vi_NOp_z) / 3

    print('Data imported in: ', str(time.time() - start_time), ' sec!')
    print(' ')
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END OF INPUT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ AWGN FUNCTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def calculate_noise_2(lat_value=-1, lon_value=-1, pressure_level=-1):
    start_time = time.time()
    lat = lat_value
    lon = lon_value
    lev = pressure_level
    
    snr = 45

    # Vertical profile case
    if lat != -1 and lon != -1:
        rms_data = np.sqrt(np.mean(Bx[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Bx[lat, lon, :].shape)
        Bx_noisy[lat, lon, :] = Bx[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(By[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, By[lat, lon, :].shape)
        By_noisy[lat, lon, :] = By[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(Bz[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Bz[lat, lon, :].shape)
        Bz_noisy[lat, lon, :] = Bz[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(Ex[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Ex[lat, lon, :].shape)
        Ex_noisy[lat, lon, :] = Ex[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(Ey[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Ey[lat, lon, :].shape)
        Ey_noisy[lat, lon, :] = Ey[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(Ez[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Ez[lat, lon, :].shape)
        Ez_noisy[lat, lon, :] = Ez[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(Unx[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Unx[lat, lon, :].shape)
        Unx_noisy[lat, lon, :] = Unx[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(Uny[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Uny[lat, lon, :].shape)
        Uny_noisy[lat, lon, :] = Uny[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(Unz[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Unz[lat, lon, :].shape)
        Unz_noisy[lat, lon, :] = Unz[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(Vi_vertx[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Vi_vertx[lat, lon, :].shape)
        Vi_vertx_noisy[lat, lon, :] = Vi_vertx[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(Vi_verty[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Vi_verty[lat, lon, :].shape)
        Vi_verty_noisy[lat, lon, :] = Vi_verty[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(Vi_vertz[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Vi_vertz[lat, lon, :].shape)
        Vi_vertz_noisy[lat, lon, :] = Vi_vertz[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(NOp[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NOp[lat, lon, :].shape)
        NOp_noisy[lat, lon, :] = NOp[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(NO2p[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NO2p[lat, lon, :].shape)
        NO2p_noisy[lat, lon, :] = NO2p[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(NNOp[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NNOp[lat, lon, :].shape)
        NNOp_noisy[lat, lon, :] = NNOp[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(NO[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NO[lat, lon, :].shape)
        NO_noisy[lat, lon, :] = NO[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(NO2[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NO2[lat, lon, :].shape)
        NO2_noisy[lat, lon, :] = NO2[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(NN2[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NN2[lat, lon, :].shape)
        NN2_noisy[lat, lon, :] = NN2[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(Ne[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Ne[lat, lon, :].shape)
        Ne_noisy[lat, lon, :] = Ne[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(Ti[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Ti[lat, lon, :].shape)
        Ti_noisy[lat, lon, :] = Ti[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(Tn[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Tn[lat, lon, :].shape)
        Tn_noisy[lat, lon, :] = Tn[lat, lon, :] + noise
    
        rms_data = np.sqrt(np.mean(Te[lat, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Te[lat, lon, :].shape)
        Te_noisy[lat, lon, :] = Te[lat, lon, :] + noise
        
    # Latitude - Longitude map profile case
    if lat == -1 and lon == -1:
        rms_data = np.sqrt(np.mean(Bx[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Bx[:, :, lev].shape)
        Bx_noisy[:, :, lev] = Bx[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(By[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, By[:, :, lev].shape)
        By_noisy[:, :, lev] = By[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(Bz[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Bz[:, :, lev].shape)
        Bz_noisy[:, :, lev] = Bz[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(Ex[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Ex[:, :, lev].shape)
        Ex_noisy[:, :, lev] = Ex[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(Ey[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Ey[:, :, lev].shape)
        Ey_noisy[:, :, lev] = Ey[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(Ez[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Ez[:, :, lev].shape)
        Ez_noisy[:, :, lev] = Ez[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(Unx[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Unx[:, :, lev].shape)
        Unx_noisy[:, :, lev] = Unx[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(Uny[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Uny[:, :, lev].shape)
        Uny_noisy[:, :, lev] = Uny[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(Unz[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Unz[:, :, lev].shape)
        Unz_noisy[:, :, lev] = Unz[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(Vi_vertx[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Vi_vertx[:, :, lev].shape)
        Vi_vertx_noisy[:, :, lev] = Vi_vertx[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(Vi_verty[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Vi_verty[:, :, lev].shape)
        Vi_verty_noisy[:, :, lev] = Vi_verty[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(Vi_vertz[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Vi_vertz[:, :, lev].shape)
        Vi_vertz_noisy[:, :, lev] = Vi_vertz[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(NOp[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NOp[:, :, lev].shape)
        NOp_noisy[:, :, lev] = NOp[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(NO2p[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NO2p[:, :, lev].shape)
        NO2p_noisy[:, :, lev] = NO2p[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(NNOp[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NNOp[:, :, lev].shape)
        NNOp_noisy[:, :, lev] = NNOp[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(NO[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NO[:, :, lev].shape)
        NO_noisy[:, :, lev] = NO[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(NO2[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NO2[:, :, lev].shape)
        NO2_noisy[:, :, lev] = NO2[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(NN2[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NN2[:, :, lev].shape)
        NN2_noisy[:, :, lev] = NN2[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(Ne[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Ne[:, :, lev].shape)
        Ne_noisy[:, :, lev] = Ne[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(Ti[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Ti[:, :, lev].shape)
        Ti_noisy[:, :, lev] = Ti[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(Tn[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Tn[:, :, lev].shape)
        Tn_noisy[:, :, lev] = Tn[:, :, lev] + noise
    
        rms_data = np.sqrt(np.mean(Te[:, :, lev] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Te[:, :, lev].shape)
        Te_noisy[:, :, lev] = Te[:, :, lev] + noise
        
        # Latitude Altitude map profile case
    if lat == -1 and lev == -1:
        rms_data = np.sqrt(np.mean(Bx[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Bx[:, lon, :].shape)
        Bx_noisy[:, lon, :] = Bx[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(By[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, By[:, lon, :].shape)
        By_noisy[:, lon, :] = By[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(Bz[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Bz[:, lon, :].shape)
        Bz_noisy[:, lon, :] = Bz[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(Ex[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Ex[:, lon, :].shape)
        Ex_noisy[:, lon, :] = Ex[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(Ey[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Ey[:, lon, :].shape)
        Ey_noisy[:, lon, :] = Ey[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(Ez[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Ez[:, lon, :].shape)
        Ez_noisy[:, lon, :] = Ez[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(Unx[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Unx[:, lon, :].shape)
        Unx_noisy[:, lon, :] = Unx[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(Uny[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Uny[:, lon, :].shape)
        Uny_noisy[:, lon, :] = Uny[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(Unz[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Unz[:, lon, :].shape)
        Unz_noisy[:, lon, :] = Unz[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(Vi_vertx[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Vi_vertx[:, lon, :].shape)
        Vi_vertx_noisy[:, lon, :] = Vi_vertx[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(Vi_verty[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Vi_verty[:, lon, :].shape)
        Vi_verty_noisy[:, lon, :] = Vi_verty[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(Vi_vertz[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Vi_vertz[:, lon, :].shape)
        Vi_vertz_noisy[:, lon, :] = Vi_vertz[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(NOp[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NOp[:, lon, :].shape)
        NOp_noisy[:, lon, :] = NOp[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(NO2p[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NO2p[:, lon, :].shape)
        NO2p_noisy[:, lon, :] = NO2p[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(NNOp[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NNOp[:, lon, :].shape)
        NNOp_noisy[:, lon, :] = NNOp[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(NO[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NO[:, lon, :].shape)
        NO_noisy[:, lon, :] = NO[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(NO2[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NO2[:, lon, :].shape)
        NO2_noisy[:, lon, :] = NO2[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(NN2[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, NN2[:, lon, :].shape)
        NN2_noisy[:, lon, :] = NN2[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(Ne[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Ne[:, lon, :].shape)
        Ne_noisy[:, lon, :] = Ne[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(Ti[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Ti[:, lon, :].shape)
        Ti_noisy[:, lon, :] = Ti[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(Tn[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Tn[:, lon, :].shape)
        Tn_noisy[:, lon, :] = Tn[:, lon, :] + noise

        rms_data = np.sqrt(np.mean(Te[:, lon, :] ** 2))
        rms_noise = rms_data / np.sqrt(10 ** (snr / 10))
        noise = np.random.normal(0, rms_noise, Te[:, lon, :].shape)
        Te_noisy[:, lon, :] = Te[:, lon, :] + noise

    print('Calculated Noise in: ', time.time() - start_time, ' sec !')
    print(" ")
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END OF AWGN FUNCTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ CREATE RANDOM NOISE $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def calculate_noise(lat_value=-1, lon_value=-1, pressure_level=-1):
    start_time = time.time()
    print('Calculating Noise.....')
    print(' ')

    # Distinguish Map from Vertical profile
    lev_range = 0
    lat_range = 0
    lon_range = 0

    lev_start = 0
    lat_start = 0
    lon_start = 0

    # Lat - Alt map profile
    if lat_value == -1 and pressure_level == -1:
        lev_range = len(glev_in) - 1
        lat_range = len(glat_in)
        lon_start = lon_value
        lon_range = lon_start + 1

    # Lat - Lon map profile
    if lat_value == -1 and lon_value == -1:
        lev_start = pressure_level
        lev_range = lev_start + 1
        lat_range = len(glat_in)
        lon_range = len(glon_in)

    # Vertical profile
    if lat_value != -1 and lon_value != -1:
        lat_start = lat_value
        lon_start = lon_value
        lat_range = lat_start + 1
        lon_range = lon_start + 1
        lev_range = len(glev_in) - 1

    # Assigning Accuracies
    # Magnetic field (in tesla)
    Bx_ac = 5 * 10 ** (-9)
    By_ac = 5 * 10 ** (-9)
    Bz_ac = 5 * 10 ** (-9)

    # Electric field (in V/m)
    Ex_ac = 2 * 10 ** (-3)
    Ey_ac = 2 * 10 ** (-3)
    Ez_ac = 2 * 10 ** (-3)

    # Neutral densities
    NO_ac = 20 / 100
    NO2_ac = 20 / 100
    NN2_ac = 20 / 100

    # Ion densities
    NOp_ac = 10 / 100
    NO2p_ac = 10 / 100
    NNOp_ac = 10 / 100
    Ne_ac = 10 / 100

    # Temperatures
    Ti_ac = 10 / 100
    Te_ac = 10 / 100
    Tn_ac = 20 / 100

    # Neutral wind (in m/s)
    Unx_ac = 10
    Uny_ac = 20
    Unz_ac = 10

    # Ion drift (in m/s)
    Vix_ac = 100
    Viy_ac = 100
    Viz_ac = 100

    for lev in range(lev_start, lev_range):
        for lat in range(lat_start, lat_range):
            for lon in range(lon_start, lon_range):
                # Magnetic field
                noise = np.random.normal(0, Bx_ac)
                Bx_noisy[lat, lon, lev] = Bx[lat, lon, lev] + noise
                noise = np.random.normal(0, By_ac)
                By_noisy[lat, lon, lev] = By[lat, lon, lev] + noise
                noise = np.random.normal(0, Bz_ac)
                Bz_noisy[lat, lon, lev] = Bz[lat, lon, lev] + noise

                # Electric field
                noise = np.random.normal(0, Ex_ac)
                Ex_noisy[lat, lon, lev] = Ex[lat, lon, lev] + noise
                noise = np.random.normal(0, Ey_ac)
                Ey_noisy[lat, lon, lev] = Ey[lat, lon, lev] + noise
                noise = np.random.normal(0, Ez_ac)
                Ez_noisy[lat, lon, lev] = Ez[lat, lon, lev] + noise

                # Neutral densities
                noise = np.random.normal(0, NO_ac * NO[lat, lon, lev])
                NO_noisy[lat, lon, lev] = NO[lat, lon, lev] + noise
                noise = np.random.normal(0, NO2_ac * NO2[lat, lon, lev])
                NO2_noisy[lat, lon, lev] = NO2[lat, lon, lev] + noise
                noise = np.random.normal(0, NN2_ac * NN2[lat, lon, lev])
                NN2_noisy[lat, lon, lev] = NN2[lat, lon, lev] + noise

                # Ion - Electron densities
                noise = np.random.normal(0, NOp_ac * NOp[lat, lon, lev])
                NOp_noisy[lat, lon, lev] = NOp[lat, lon, lev] + noise
                noise = np.random.normal(0, NO2p_ac * NO2p[lat, lon, lev])
                NO2p_noisy[lat, lon, lev] = NO2p[lat, lon, lev] + noise
                noise = np.random.normal(0, NNOp_ac * NNOp[lat, lon, lev])
                NNOp_noisy[lat, lon, lev] = NNOp[lat, lon, lev] + noise
                noise = np.random.normal(0, Ne_ac * Ne[lat, lon, lev])
                Ne_noisy[lat, lon, lev] = Ne[lat, lon, lev] + noise

                # Temperatures
                noise = np.random.normal(0, Ti_ac * Ti[lat, lon, lev])
                Ti_noisy[lat, lon, lev] = Ti[lat, lon, lev] + noise
                noise = np.random.normal(0, Tn_ac * Tn[lat, lon, lev])
                Tn_noisy[lat, lon, lev] = Tn[lat, lon, lev] + noise
                noise = np.random.normal(0, Te_ac * Te[lat, lon, lev])
                Te_noisy[lat, lon, lev] = Te[lat, lon, lev] + noise

                # Neutral wind
                noise = np.random.normal(0, Unx_ac)
                Unx_noisy[lat, lon, lev] = Unx[lat, lon, lev] + noise
                noise = np.random.normal(0, Uny_ac)
                Uny_noisy[lat, lon, lev] = Uny[lat, lon, lev] + noise
                noise = np.random.normal(0, Unz_ac)
                Unz_noisy[lat, lon, lev] = Unz[lat, lon, lev] + noise

                # Ion drift
                noise = np.random.normal(0, Vix_ac)
                Vi_vertx_noisy[lat, lon, lev] = Vi_vertx[lat, lon, lev] + noise
                noise = np.random.normal(0, Viy_ac)
                Vi_verty_noisy[lat, lon, lev] = Vi_verty[lat, lon, lev] + noise
                noise = np.random.normal(0, Viz_ac)
                Vi_vertz_noisy[lat, lon, lev] = Vi_vertz[lat, lon, lev] + noise

    print('Calculated Noise in: ', time.time() - start_time, ' sec !')
    print(" ")
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Noise Calculation End $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Calculate Products $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def calculate_products(Bx_in, By_in, Bz_in, Ex_in, Ey_in, Ez_in, Unx_in, Uny_in, Unz_in, Vix_in, Viy_in, Viz_in, NO_in, NO2_in, NN2_in, NOp_in,
                       NO2p_in, NNOp_in, Ne_in, Ti_in, Te_in, Tn_in, lat_value=-1, lon_value=-1, pressure_level=-1):

    start_time = time.time()
    print('Calculating Products.....')
    print(' ')

    # Create temporary variables for output
    # Collision frequencies
    nu_Op_sum_temp = np.zeros((72, 144, 57), order='F')
    nu_O2p_sum_temp = np.zeros((72, 144, 57), order='F')
    nu_NOp_sum_temp = np.zeros((72, 144, 57), order='F')
    nu_e_sum_temp = np.zeros((72, 144, 57), order='F')

    # Conductivities
    pedersen_con_temp = np.zeros((72, 144, 57), order='F')
    hall_con_temp = np.zeros((72, 144, 57), order='F')
    parallel_con_temp = np.zeros((72, 144, 57), order='F')

    # Heating rates
    Joule_Heating_temp = np.zeros((72, 144, 57), order='F')
    Frictional_Heating_temp = np.zeros((72, 144, 57), order='F')
    Ohmic_Heating_temp = np.zeros((72, 144, 57), order='F')

    # Cross sections
    C_Op_temp = np.zeros((72, 144, 57), order='F')
    C_O2p_temp = np.zeros((72, 144, 57), order='F')
    C_NOp_temp = np.zeros((72, 144, 57), order='F')
    C_ion_temp = np.zeros((72, 144, 57), order='F')

    # Perpendicular currents
    J_pedersen_temp = np.zeros((72, 144, 57), order='F')
    J_hall_temp = np.zeros((72, 144, 57), order='F')
    J_ohmic_temp = np.zeros((72, 144, 57), order='F')
    J_dens_temp = np.zeros((72, 144, 57), order='F')

    # Distinguish Map from Vertical profile
    lev_range = 0
    lat_range = 0
    lon_range = 0

    lev_start = 0
    lat_start = 0
    lon_start = 0

    # Lat - Alt map profile
    if lat_value == -1 and pressure_level == -1:
        lev_range = len(glev_in) - 1
        lat_range = len(glat_in)
        lon_start = lon_value
        lon_range = lon_start + 1

    # Lat - Lon map profile
    if lat_value == -1 and lon_value == -1:
        lev_start = pressure_level
        lev_range = lev_start + 1
        lat_range = len(glat_in)
        lon_range = len(glon_in)

    # Vertical profile
    if lat_value != -1 and lon_value != -1:
        lat_start = lat_value
        lon_start = lon_value
        lat_range = lat_start + 1
        lon_range = lon_start + 1
        lev_range = len(glev_in) - 1

    for lev in range(lev_start, lev_range):
        for lat in range(lat_start, lat_range):
            for lon in range(lon_start, lon_range):
                # ########################## COLLISION FREQUENCIES ##########################
                # nu-Op = nu(Op-N2) + nu(Op-O2) + nu(Op-O) (in Hz)
                # densities in cm^(-3)
                nu_Op_N2 = 6.82 * NN2_in[lat, lon, lev] * 10 ** (-10)
                nu_Op_O2 = 6.64 * NO2_in[lat, lon, lev] * 10 ** (-10)

                Tr = (Ti_in[lat, lon, lev] + Tn_in[lat, lon, lev]) / 2  # in kelvin
                nu_Op_O = fb * (3.67 * NO_in[lat, lon, lev] * 10 ** (-11) * Tr ** (1 / 2) * (1 - 0.064 * np.log10(Tr)) ** 2)

                nu_Op_sum_temp[lat, lon, lev] = nu_Op_N2 + nu_Op_O2 + nu_Op_O

                # nu-O2p = nu(O2p-N2) + nu(O2p-O) + nu(O2p-O2) (in Hz)
                # densities in cm^(-3)
                nu_O2p_N2 = 4.13 * NN2_in[lat, lon, lev] * 10 ** (-10)
                nu_O2p_O = 2.31 * NO_in[lat, lon, lev] * 10 ** (-10)
                nu_O2p_O2 = 2.59 * NO2_in[lat, lon, lev] * 10 ** (-11) * Tr ** (1 / 2) * (1 - 0.073 * np.log10(Tr)) ** 2

                nu_O2p_sum_temp[lat, lon, lev] = nu_O2p_N2 + nu_O2p_O + nu_O2p_O2

                # nu-NOp = nu(NOp-N2) + nu(NOp-O) + nu(NOp-O2) (in Hz)
                # densities in cm^(-3)
                nu_NOp_N2 = 4.34 * NN2_in[lat, lon, lev] * 10 ** (-10)
                nu_NOp_O = 2.44 * NO_in[lat, lon, lev] * 10 ** (-10)
                nu_NOp_O2 = 4.27 * NO2_in[lat, lon, lev] * 10 ** (-10)

                nu_NOp_sum_temp[lat, lon, lev] = nu_NOp_N2 + nu_NOp_O + nu_NOp_O2

                # nu-e = nu(e-N2) + nu(e-O) + nu(e-O2) (in Hz)
                # densities in cm^(-3)
                nu_e_N2 = 2.33 * 10 ** (-11) * NN2_in[lat, lon, lev] * Te_in[lat, lon, lev] * (1 - 1.21 * 10 ** (-4) * Te_in[lat, lon, lev])
                nu_e_O2 = 1.82 * 10 ** (-10) * NO2_in[lat, lon, lev] * Te_in[lat, lon, lev] ** (1 / 2) * (1 + 3.6 * 10 ** (-2) *
                                                                                                          Te_in[lat, lon, lev] ** (1 / 2))
                nu_e_O = 8.9 * 10 ** (-11) * NO_in[lat, lon, lev] * Te_in[lat, lon, lev] ** (1 / 2) * (1 + 5.7 * 10 ** (-4) * Te_in[lat, lon, lev])

                nu_e_sum_temp[lat, lon, lev] = nu_e_N2 + nu_e_O2 + nu_e_O
                # ################ GYRO-FREQUENCIES(OMEGAS) ################
                # Magnetic field vector (in tesla)
                B = [Bx_in[lat, lon, lev], By_in[lat, lon, lev], Bz_in[lat, lon, lev]]
                # Magnetic field magnitude (in tesla)
                Bnorm = np.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)
                # Magnetic field unit vector
                bunit = [B[0] / Bnorm, B[1] / Bnorm, B[2] / Bnorm]

                # qe(in coulomb), mk(masses in kg), omegas(in Hz)
                omega_Op = (qe * Bnorm) / mkO
                omega_O2p = (qe * Bnorm) / mkO2
                omega_NOp = (qe * Bnorm) / mkNO
                omega_e = (qe * Bnorm) / me
                # ################## RATIOS ##################
                # dimensionless
                r_Op = nu_Op_sum_temp[lat, lon, lev] / omega_Op
                r_O2p = nu_O2p_sum_temp[lat, lon, lev] / omega_O2p
                r_NOp = nu_NOp_sum_temp[lat, lon, lev] / omega_NOp
                r_e = nu_e_sum_temp[lat, lon, lev] / omega_e
                # ############################# CONDUCTIVITIES #############################
                # Pedersen conductivity (in siemens/meter)
                # qe(in coulomb), B_norm(in tesla), N(densities in m^(-3)), ratios(dimensionless)
                term_a_ped = (Ne_in[lat, lon, lev] * ccm) * (r_e / (1 + r_e ** 2))
                term_b_ped = (NOp_in[lat, lon, lev] * ccm) * (r_Op / (1 + r_Op ** 2))
                term_c_ped = (NO2p_in[lat, lon, lev] * ccm) * (r_O2p / (1 + r_O2p ** 2))
                term_d_ped = (NNOp_in[lat, lon, lev] * ccm) * (r_NOp / (1 + r_NOp ** 2))

                pedersen_con_temp[lat, lon, lev] = (qe / Bnorm) * (term_a_ped + term_b_ped + term_c_ped + term_d_ped)

                # Hall conductivity (in siemens/meter)
                # qe(in coulomb), B_norm(in tesla), N(densities in m^(-3)), ratios(dimensionless)
                term_a_hall = (Ne_in[lat, lon, lev] * ccm) / (1 + r_e ** 2)
                term_b_hall = (NOp_in[lat, lon, lev] * ccm) / (1 + r_Op ** 2)
                term_c_hall = (NO2p_in[lat, lon, lev] * ccm) / (1 + r_O2p ** 2)
                term_d_hall = (NNOp_in[lat, lon, lev] * ccm) / (1 + r_NOp ** 2)

                hall_con_temp[lat, lon, lev] = (qe / Bnorm) * (term_a_hall - term_b_hall - term_c_hall - term_d_hall)

                # Parallel conductivity (in siemens/meter)
                # qe(in coulomb), me(mass in tesla), N(density) (in m^(-3)), collision frequency(in Hz)
                parallel_con_temp[lat, lon, lev] = (Ne_in[lat, lon, lev] * ccm * qe ** 2) / (me * nu_e_sum_temp[lat, lon, lev])

                # ################################ HEATING RATES ################################
                # Electric field vector(in volt/meter)
                E = [Ex_in[lat, lon, lev], Ey_in[lat, lon, lev], Ez_in[lat, lon, lev]]

                # Electric field perpendicular to magnetic field
                # Evert = E cross bunit
                Evertx = E[1] * bunit[2] - E[2] * bunit[1]
                Everty = E[2] * bunit[0] - E[0] * bunit[2]
                Evertz = E[0] * bunit[1] - E[1] * bunit[0]

                # E vertical vector
                Evert = [Evertx, Everty, Evertz]

                # Neutral wind vector(in meter/sec)
                Un = [Unx_in[lat, lon, lev], Uny_in[lat, lon, lev], Unz_in[lat, lon, lev]]

                # Vi perpendicular to magnetic field vector
                Vi_vert = [Vix_in[lat, lon, lev], Viy_in[lat, lon, lev], Viz_in[lat, lon, lev]]

                # Neutral wind perpendicular to magnetic field
                # Unvert = Un cross bunit
                Un_vertx = Un[1] * bunit[2] - Un[2] * bunit[1]
                Un_verty = Un[2] * bunit[0] - Un[0] * bunit[2]
                Un_vertz = Un[0] * bunit[1] - Un[1] * bunit[0]

                # Un perpendicular to magnetic field vector
                Un_vert = [Un_vertx, Un_verty, Un_vertz]

                # Unvert cross B vector
                UnvertXBx = Un_vert[1] * B[2] - Un_vert[2] * B[1]
                UnvertXBy = Un_vert[2] * B[0] - Un_vert[0] * B[2]
                UnvertXBz = Un_vert[0] * B[1] - Un_vert[1] * B[0]
                UnvXB = [UnvertXBx, UnvertXBy, UnvertXBz]

                # Estar: perpendicular electric field in the neutral frame Estar = Evert + Unvert cross B
                # vector addition
                Estar_x = Evert[0] + UnvXB[0]
                Estar_y = Evert[1] + UnvXB[1]
                Estar_z = Evert[2] + UnvXB[2]

                Estar = [Estar_x, Estar_y, Estar_z]

                # ################################## JOULE HEATING ##################################
                # ###################################################################################
                # Joule Heating = qeNe(Vi_vert - Un_vert)dot(E_vert + Un_vert cross B)(in watt/m^3)
                # qe(in coulomb), B(in tesla), E(in volt/meter), Vi,Un(in meter/sec)
                Joule_Heating_temp[lat, lon, lev] = qe * (Ne_in[lat, lon, lev] * ccm) * \
                                               (Vi_vert[0] * Estar[0] - Un_vert[0] * Evert[0] +
                                                Vi_vert[1] * Estar[1] - Un_vert[1] * Evert[1] +
                                                Vi_vert[2] * Estar[2] - Un_vert[2] * Evert[2])
                # ################################# OHMIC HEATING ###################################
                # ###################################################################################
                # Ohmic Heating = sigmaPedersen * |Evert + Unvert cross B|^2 (in watt/m^3)
                # sigmaPedersen(in siemens/meter), E(in volt/meter), Un(in meter/sec), B(in tesla)
                term_ohm_x = (Evert[0] + UnvXB[0]) ** 2
                term_ohm_y = (Evert[1] + UnvXB[1]) ** 2
                term_ohm_z = (Evert[2] + UnvXB[2]) ** 2
                Ohmic_Heating_temp[lat, lon, lev] = pedersen_con_temp[lat, lon, lev] * (term_ohm_x + term_ohm_y + term_ohm_z)
                # ############################### FRICTIONAL HEATING ################################
                # ###################################################################################
                # Frictional Heating = m_ion * nu_ion * N_ion * |Vi_vert - Un_vert|^2 (in watt/m^3)
                # m_ion(in kg), nu_ion(in Hz), N_ion(in m^(-3)), Vi,Un(in meter/sec)
                term_fric_x = (Vi_vert[0] - Un_vert[0]) ** 2
                term_fric_y = (Vi_vert[1] - Un_vert[1]) ** 2
                term_fric_z = (Vi_vert[2] - Un_vert[2]) ** 2
                term_Op = mkO * nu_Op_sum_temp[lat, lon, lev] * (NOp_in[lat, lon, lev] * ccm)
                term_O2p = mkO2 * nu_O2p_sum_temp[lat, lon, lev] * (NO2p_in[lat, lon, lev] * ccm)
                term_NOp = mkNO * nu_NOp_sum_temp[lat, lon, lev] * (NNOp_in[lat, lon, lev] * ccm)
                Frictional_Heating_temp[lat, lon, lev] = (term_Op + term_O2p + term_NOp) * (term_fric_x + term_fric_y + term_fric_z)

                # ############################ CROSS SECTIONS ############################
                # C = (nu_ion / N_neutral) / (sqrt(2 * boltzmann * T_i / m_ion)) (in m^2)
                # nu(in Hz), N(in m^(-3)), T(in kelvin), mass(in kg)
                N_neutral = NO_in[lat, lon, lev] + NO2_in[lat, lon, lev] + NN2_in[lat, lon, lev]
                N_neutral = N_neutral * ccm
                # ####### O+ #######
                C_Op_temp[lat, lon, lev] = (nu_Op_sum_temp[lat, lon, lev] / N_neutral) / (np.sqrt(2 * boltzmann * Ti_in[lat, lon, lev] / mkO))
                # ####### O2+ #######
                C_O2p_temp[lat, lon, lev] = (nu_O2p_sum_temp[lat, lon, lev] / N_neutral) / (np.sqrt(2 * boltzmann * Ti_in[lat, lon, lev] / mkO2))
                # ####### NO+ #######
                C_NOp_temp[lat, lon, lev] = (nu_NOp_sum_temp[lat, lon, lev] / N_neutral) / (np.sqrt(2 * boltzmann * Ti_in[lat, lon, lev] / mkNO))
                # ####### ION #######
                nu_ion = nu_Op_sum_temp[lat, lon, lev] + nu_O2p_sum_temp[lat, lon, lev] + nu_NOp_sum_temp[lat, lon, lev]
                # Average collision frequency
                nu_ion = nu_ion / 3  # in Hz
                m_ion = mkO + mkO2 + mkNO
                # Average mass
                m_ion = m_ion / 3  # in kg
                # Because measurements for each species cannot be made we assume an average ion cross section
                C_ion_temp[lat, lon, lev] = (nu_ion / N_neutral) / (np.sqrt(2 * boltzmann * Ti_in[lat, lon, lev] / m_ion))

                # ################################# PERPENDICULAR CURRENTS ##############################
                # ###############################  1st Methodology - Ohms law ###########################
                # Pedersen current = sigmaPedersen * E_star = (Evert + Unvert cross B) (in ampere/meter^2)
                # sigmaPedersen(in siemens/meter), E(in volt/meter)
                J_px = pedersen_con_temp[lat, lon, lev] * Estar[0]
                J_py = pedersen_con_temp[lat, lon, lev] * Estar[1]
                J_pz = pedersen_con_temp[lat, lon, lev] * Estar[2]

                J_pedersen_temp[lat, lon, lev] = np.sqrt(J_px ** 2 + J_py ** 2 + J_pz ** 2)

                # b_unit cross E_star
                x1 = bunit[1] * Estar[2] - bunit[2] * Estar[1]
                y1 = bunit[2] * Estar[0] - bunit[0] * Estar[2]
                z1 = bunit[0] * Estar[1] - bunit[1] * Estar[0]

                # Hall current = sigmaHall * (b_unit cross E_star) (in ampere/meter^2)
                # sigmaHall(in siemens/meter), E(in volt/meter)
                J_hx = hall_con_temp[lat, lon, lev] * x1
                J_hy = hall_con_temp[lat, lon, lev] * y1
                J_hz = hall_con_temp[lat, lon, lev] * z1

                J_hall_temp[lat, lon, lev] = np.sqrt(J_hx ** 2 + J_hy ** 2 + J_hz ** 2)

                # Ohmic current = |JPedersen + JHall|
                J_ohmic_temp[lat, lon, lev] = np.sqrt((J_px + J_hx) ** 2 + (J_py + J_hy) ** 2 + (J_pz + J_hz) ** 2)

                # ####################### 2nd Methodology - Current definition #######################
                # J_density = qe * Ne * (Vi_vert_star - Ve_vert_star) (in neutral frame) or
                # J_density = qe * Ne * (Vi_vert - Ve_vert) (in ECEF frame) in ampere/meter^2

                # Ve_vert (electron velocity perpendicular to magnetic field) in meter/sec
                # Ve_vert_star = (E_star cross B) / |B|^2
                # Ve_vert = (E_star cross B) / |B|^2 + Un_vert
                # ####################################################################################
                # Estar cross B
                x2 = Estar[1] * B[2] - Estar[2] * B[1]
                y2 = Estar[2] * B[0] - Estar[0] * B[2]
                z2 = Estar[0] * B[1] - Estar[1] * B[0]

                x2 = x2 / Bnorm ** 2
                y2 = y2 / Bnorm ** 2
                z2 = z2 / Bnorm ** 2

                Ve_vertx = x2 + Un_vert[0]
                Ve_verty = y2 + Un_vert[1]
                Ve_vertz = z2 + Un_vert[2]

                Ve_vert = [Ve_vertx, Ve_verty, Ve_vertz]

                # qe(in coulomb), Ne(in m^(-3)), Vi,Ve(in meter/sec)
                J_denx = qe * (Ne_in[lat, lon, lev] * ccm) * (Vi_vert[0] - Ve_vert[0])
                J_deny = qe * (Ne_in[lat, lon, lev] * ccm) * (Vi_vert[1] - Ve_vert[1])
                J_denz = qe * (Ne_in[lat, lon, lev] * ccm) * (Vi_vert[2] - Ve_vert[2])
                J_dens_temp[lat, lon, lev] = np.sqrt(J_denx ** 2 + J_deny ** 2 + J_denz ** 2)

    print('Products calculated in: ', time.time() - start_time, ' sec!')
    print(' ')

    return Joule_Heating_temp, Ohmic_Heating_temp, Frictional_Heating_temp, pedersen_con_temp, hall_con_temp, parallel_con_temp, nu_Op_sum_temp, \
           nu_O2p_sum_temp, nu_NOp_sum_temp, nu_e_sum_temp, C_Op_temp, C_O2p_temp, C_NOp_temp, C_ion_temp, J_pedersen_temp, J_hall_temp, \
           J_ohmic_temp, J_dens_temp
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END OF PRODUCTS CALCULATION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ VERTICAL PROFILE PLOTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def plot_heating_rates(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_"*40, "Plotting.....", "_"*40, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=Ohmic_Heating[lat, lon, :-1], y=heights[lat, lon, :-1], name="Ohmic Heating", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=Frictional_Heating[lat, lon, :-1], y=heights[lat, lon, :-1], name="Frictional Heating", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=Joule_Heating[lat, lon, :-1], y=heights[lat, lon, :-1], name="Joule Heating", mode='lines',
                             line=dict(shape='spline', color='green')))

    x_range = max(Joule_Heating[lat, lon, :-1])

    # updating the layout of the figure
    fig.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 10)), xaxis=dict(range=[0, x_range + x_range/4]),
                      xaxis_title="$(W/m^{3})$", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Heating Rates' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig.show()


def plot_collisions(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    mean_ion_collisions = (nu_Op_sum[lat, lon, :-1] + nu_O2p_sum[lat, lon, :-1] + nu_NOp_sum[lat, lon, :-1]) / 3

    fig1 = go.Figure()

    # adding the various plots
    fig1.add_trace(go.Scatter(x=nu_Op_sum[lat, lon, :-1], y=heights[lat, lon, :-1], name="$ν_{O^{+}}$", mode='lines',
                              line=dict(shape='spline', color='red')))
    fig1.add_trace(go.Scatter(x=nu_O2p_sum[lat, lon, :-1], y=heights[lat, lon, :-1], name="$ν_{O2^{+}}$", mode='lines',
                              line=dict(shape='spline', color='blue')))
    fig1.add_trace(go.Scatter(x=nu_NOp_sum[lat, lon, :-1], y=heights[lat, lon, :-1], name="$ν_{NO^{+}}$", mode='lines',
                              line=dict(shape='spline', color='yellow')))
    fig1.add_trace(go.Scatter(x=mean_ion_collisions, y=heights[lat, lon, :-1], name="$ν_{in}$", mode='lines',
                              line=dict(shape='spline', color='orange')))
    fig1.add_trace(go.Scatter(x=nu_e_sum[lat, lon, :-1], y=heights[lat, lon, :-1], name="$ν_{en}$", mode='lines',
                              line=dict(shape='spline', color='purple')))

    # updating the layout of the figure
    fig1.update_layout(xaxis_type="log", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                       tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 10)), xaxis_title="$Frequency \ (Hz)$",
                       yaxis_title="$Altitude \ (km)$", width=800, height=650,
                       title={'text': 'Collision-Gyro Frequencies' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig1.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig1.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig1.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig1.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig1.show()


def plot_conductivities(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=pedersen_con[lat, lon, :-1], y=heights[lat, lon, :-1], name="σPedersen", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=hall_con[lat, lon, :-1], y=heights[lat, lon, :-1], name="σHall", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=parallel_con[lat, lon, :-1], y=heights[lat, lon, :-1], name="σParallel", mode='lines',
                             line=dict(shape='spline', color='green')))

    # updating the layout of the figure
    fig.update_layout(xaxis_type="log", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 10)), xaxis_title="$(S/m)$", yaxis_title="$Altitude \ (km)$",
                      width=800, height=650, title={'text': 'Conductivities' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig.show()


def plot_currents(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=J_pedersen[lat, lon, :-1], y=heights[lat, lon, :-1], name="Pedersen Current", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=J_hall[lat, lon, :-1], y=heights[lat, lon, :-1], name="Hall Current", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=J_ohmic[lat, lon, :-1], y=heights[lat, lon, :-1], name="Ohmic Current", mode='lines',
                             line=dict(shape='spline', color='green')))
    fig.add_trace(go.Scatter(x=J_dens[lat, lon, :-1], y=heights[lat, lon, :-1], name="Densities Current", mode='lines',
                             line=dict(shape='spline', color='black')))
    x_range = max(J_ohmic[lat, lon, :-1])
    # updating the layout of the figure
    fig.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 10)), xaxis=dict(range=[0, x_range + x_range/4]),
                      xaxis_title="$(A/m^{2})$", yaxis_title="$Altitude \ (km)$", width=900, height=750,
                      title={'text': 'Perpendicular Currents' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig.show()


def plot_cross_sections(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=C_Op[lat, lon, :-1], y=heights[lat, lon, :-1], name="$O^{+}$", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=C_O2p[lat, lon, :-1], y=heights[lat, lon, :-1], name="$O_{2}^{+}$", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=C_NOp[lat, lon, :-1], y=heights[lat, lon, :-1], name="$NO^{+}$", mode='lines',
                             line=dict(shape='spline', color='green')))
    fig.add_trace(go.Scatter(x=C_ion[lat, lon, :-1], y=heights[lat, lon, :-1], name="$Avg$", mode='lines',
                             line=dict(shape='spline', color='black')))

    x_range = max(C_Op[lat, lon, :-1])

    # updating the layout of the figure
    fig.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 10)), xaxis=dict(range=[0, x_range + x_range/4]),
                      xaxis_title="$(m^{2})$", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Cross Sections' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig.show()


def plot_heating_rates_error(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=Ohmic_Heating_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="Ohmic Heating with error", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=Frictional_Heating_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="Frictional Heating with error", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=Joule_Heating_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="Joule Heating with error", mode='lines',
                             line=dict(shape='spline', color='green')))

    # updating the layout of the figure
    fig.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                      xaxis_title="$(W/m^{3})$", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Heating Rates With Error' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig.show()


def plot_heating_rates_rel_error(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    Ohmic_rel = abs((Ohmic_Heating[lat, lon, :-1] - Ohmic_Heating_noisy[lat, lon, :-1]) / Ohmic_Heating[lat, lon, :-1])
    Frict_rel = abs((Frictional_Heating[lat, lon, :-1] - Frictional_Heating_noisy[lat, lon, :-1]) / Frictional_Heating[lat, lon, :-1])
    Joule_rel = abs((Joule_Heating[lat, lon, :-1] - Joule_Heating_noisy[lat, lon, :-1]) / Joule_Heating[lat, lon, :-1])

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=Ohmic_rel, y=heights[lat, lon, :-1], name="Ohmic Heating error", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=Frict_rel, y=heights[lat, lon, :-1], name="Frictional Heating error", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=Joule_rel, y=heights[lat, lon, :-1], name="Joule Heating error", mode='lines',
                             line=dict(shape='spline', color='green')))

    # updating the layout of the figure
    fig.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', xaxis=dict(range=[0, 0.4]),
                      yaxis=dict(range=[min_alt, max_alt], tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                      xaxis_title="", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Heating Rates Relative Error' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig.show()


def plot_collisions_error(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=nu_Op_sum_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="$ν_{O^{+}}$ with error", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=nu_O2p_sum_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="$ν_{O2^{+}}$ with error", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=nu_NOp_sum_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="$ν_{NO^{+}}$ with error", mode='lines',
                             line=dict(shape='spline', color='green')))
    fig.add_trace(go.Scatter(x=nu_e_sum_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="$ν_{e}$ with error", mode='lines',
                             line=dict(shape='spline', color='purple')))

    # updating the layout of the figure
    fig.update_layout(xaxis_type="log", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                      xaxis_title="$Frequency \ (Hz)$", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Collision Frequencies with Error' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig.show()


def plot_collisions_rel_error(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    nuOp_rel = abs((nu_Op_sum[lat, lon, :-1] - nu_Op_sum_noisy[lat, lon, :-1]) / nu_Op_sum[lat, lon, :-1])
    nuO2p_rel = abs((nu_O2p_sum[lat, lon, :-1] - nu_O2p_sum_noisy[lat, lon, :-1]) / nu_O2p_sum[lat, lon, :-1])
    nuNOp_rel = abs((nu_NOp_sum[lat, lon, :-1] - nu_NOp_sum_noisy[lat, lon, :-1]) / nu_NOp_sum[lat, lon, :-1])
    nu_ion = (nu_Op_sum[lat, lon, :-1] + nu_O2p_sum[lat, lon, :-1] + nu_NOp_sum[lat, lon, :-1]) / 3
    nuion_error = (nu_Op_sum_noisy[lat, lon, :-1] + nu_O2p_sum_noisy[lat, lon, :-1] + nu_NOp_sum_noisy[lat, lon, :-1]) / 3
    nuion_rel = abs((nu_ion - nuion_error) / nu_ion)
    nue_rel = abs((nu_e_sum[lat, lon, :-1] - nu_e_sum_noisy[lat, lon, :-1]) / nu_e_sum[lat, lon, :-1])

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=nuOp_rel, y=heights[lat, lon, :-1], name="$ν_{O^{+}}$ error", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=nuO2p_rel, y=heights[lat, lon, :-1], name="$ν_{O2^{+}}$ error", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=nuNOp_rel, y=heights[lat, lon, :-1], name="$ν_{NO^{+}}$ error", mode='lines',
                             line=dict(shape='spline', color='green')))
    fig.add_trace(go.Scatter(x=nue_rel, y=heights[lat, lon, :-1], name="$ν_{e}$ error", mode='lines',
                             line=dict(shape='spline', color='purple')))
    fig.add_trace(go.Scatter(x=nuion_rel, y=heights[lat, lon, :-1], name="$ν_{ion}$ error", mode='lines',
                             line=dict(shape='spline', color='yellow')))

    # updating the layout of the figure
    fig.update_layout(xaxis_type="log", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt], tickmode='array',
                      tickvals=np.arange(min_alt, max_alt + 5, 5)), xaxis_title="", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Collision Frequencies Relative Error' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig.show()


def plot_conductivities_error(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=pedersen_con_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="σPedersen with error", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=hall_con_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="σHall with error", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=parallel_con_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="σParallel with error",
                             mode='lines', line=dict(shape='spline', color='green'), visible="legendonly"))
    # updating the layout of the figure
    fig.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                      xaxis_title="$(S/m)$", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Conductivities With Error' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig.show()


def plot_conductivities_rel_error(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    pedersen_rel = abs((pedersen_con[lat, lon, :-1] - pedersen_con_noisy[lat, lon, :-1]) / pedersen_con[lat, lon, :-1])
    hall_rel = abs((hall_con[lat, lon, :-1] - hall_con_noisy[lat, lon, :-1]) / hall_con[lat, lon, :-1])
    parallel_rel = abs((parallel_con[lat, lon, :-1] - parallel_con_noisy[lat, lon, :-1]) / parallel_con[lat, lon, :-1])

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=pedersen_rel, y=heights[lat, lon, :-1], name="σPedersen error", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=hall_rel, y=heights[lat, lon, :-1], name="σHall error", mode='lines', line=dict(shape='spline', color='blue'),
                             visible="legendonly"))
    fig.add_trace(go.Scatter(x=parallel_rel, y=heights[lat, lon, :-1], name="σParallel error", mode='lines',
                             line=dict(shape='spline', color='green')))

    # updating the layout of the figure
    fig.update_layout(xaxis_type="log", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                      xaxis_title="", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Conductivities Relative Error' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig.show()


def plot_currents_error(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=J_pedersen_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="Pedersen Current with error", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=J_hall_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="Hall Current with error", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=J_ohmic_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="Ohmic current with error", mode='lines',
                             line=dict(shape='spline', color='green')))
    fig.add_trace(go.Scatter(x=J_dens_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="Densities current with error", mode='lines',
                             line=dict(shape='spline', color='black')))

    x_range = max(J_dens_noisy[lat, lon, :-1])

    # updating the layout of the figure
    fig.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)), xaxis=dict(range=[0, x_range + x_range/8]),
                      xaxis_title="$(A/m^{2})$", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Perpendicular Currents Absolute Error' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig.show()


def plot_currents_rel_error(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    Jp_rel = abs((J_pedersen[lat, lon, :-1] - J_pedersen_noisy[lat, lon, :-1]) / J_pedersen[lat, lon, :-1])
    Jh_rel = abs((J_hall[lat, lon, :-1] - J_hall_noisy[lat, lon, :-1]) / J_hall[lat, lon, :-1])
    Johmic_rel = abs((J_ohmic[lat, lon, :-1] - J_ohmic_noisy[lat, lon, :-1]) / J_ohmic[lat, lon, :-1])
    Jdens_rel = abs((J_dens[lat, lon, :-1] - J_dens_noisy[lat, lon, :-1]) / J_dens[lat, lon, :-1])

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=Jp_rel, y=heights[lat, lon, :-1], name="Pedersen Current error", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=Jh_rel, y=heights[lat, lon, :-1], name="Hall Current error", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=Johmic_rel, y=heights[lat, lon, :-1], name="Ohmic current error", mode='lines',
                             line=dict(shape='spline', color='green')))
    fig.add_trace(go.Scatter(x=Jdens_rel, y=heights[lat, lon, :-1], name="Densities current error", mode='lines',
                             line=dict(shape='spline', color='black')))

    # updating the layout of the figure
    fig.update_layout(xaxis_type="log", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                      xaxis_title="", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Perpendicular Currents Relative Error' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig.show()


def plot_csections_error(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=C_Op_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="$O^{+}$", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=C_O2p_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="$O_{2}^{+}$", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=C_NOp_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="$NO^{+}$", mode='lines',
                             line=dict(shape='spline', color='green')))
    fig.add_trace(go.Scatter(x=C_ion_noisy[lat, lon, :-1], y=heights[lat, lon, :-1], name="$Avg$", mode='lines',
                             line=dict(shape='spline', color='black')))

    # updating the layout of the figure
    fig.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)), xaxis=dict(range=[0, max(C_Op_noisy[lat, lon, :-1])]),
                      xaxis_title="$(m^{2})$", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Cross Sections With Error' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig.show()


def plot_csections_rel_error(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    COp_rel = abs((C_Op[lat, lon, :-1] - C_Op_noisy[lat, lon, :-1]) / C_Op[lat, lon, :-1])
    CO2p_rel = abs((C_O2p[lat, lon, :-1] - C_O2p_noisy[lat, lon, :-1]) / C_O2p[lat, lon, :-1])
    CNOp_rel = abs((C_NOp[lat, lon, :-1] - C_NOp_noisy[lat, lon, :-1]) / C_NOp[lat, lon, :-1])
    Cion_rel = abs((C_ion[lat, lon, :-1] - C_ion_noisy[lat, lon, :-1]) / C_ion[lat, lon, :-1])

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=COp_rel, y=heights[lat, lon, :-1], name="$O^{+}$", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=CO2p_rel, y=heights[lat, lon, :-1], name="$O_{2}^{+}$", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=CNOp_rel, y=heights[lat, lon, :-1], name="$NO^{+}$", mode='lines',
                             line=dict(shape='spline', color='green')))
    fig.add_trace(go.Scatter(x=Cion_rel, y=heights[lat, lon, :-1], name="$Avg$", mode='lines',
                             line=dict(shape='spline', color='black')))

    # updating the layout of the figure
    fig.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                      xaxis_title="", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Cross Sections Relative Error' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig.show()
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END OF VERTICAL PROFILE PLOTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ LATITUDE-LONGITUDE MAP PROFILE PLOTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def mapll_heating_rates_plot(pressure_level, night_shade):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lev = pressure_level

    # Joule Heating
    fig1 = plt.figure(figsize=(13, 13))
    ax1 = fig1.add_subplot(1, 1, 1, aspect='equal')

    m1 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m1.drawcoastlines()

    sc1 = m1.imshow(Joule_Heating[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m1.nightshade(map_time, alpha=0.3)

    m1.drawparallels(np.arange(-90., 91., 5.))
    m1.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Joule Heating' + title)

    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    cbar = plt.colorbar(sc1, cax=cax1)
    cbar.set_label(label='$(W/m^{3})$', size='large', weight='bold', rotation=270, labelpad=30)

    # Ohmic Heating
    fig2 = plt.figure(figsize=(13, 13))
    ax2 = fig2.add_subplot(1, 1, 1, aspect='equal')

    m2 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m2.drawcoastlines()

    sc2 = m2.imshow(Ohmic_Heating[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m2.nightshade(map_time, alpha=0.3)

    m2.drawparallels(np.arange(-90., 91., 5.))
    m2.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Ohmic Heating' + title)

    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    cbar = plt.colorbar(sc2, cax=cax2)
    cbar.set_label(label='$(W/m^{3})$', size='large', weight='bold', rotation=270, labelpad=30)

    # Frictional Heating
    fig3 = plt.figure(figsize=(13, 13))
    ax3 = fig3.add_subplot(1, 1, 1, aspect='equal')

    m3 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m3.drawcoastlines()

    sc3 = m3.imshow(Frictional_Heating[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m3.nightshade(map_time, alpha=0.3)

    m3.drawparallels(np.arange(-90., 91., 5.))
    m3.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Frictional Heating' + title)

    divider = make_axes_locatable(ax3)
    cax3 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    cbar = plt.colorbar(sc3, cax=cax3)
    cbar.set_label(label='$(W/m^{3})$', size='large', weight='bold', rotation=270, labelpad=30)

    plt.show()


def mapll_collisions_plot(pressure_level, night_shade):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lev = pressure_level

    # Average Ion Collision Frequency
    fig1 = plt.figure(figsize=(13, 13))
    ax1 = fig1.add_subplot(1, 1, 1, aspect='equal')

    m1 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m1.drawcoastlines()

    sc1 = m1.imshow((nu_Op_sum[:, :, lev] + nu_O2p_sum[:, :, lev] + nu_NOp_sum[:, :, lev]) / 3, cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m1.nightshade(map_time, alpha=0.3)

    m1.drawparallels(np.arange(-90., 91., 5.))
    m1.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Average Ion Collision Frequency' + title)

    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    cbar = plt.colorbar(sc1, cax=cax1)
    cbar.set_label(label='$(Hz)$', size='large', weight='bold', rotation=270, labelpad=30)

    # Electron Collision Frequency
    fig2 = plt.figure(figsize=(13, 13))
    ax2 = fig2.add_subplot(1, 1, 1, aspect='equal')

    m2 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m2.drawcoastlines()

    sc2 = m2.imshow(nu_e_sum[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m2.nightshade(map_time, alpha=0.3)

    m2.drawparallels(np.arange(-90., 91., 5.))
    m2.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Electron Collision Frequency' + title)

    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    cbar = plt.colorbar(sc2, cax=cax2)
    cbar.set_label(label='$(Hz)$', size='large', weight='bold', rotation=270, labelpad=30)

    plt.show()


def mapll_conductivities_plot(pressure_level, night_shade):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lev = pressure_level

    # Pedersen Conductivity
    fig1 = plt.figure(figsize=(13, 13))
    ax1 = fig1.add_subplot(1, 1, 1, aspect='equal')

    m1 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m1.drawcoastlines()

    sc1 = m1.imshow(pedersen_con[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m1.nightshade(map_time, alpha=0.3)

    m1.drawparallels(np.arange(-90., 91., 5.))
    m1.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Pedersen Conductivity' + title)

    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    cbar = plt.colorbar(sc1, cax=cax1)
    cbar.set_label(label='$(S/m)$', size='large', weight='bold', rotation=270, labelpad=30)

    # Hall Conductivity
    fig2 = plt.figure(figsize=(13, 13))
    ax2 = fig2.add_subplot(1, 1, 1, aspect='equal')

    m2 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m2.drawcoastlines()

    sc2 = m2.imshow(hall_con[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m2.nightshade(map_time, alpha=0.3)

    m2.drawparallels(np.arange(-90., 91., 5.))
    m2.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Hall Conductivity' + title)

    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    cbar = plt.colorbar(sc2, cax=cax2)
    cbar.set_label(label='$(S/m)$', size='large', weight='bold', rotation=270, labelpad=30)

    # Parallel Conductivity
    fig3 = plt.figure(figsize=(13, 13))
    ax3 = fig3.add_subplot(1, 1, 1, aspect='equal')

    m3 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m3.drawcoastlines()

    sc3 = m3.imshow(parallel_con[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m3.nightshade(map_time, alpha=0.3)

    m3.drawparallels(np.arange(-90., 91., 5.))
    m3.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Parallel Conductivity' + title)

    divider = make_axes_locatable(ax3)
    cax3 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    cbar = plt.colorbar(sc3, cax=cax3)
    cbar.set_label(label='$(S/m)$', size='large', weight='bold', rotation=270, labelpad=30)

    plt.show()


def mapll_currents_plot(pressure_level, night_shade):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lev = pressure_level

    # Ohmic Current
    fig1 = plt.figure(figsize=(13, 13))
    ax1 = fig1.add_subplot(1, 1, 1, aspect='equal')

    m1 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m1.drawcoastlines()

    sc1 = m1.imshow(J_ohmic[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m1.nightshade(map_time, alpha=0.3)

    m1.drawparallels(np.arange(-90., 91., 5.))
    m1.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Ohmic Current' + title)

    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    cbar = plt.colorbar(sc1, cax=cax1)
    cbar.set_label(label='$(A/m^{2})$', size='large', weight='bold', rotation=270, labelpad=30)

    # Densities Current
    fig2 = plt.figure(figsize=(13, 13))
    ax2 = fig2.add_subplot(1, 1, 1, aspect='equal')

    m2 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m2.drawcoastlines()

    sc2 = m2.imshow(J_dens[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m2.nightshade(map_time, alpha=0.3)

    m2.drawparallels(np.arange(-90., 91., 5.))
    m2.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Densities Current' + title)

    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    cbar = plt.colorbar(sc2, cax=cax2)
    cbar.set_label(label='$(A/m^{2})$', size='large', weight='bold', rotation=270, labelpad=30)

    plt.show()


def mapll_csection_plot(pressure_level, night_shade):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lev = pressure_level

    fig1 = plt.figure(figsize=(13, 13))
    ax1 = fig1.add_subplot(1, 1, 1, aspect='equal')

    m1 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m1.drawcoastlines()

    sc1 = m1.imshow(C_ion[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m1.nightshade(map_time, alpha=0.3)

    m1.drawparallels(np.arange(-90., 91., 5.))
    m1.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Average Ion Cross Section' + title)

    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    cbar = plt.colorbar(sc1, cax=cax1)
    cbar.set_label(label='$(m^{2})$', size='large', weight='bold', rotation=270, labelpad=30)

    plt.show()


def mapll_heating_rates_rel_error_plot(pressure_level, night_shade):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lev = pressure_level

    # Joule Heating
    fig1 = plt.figure(figsize=(13, 13))
    ax1 = fig1.add_subplot(1, 1, 1, aspect='equal')

    m1 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m1.drawcoastlines()

    sc1 = m1.imshow(abs((Joule_Heating_noisy[:, :, lev] - Joule_Heating[:, :, lev]) / Joule_Heating[:, :, lev]), cmap=cm.batlow,
                    interpolation='bicubic')

    if night_shade:
        m1.nightshade(map_time, alpha=0.3)

    m1.drawparallels(np.arange(-90., 91., 5.))
    m1.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Joule Heating Relative Error' + title)

    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    plt.colorbar(sc1, cax=cax1)

    # Ohmic Heating
    fig2 = plt.figure(figsize=(13, 13))
    ax2 = fig2.add_subplot(1, 1, 1, aspect='equal')

    m2 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m2.drawcoastlines()

    sc2 = m2.imshow(abs((Ohmic_Heating_noisy[:, :, lev] - Ohmic_Heating[:, :, lev]) / Ohmic_Heating[:, :, lev]), cmap=cm.batlow,
                    interpolation='bicubic')

    if night_shade:
        m2.nightshade(map_time, alpha=0.3)

    m2.drawparallels(np.arange(-90., 91., 5.))
    m2.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Ohmic Heating Relative Error' + title)

    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    plt.colorbar(sc2, cax=cax2)

    # Frictional Heating
    fig3 = plt.figure(figsize=(13, 13))
    ax3 = fig3.add_subplot(1, 1, 1, aspect='equal')

    m3 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m3.drawcoastlines()

    sc3 = m3.imshow(abs((Frictional_Heating_noisy[:, :, lev] - Frictional_Heating[:, :, lev]) / Frictional_Heating[:, :, lev]), cmap=cm.batlow,
                    interpolation='bicubic')

    if night_shade:
        m3.nightshade(map_time, alpha=0.3)

    m3.drawparallels(np.arange(-90., 91., 5.))
    m3.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Frictional Heating Relative Error' + title)

    divider = make_axes_locatable(ax3)
    cax3 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    plt.colorbar(sc3, cax=cax3)

    plt.show()


def mapll_collisions_rel_error_plot(pressure_level, night_shade):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lev = pressure_level

    # Average Ion Collisions Relative Error
    fig1 = plt.figure(figsize=(13, 13))
    ax1 = fig1.add_subplot(1, 1, 1, aspect='equal')

    m1 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m1.drawcoastlines()

    sc1 = m1.imshow(abs((((nu_Op_sum_noisy[:, :, lev] + nu_O2p_sum_noisy[:, :, lev] + nu_NOp_sum_noisy[:, :, lev]) / 3) -
                    ((nu_Op_sum[:, :, lev] + nu_O2p_sum[:, :, lev] + nu_NOp_sum[:, :, lev]) / 3)) /
                    ((nu_Op_sum[:, :, lev] + nu_O2p_sum[:, :, lev] + nu_NOp_sum[:, :, lev]) / 3)),
                    cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m1.nightshade(map_time, alpha=0.3)

    m1.drawparallels(np.arange(-90., 91., 5.))
    m1.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Average Ion Collision Frequency Relative Error' + title)

    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    plt.colorbar(sc1, cax=cax1)

    # Electron Collision Frequency Relative Error
    fig2 = plt.figure(figsize=(13, 13))
    ax2 = fig2.add_subplot(1, 1, 1, aspect='equal')

    m2 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m2.drawcoastlines()

    sc2 = m2.imshow(abs((nu_e_sum_noisy[:, :, lev] - nu_e_sum[:, :, lev]) / nu_e_sum[:, :, lev]), cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m2.nightshade(map_time, alpha=0.3)

    m2.drawparallels(np.arange(-90., 91., 5.))
    m2.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Electron Collision Frequency Relative Error' + title)

    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    plt.colorbar(sc2, cax=cax2)

    plt.show()


def mapll_conductivities_rel_error_plot(pressure_level, night_shade):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lev = pressure_level

    # Pedersen Conductivity
    fig1 = plt.figure(figsize=(13, 13))
    ax1 = fig1.add_subplot(1, 1, 1, aspect='equal')

    m1 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m1.drawcoastlines()

    sc1 = m1.imshow(abs((pedersen_con_noisy[:, :, lev] - pedersen_con[:, :, lev]) / pedersen_con[:, :, lev]), cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m1.nightshade(map_time, alpha=0.3)

    m1.drawparallels(np.arange(-90., 91., 5.))
    m1.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Pedersen Conductivity Relative Error' + title)

    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    plt.colorbar(sc1, cax=cax1)

    # Hall Conductivity
    fig2 = plt.figure(figsize=(13, 13))
    ax2 = fig2.add_subplot(1, 1, 1, aspect='equal')

    m2 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m2.drawcoastlines()

    sc2 = m2.imshow(abs((hall_con_noisy[:, :, lev] - hall_con[:, :, lev]) / hall_con[:, :, lev]), cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m2.nightshade(map_time, alpha=0.3)

    m2.drawparallels(np.arange(-90., 91., 5.))
    m2.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Hall Conductivity Relative Error' + title)

    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    plt.colorbar(sc2, cax=cax2)

    # Parallel Conductivity
    fig3 = plt.figure(figsize=(13, 13))
    ax3 = fig3.add_subplot(1, 1, 1, aspect='equal')

    m3 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m3.drawcoastlines()

    sc3 = m3.imshow(abs((parallel_con_noisy[:, :, lev] - parallel_con[:, :, lev]) / parallel_con[:, :, lev]), cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m3.nightshade(map_time, alpha=0.3)

    m3.drawparallels(np.arange(-90., 91., 5.))
    m3.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Parallel Conductivity Relative Error' + title)

    divider = make_axes_locatable(ax3)
    cax3 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    plt.colorbar(sc3, cax=cax3)

    plt.show()


def mapll_currents_rel_error_plot(pressure_level, night_shade):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lev = pressure_level

    # Ohmic Current
    fig1 = plt.figure(figsize=(13, 13))
    ax1 = fig1.add_subplot(1, 1, 1, aspect='equal')

    m1 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m1.drawcoastlines()

    sc1 = m1.imshow(abs((J_ohmic_noisy[:, :, lev] - J_ohmic[:, :, lev]) / J_ohmic[:, :, lev]), cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m1.nightshade(map_time, alpha=0.3)

    m1.drawparallels(np.arange(-90., 91., 5.))
    m1.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Ohmic Current Relative Error' + title)

    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    plt.colorbar(sc1, cax=cax1)

    # Densities Current
    fig2 = plt.figure(figsize=(13, 13))
    ax2 = fig2.add_subplot(1, 1, 1, aspect='equal')

    m2 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m2.drawcoastlines()

    sc2 = m2.imshow(abs((J_dens_noisy[:, :, lev] - J_dens[:, :, lev]) / J_dens[:, :, lev]), cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m2.nightshade(map_time, alpha=0.3)

    m2.drawparallels(np.arange(-90., 91., 5.))
    m2.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Densities Current Relative Error' + title)

    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    plt.colorbar(sc2, cax=cax2)

    plt.show()


def mapll_csection_rel_error_plot(pressure_level, night_shade):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lev = pressure_level

    fig1 = plt.figure(figsize=(13, 13))
    ax1 = fig1.add_subplot(1, 1, 1, aspect='equal')

    m1 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m1.drawcoastlines()

    sc1 = m1.imshow(abs((C_ion_noisy[:, :, lev] - C_ion[:, :, lev]) / C_ion[:, :, lev]), cmap=cm.batlow, interpolation='bicubic')

    if night_shade:
        m1.nightshade(map_time, alpha=0.3)

    m1.drawparallels(np.arange(-90., 91., 5.))
    m1.drawmeridians(np.arange(-180., 181., 10.))

    plt.xticks(np.arange(-180., 181., 60.))
    plt.yticks(np.arange(-90., 91., 30.))
    plt.xlabel('$Longitude \ (deg)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Average Ion Cross Section Relative Error' + title)

    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad=0.2, aspect=15)
    plt.colorbar(sc1, cax=cax1)

    plt.show()
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END OF LATITUDE-LONGITUDE MAP PROFILE PLOTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ LATITUDE-ALTITUDE PROFILE PLOTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def mapla_heating_rates_plot(lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    # Joule Heating
    plt.figure(figsize=(12, 12))
    cp1 = plt.contourf(heights_la[:-1], glat_in[:], Joule_Heating[:, lon, :-1], cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Joule Heating' + title)
    cbar = plt.colorbar(cp1)
    cbar.set_label(label='$(W/m^{3})$', size='large', weight='bold', rotation=270, labelpad=30)

    # Ohmic Heating
    plt.figure(figsize=(12, 12))
    cp2 = plt.contourf(heights_la[:-1], glat_in[:], Ohmic_Heating[:, lon, :-1], cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Ohmic Heating' + title)
    cbar = plt.colorbar(cp2)
    cbar.set_label(label='$(W/m^{3})$', size='large', weight='bold', rotation=270, labelpad=30)

    # Frictional Heating
    plt.figure(figsize=(12, 12))
    cp3 = plt.contourf(heights_la[:-1], glat_in[:], Frictional_Heating[:, lon, :-1], cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Frictional Heating' + title)
    cbar = plt.colorbar(cp3)
    cbar.set_label(label='$(W/m^{3})$', size='large', weight='bold', rotation=270, labelpad=30)

    plt.show()


def mapla_collisions_plot(lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    # Ion Average Collision Frequency
    plt.figure(figsize=(12, 12))
    cp1 = plt.contourf(heights_la[:-1], glat_in[:], (nu_Op_sum[:, lon, :-1] + nu_O2p_sum[:, lon, :-1] + nu_NOp_sum[:, lon, :-1]) / 3,
                       locator=ticker.LogLocator(), cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Average Ion Collision Frequency' + title)
    cbar = plt.colorbar(cp1)
    cbar.set_label(label='$(Hz)$', size='large', weight='bold', rotation=270, labelpad=30)

    # Electron Collision Frequency
    plt.figure(figsize=(12, 12))
    cp2 = plt.contourf(heights_la[:-1], glat_in[:], nu_e_sum[:, lon, :-1], locator=ticker.LogLocator(), cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Electron Collision Frequency' + title)
    cbar = plt.colorbar(cp2)
    cbar.set_label(label='$(Hz)$', size='large', weight='bold', rotation=270, labelpad=30)

    plt.show()


def mapla_conductivities_plot(lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    # Pedersen
    plt.figure(figsize=(12, 12))
    cp1 = plt.contourf(heights_la[:-1], glat_in[:], pedersen_con[:, lon, :-1], cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Pedersen Conductivity' + title)
    cbar = plt.colorbar(cp1)
    cbar.set_label(label='$(S/m)$', size='large', weight='bold', rotation=270, labelpad=30)

    # Hall Conductivity
    plt.figure(figsize=(12, 12))
    cp2 = plt.contourf(heights_la[:-1], glat_in[:], hall_con[:, lon, :-1], cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Hall Conductivity' + title)
    cbar = plt.colorbar(cp2)
    cbar.set_label(label='$(S/m)$', size='large', weight='bold', rotation=270, labelpad=30)

    # Parallel Conductivity
    plt.figure(figsize=(12, 12))
    cp3 = plt.contourf(heights_la[:-1], glat_in[:], parallel_con[:, lon, :-1], locator=ticker.LogLocator(), cmap=cm.batlow,
                       interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')
    plt.title('Parallel Conductivity' + title)
    cbar = plt.colorbar(cp3)
    cbar.set_label(label='$(S/m)$', size='large', weight='bold', rotation=270, labelpad=30)

    plt.show()


def mapla_currents_plot(lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    # Ohmic Current
    plt.figure(figsize=(12, 12))
    cp1 = plt.contourf(heights_la[:-1], glat_in[:], J_ohmic[:, lon, :-1], cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Ohmic Current' + title)
    cbar = plt.colorbar(cp1)
    cbar.set_label(label='$(A/m^{2})$', size='large', weight='bold', rotation=270, labelpad=30)

    # Densities Current
    plt.figure(figsize=(12, 12))
    cp2 = plt.contourf(heights_la[:-1], glat_in[:], J_dens[:, lon, :-1], cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Densities Current' + title)
    cbar = plt.colorbar(cp2)
    cbar.set_label(label='$(A/m^{2})$', size='large', weight='bold', rotation=270, labelpad=30)

    plt.show()


def mapla_cross_section_plot(lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    # Average Ion Cross Section
    plt.figure(figsize=(12, 12))
    cp1 = plt.contourf(heights_la[:-1], glat_in[:], C_ion[:, lon, :-1], cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Average Ion Cross Section' + title)
    cbar = plt.colorbar(cp1)
    cbar.set_label(label='$(m^{2})$', size='large', weight='bold', rotation=270, labelpad=30)

    plt.show()


def mapla_heating_rates_rel_error_plot(lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    # Joule Heating
    plt.figure(figsize=(12, 12))
    cp1 = plt.contourf(heights_la[:-1], glat_in[:], abs((Joule_Heating_noisy[:, lon, :-1] - Joule_Heating[:, lon, :-1]) / Joule_Heating[:, lon, :-1]),
                       locator=ticker.LogLocator(), cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Joule Heating Relative Error' + title)
    cbar = plt.colorbar(cp1)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    # Ohmic Heating
    plt.figure(figsize=(12, 12))
    cp2 = plt.contourf(heights_la[:-1], glat_in[:], abs((Ohmic_Heating_noisy[:, lon, :-1] - Ohmic_Heating[:, lon, :-1]) / Ohmic_Heating[:, lon, :-1]),
                       locator=ticker.LogLocator(), cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Ohmic Heating Relative Error' + title)
    cbar = plt.colorbar(cp2)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    # Frictional Heating
    plt.figure(figsize=(12, 12))
    cp3 = plt.contourf(heights_la[:-1], glat_in[:], abs((Frictional_Heating_noisy[:, lon, :-1] - Frictional_Heating[:, lon, :-1]) /
                                                        Frictional_Heating[:, lon, :-1]),
                       locator=ticker.LogLocator(), cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Frictional Heating Relative Error' + title)
    cbar = plt.colorbar(cp3)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    plt.show()


def mapla_collisions_rel_error_plot(lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    # Ion Average Collision Frequency
    plt.figure(figsize=(12, 12))
    cp1 = plt.contourf(heights_la[:-1], glat_in[:], abs((((nu_Op_sum_noisy[:, lon, :-1] + nu_O2p_sum_noisy[:, lon, :-1] +
                                                           nu_NOp_sum_noisy[:, lon, :-1]) / 3) - ((nu_Op_sum[:, lon, :-1] + nu_O2p_sum[:, lon, :-1] +
                                                                                                   nu_NOp_sum[:, lon, :-1]) / 3)) /
                                                           ((nu_Op_sum[:, lon, :-1] + nu_O2p_sum[:, lon, :-1] + nu_NOp_sum[:, lon, :-1]) / 3)),
                       locator=ticker.LogLocator(), cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Average Ion Collision Frequency Relative Error' + title)
    cbar = plt.colorbar(cp1)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    # Electron Collision Frequency
    plt.figure(figsize=(12, 12))
    cp2 = plt.contourf(heights_la[:-1], glat_in[:], abs((nu_e_sum_noisy[:, lon, :-1] - nu_e_sum[:, lon, :-1]) / nu_e_sum[:, lon, :-1]),
                       locator=ticker.LogLocator(), cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Electron Collision Frequency Relative Error' + title)
    cbar = plt.colorbar(cp2)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    plt.show()


def mapla_conductivities_rel_error_plot(lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    # Pedersen
    plt.figure(figsize=(12, 12))
    cp1 = plt.contourf(heights_la[:-1], glat_in[:], abs((pedersen_con_noisy[:, lon, :-1] - pedersen_con[:, lon, :-1]) / pedersen_con[:, lon, :-1]),
                       locator=ticker.LogLocator(), cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Pedersen Conductivity Relative Error' + title)
    cbar = plt.colorbar(cp1)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    # Hall Conductivity
    plt.figure(figsize=(12, 12))
    cp2 = plt.contourf(heights_la[:-1], glat_in[:], abs((hall_con_noisy[:, lon, :-1] - hall_con[:, lon, :-1]) / hall_con[:, lon, :-1]),
                       locator=ticker.LogLocator(), cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Hall Conductivity Relative Error' + title)
    cbar = plt.colorbar(cp2)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    # Parallel Conductivity
    plt.figure(figsize=(12, 12))
    cp3 = plt.contourf(heights_la[:-1], glat_in[:], abs((parallel_con_noisy[:, lon, :-1] - parallel_con[:, lon, :-1]) / parallel_con[:, lon, :-1]),
                       locator=ticker.LogLocator(), cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Parallel Conductivity Relative Error' + title)
    cbar = plt.colorbar(cp3)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    plt.show()


def mapla_currents_rel_error_plot(lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    # Ohmic Current
    plt.figure(figsize=(12, 12))
    cp1 = plt.contourf(heights_la[:-1], glat_in[:], abs((J_ohmic_noisy[:, lon, :-1] - J_ohmic[:, lon, :-1]) / J_ohmic[:, lon, :-1]),
                       locator=ticker.LogLocator(), cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Ohmic Current Relative Error' + title)
    cbar = plt.colorbar(cp1)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    # Densities Current
    plt.figure(figsize=(12, 12))
    cp2 = plt.contourf(heights_la[:-1], glat_in[:], abs((J_dens_noisy[:, lon, :-1] - J_dens[:, lon, :-1]) / J_dens[:, lon, :-1]),
                       locator=ticker.LogLocator(), cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Densities Current Relative Error' + title)
    cbar = plt.colorbar(cp2)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    plt.show()


def mapla_cross_section_rel_error_plot(lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    # Average Ion Cross Section
    plt.figure(figsize=(12, 12))
    cp1 = plt.contourf(heights_la[:-1], glat_in[:], abs((C_ion_noisy[:, lon, :-1] - C_ion[:, lon, :-1]) / C_ion[:, lon, :-1]),
                       locator=ticker.LogLocator(), cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Average Ion Cross Section Relative Error' + title)
    cbar = plt.colorbar(cp1)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    plt.show()


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ GUI CREATION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def gui():
    # function used to define buttons in gui
    def button(name, tooltip):
        return Sg.Button(name, tooltip=tooltip)

    # dictionaries for getting the necessary indexes needed in TIEGCM file out of user's gui selections
    lat_dictionary = {"-88.75": 0, "-86.25": 1, "-83.75": 2, "-81.25": 3, "-78.75": 4, "-76.25": 5, "-73.75": 6,
                      "-71.25": 7, "-68.75": 8, "-66.25": 9, "-63.75": 10, "-61.25": 11, "-58.75": 12, "-56.25": 13,
                      "-53.75": 14, "-51.25": 15, "-48.75": 16, "-46.25": 17, "-43.75": 18, "-41.25": 19, "-38.75": 20,
                      "-36.25": 21, "-33.75": 22, "-31.25": 23, "-28.75": 24, "-26.25": 25, "-23.75": 26, "-21.25": 27,
                      "-18.75": 28, "-16.25": 29, "-13.75": 30, "-11.25": 31, "-8.75": 32, "-6.25": 33, "-3.75": 34,
                      "-1.25": 35, "1.25": 36, "3.75": 37, "6.25": 38, "8.75": 39, "11.25": 40, "13.75": 41,
                      "16.25": 42, "18.75": 43, "21.25": 44, "23.75": 45, "26.25": 46, "28.75": 47, "31.25": 48,
                      "33.75": 49, "36.25": 50, "38.75": 51, "41.25": 52, "43.75": 53, "46.25": 54, "48.75": 55,
                      "51.25": 56, "53.75": 57, "56.25": 58, "58.75": 59, "61.25": 60, "63.75": 61, "66.25": 62,
                      "68.75": 63, "71.25": 64, "73.75": 65, "76.25": 66, "78.75": 67, "81.25": 68, "83.75": 69,
                      "84.25": 70, "87.75": 71}

    lon_dictionary = {"-180.0": 0, "-177.5": 1, "-175.0": 2, "-172.5": 3, "-170.0": 4, "-167.5": 5, "-165.0": 6,
                      "-162.5": 7, "-160.0": 8, "-157.5": 9, "-155.0": 10, "-152.5": 11, "-150.0": 12, "-147.5": 13,
                      "-145.0": 14, "-142.5": 15, "-140.0": 16, "-137.5": 17, "-135.0": 18, "-132.5": 19, "-130.0": 20,
                      "-127.5": 21, "-125.0": 22, "-122.5": 23, "-120.0": 24, "-117.5": 25, "-115.0": 26, "-112.5": 27,
                      "-110.0": 28, "-107.5": 29, "-105.0": 30, "-102.5": 31, "-100.0": 32, "-97.5": 33, "-95.0": 34,
                      "-92.5": 35, "-90.0": 36, "-87.5": 37, "-85.0": 38, "-82.5": 39, "-80.0": 40, "-77.5": 41,
                      "-75.0": 42, "-72.5": 43, "-70.0": 44, "-67.5": 45, "-65.0": 46, "-62.5": 47, "-60.0": 48,
                      "-57.5": 49, "-55.0": 50, "-52.5": 51, "-50.0": 52, "-47.5": 53, "-45.0": 54, "-42.5": 55,
                      "-40.0": 56, "-37.5": 57, "-35.0": 58, "-32.5": 59, "-30.0": 60, "-27.5": 61, "-25.0": 62,
                      "-22.5": 63, "-20.0": 64, "-17.5": 65, "-15.0": 66, "-12.5": 67, "-10.0": 68, "-7.5": 69,
                      "-5.0": 70, "-2.5": 71, "0.0": 72, "2.5": 73, "5.0": 74, "7.5": 75, "10.0": 76, "12.5": 77,
                      "15.0": 78, "17.5": 79, "20.0": 80, "22.5": 81, "25.0": 82, "27.5": 83, "30.0": 84, "32.5": 85,
                      "35.0": 86, "37.5": 87, "40.0": 88, "42.5": 89, "45.0": 90, "47.5": 91, "50.0": 92, "52.5": 93,
                      "55.0": 94, "57.5": 95, "60.0": 96, "62.5": 97, "65.0": 98, "67.5": 99, "70.0": 100, "72.5": 101,
                      "75.0": 102, "77.5": 103, "80.0": 104, "82.5": 105, "85.0": 106, "87.5": 107, "90.0": 108,
                      "92.5": 109, "95.0": 110, "97.5": 111, "100.0": 112, "102.5": 113, "105.0": 114, "107.5": 115,
                      "110.0": 116, "112.5": 117, "115.0": 118, "117.5": 119, "120.0": 120, "122.5": 121, "125.0": 122,
                      "127.5": 123, "130.0": 124, "132.5": 125, "135.0": 126, "137.5": 127, "140.0": 128, "142.5": 129,
                      "145.0": 130, "147.5": 131, "150.0": 132, "152.5": 133, "155.0": 134, "157.5": 135, "160.0": 136,
                      "162.5": 137, "165.0": 138, "167.5": 139, "170.0": 140, "172.5": 141, "175.0": 142, "177.5": 143}

    # used in plotting choices frame
    col1 = [[Sg.Checkbox("Plot Heating Rates With Error", default=False, tooltip="Plots heating rates with error", key="-HR_abs-")],
            [Sg.Checkbox("Plot Heating Rates Relative Error", default=False, tooltip="Plots heating rates relative error", key="-HR_rel-")],
            [Sg.Checkbox("Plot Collision Frequencies With Error", default=False, tooltip="Plots collision frequencies with error",
                         key="-COL_abs-")],
            [Sg.Checkbox("Plot Collision Frequencies Relative Error", default=False, tooltip="Plots collision frequencies relative error",
                         key="-COL_rel-")],
            [Sg.Checkbox("Plot Conductivities With Error", default=False, tooltip="Plots conductivities with error", key="-CON_abs-")],
            [Sg.Checkbox("Plot Conductivities Relative Error", default=False, tooltip="Plots conductivities relative error", key="-CON_rel-")],
            [Sg.Checkbox("Plot Currents With Error", default=False, tooltip="Plots currents with error", key="-CUR_abs-")],
            [Sg.Checkbox("Plot Currents Relative Error", default=False, tooltip="Plots currents relative error", key="-CUR_rel-")],
            [Sg.Checkbox("Plot Cross Sections With Error", default=False, tooltip="Plots cross sections with error", key="-CR_abs-")],
            [Sg.Checkbox("Plot Cross Sections Relative Error", default=False, tooltip="Plots cross sections relative error", key="-CR_rel-")]]

    col2 = [[Sg.Checkbox("Plot Heating Rates", default=False, tooltip="Plots heating rates", key="-HR-")],
            [Sg.Checkbox("Plot Collision Frequencies", default=False, tooltip="Plots collision frequencies", key="-COL-")],
            [Sg.Checkbox("Plot Conductivities", default=False, tooltip="Plots conductivities", key="-CON-")],
            [Sg.Checkbox("Plot Currents", default=False, tooltip="Plots currents", key="-CUR-")],
            [Sg.Checkbox("Plot Cross Sections", default=False, tooltip="Plots cross sections", key="-CR-")]]

    col3 = [[Sg.Checkbox("Plot Heating Rates", default=False, tooltip="Plots heating rates", key="-HR_mapll-")],
            [Sg.Checkbox("Plot Collision Frequencies", default=False, tooltip="Plots collision frequencies", key="-COL_mapll-")],
            [Sg.Checkbox("Plot Conductivities", default=False, tooltip="Plots conductivities", key="-CON_mapll-")],
            [Sg.Checkbox("Plot Currents", default=False, tooltip="Plots currents", key="-CUR_mapll-")],
            [Sg.Checkbox("Plot Cross Sections", default=False, tooltip="Plots cross sections", key="-CR_mapll-")]]

    col4 = [[Sg.Checkbox("Plot Heating Rates Relative Error", default=False, tooltip="Plots heating rates relative error", key="-HR_mapll_error-")],
            [Sg.Checkbox("Plot Collision Frequencies Relative Error", default=False, tooltip="Plots collision frequencies relative error",
                key="-COL_mapll_error-")],
            [Sg.Checkbox("Plot Conductivities Relative Error", default=False, tooltip="Plots conductivities relative error",
                key="-CON_mapll_error-")],
            [Sg.Checkbox("Plot Currents Relative Error", default=False, tooltip="Plots currents relative error", key="-CUR_mapll_error-")],
            [Sg.Checkbox("Plot Cross Sections Relative Error", default=False, tooltip="Plots cross sections relative error", key="-CR_mapll_error-")]]

    col5 = [[Sg.Checkbox("Plot Heating Rates", default=False, tooltip="Plots heating rates", key="-HR_mapla-")],
            [Sg.Checkbox("Plot Collision Frequencies", default=False, tooltip="Plots collision frequencies", key="-COL_mapla-")],
            [Sg.Checkbox("Plot Conductivities", default=False, tooltip="Plots conductivities", key="-CON_mapla-")],
            [Sg.Checkbox("Plot Currents", default=False, tooltip="Plots currents", key="-CUR_mapla-")],
            [Sg.Checkbox("Plot Cross Sections", default=False, tooltip="Plots cross sections", key="-CR_mapla-")]]

    col6 = [[Sg.Checkbox("Plot Heating Rates Relative Error", default=False, tooltip="Plots heating rates relative error", key="-HR_mapla_error-")],
            [Sg.Checkbox("Plot Collision Frequencies Relative Error", default=False, tooltip="Plots collision frequencies relative error",
                key="-COL_mapla_error-")],
            [Sg.Checkbox("Plot Conductivities Relative Error", default=False, tooltip="Plots conductivities relative error",
                key="-CON_mapla_error-")],
            [Sg.Checkbox("Plot Currents Relative Error", default=False, tooltip="Plots currents relative error", key="-CUR_mapla_error-")],
            [Sg.Checkbox("Plot Cross Sections Relative Error", default=False, tooltip="Plots cross sections relative error", key="-CR_mapla_error-")]]

    templay1 = [[Sg.Text("min altitude(km)", pad=((30, 0), (15, 0))), Sg.Text("max altitude(km)", pad=((30, 10), (15, 0)))],
                [Sg.InputCombo(values=[i for i in range(100, 601, 10)], pad=((40, 0), (10, 20)), size=(10, 1), default_value="110", key="-min_alt-"),
                 Sg.InputCombo(values=[i for i in range(100, 601, 10)], pad=((45, 0), (10, 20)), size=(11, 1), default_value="200", key="-max_alt-")],
                [Sg.Text("Products", pad=((30, 0), (0, 0))), Sg.Text("Errors", pad=((200, 0), (0, 0)))],
                [Sg.Text("_" * 15, pad=((20, 0), (0, 0))), Sg.Text("_" * 30, pad=((105, 0), (0, 0)))],
                [Sg.Column(col2, pad=((0, 0), (0, 100)), scrollable=False), Sg.Column(col1, size=(360, 270), pad=((15, 0), (0, 0)), scrollable=True)]]

    templay2 = [[Sg.Text("Products", pad=((30, 0), (30, 0))), Sg.Text("Errors", pad=((200, 0), (30, 0)))],
                [Sg.Text("_" * 15, pad=((20, 0), (0, 0))), Sg.Text("_" * 30, pad=((105, 0), (0, 0)))],
                [Sg.Column(col3, scrollable=False), Sg.Column(col4, scrollable=False)]]

    templay3 = [[Sg.Text("min altitude(km)", pad=((30, 0), (15, 0))), Sg.Text("max altitude(km)", pad=((30, 10), (15, 0)))],
                [Sg.InputCombo(values=[i for i in range(100, 601, 10)], pad=((40, 0), (10, 20)), size=(10, 1), default_value="110",
                    key="-min_alt_la-"),
                 Sg.InputCombo(values=[i for i in range(100, 601, 10)], pad=((45, 0), (10, 20)), size=(11, 1), default_value="200",
                     key="-max_alt_la-")],
                [Sg.Text("Products", pad=((30, 0), (30, 0))), Sg.Text("Errors", pad=((200, 0), (30, 0)))],
                [Sg.Text("_" * 15, pad=((20, 0), (0, 0))), Sg.Text("_" * 30, pad=((105, 0), (0, 0)))],
                [Sg.Column(col5, scrollable=False), Sg.Column(col6, scrollable=False)]]

    templay4 = [[Sg.TabGroup([[Sg.Tab("Vertical Profile Plots", templay1), Sg.Tab("Map Profile (Lat-Lon) Plots", templay2),
                               Sg.Tab("Map Profile (Lat-Alt) Plots", templay3)]], key="-TABGROUP1-")]]

    lat_values = ["-88.75", "-86.25", "-83.75", "-81.25", "-78.75", "-76.25", "-73.75", "-71.25", "-68.75", "-66.25", "-63.75", "-61.25", "-58.75",
                  "-56.25", "-53.75", "-51.25", "-48.75", "-46.25", "-43.75", "-41.25", "-38.75", "-36.25", "-33.75", "-31.25", "-28.75", "-26.25",
                  "-23.75", "-21.25", "-18.75", "-16.25", "-13.75", "-11.25", "-8.75", "-6.25", "-3.75", "-1.25", "1.25", "3.75", "6.25", "8.75",
                  "11.25", "13.75", "16.25", "18.75", "21.25", "23.75", "26.25", "28.75", "31.25", "33.75", "36.25", "38.75", "41.25", "43.75",
                  "46.25", "48.75", "51.25", "53.75", "56.25", "58.75", "61.25", "63.75", "66.25", "68.75", "71.25", "73.75", "76.25", "78.75",
                  "81.25", "83.75", "84.25", "87.75"]

    lon_values = ["-180.0", "-177.5", "-175.0", "-172.5", "-170.0", "-167.5", "-165.0", "-162.5", "-160.0", "-157.5", "-155.0", "-152.5", "-150.0",
                  "-147.5", "-145.0", "-142.5", "-140.0", "-137.5", "-135.0", "-132.5", "-130.0", "-127.5", "-125.0", "-122.5", "-120.0", "-117.5",
                  "-115.0", "-112.5", "-110.0", "-107.5", "-105.0", "-102.5", "-100.0", "-97.5", "-95.0", "-92.5", "-90.0", "-87.5", "-85.0", "-82.5",
                  "-80.0", "-77.5", "-75.0", "-72.5", "-70.0", "-67.5", "-65.0", "-62.5", "-60.0", "-57.5", "-55.0", "-52.5", "-50.0", "-47.5",
                  "-45.0", "-42.5", "-40.0", "-37.5", "-35.0", "-32.5", "-30.0", "-27.5", "-25.0", "-22.5", "-20.0", "-17.5", "-15.0", "-12.5",
                  "-10.0", "-7.5", "-5.0", "-2.5", "0.0", "2.5", "5.0", "7.5", "10.0", "12.5", "15.0", "17.5", "20.0", "22.5", "25.0", "27.5", "30.0",
                  "32.5", "35.0", "37.5", "40.0", "42.5", "45.0", "47.5", "50.0", "52.5", "55.0", "57.5", "60.0", "62.5", "65.0", "67.5", "70.0",
                  "72.5", "75.0", "77.5", "80.0", "82.5", "85.0", "87.5", "90.0", "92.5", "95.0", "97.5", "100.0", "102.5", "105.0", "107.5", "110.0",
                  "112.5", "115.0", "117.5", "120.0", "122.5", "125.0", "127.5", "130.0", "132.5", "135.0", "137.5", "140.0", "142.5", "145.0",
                  "147.5", "150.0", "152.5", "155.0", "157.5", "160.0", "162.5", "165.0", "167.5", "170.0", "172.5", "175.0", "177.5"]

    # ############################################ VERTICAL PROFILE LAYOUT #############################################
    # ##################################################################################################################
    vert_layout = [[Sg.Text("Latitude(degrees)", pad=((30, 20), (30, 0))), Sg.Text("Longitude(degrees)", pad=((30, 20), (30, 0))),
                    Sg.Text("Timestep", pad=((30, 20), (30, 0)))],
                   [Sg.InputCombo(values=lat_values, default_value="63.75", pad=((30, 0), (0, 20)), size=(14, 1), key="-LAT-"),
                    Sg.InputCombo(values=lon_values, default_value="-57.5", pad=((45, 0), (0, 20)), size=(15, 1), key="-LON-"),
                    Sg.InputCombo(values=[i for i in range(0, 60)], pad=((45, 0), (0, 20)), size=(7, 1), default_value="9", key="-TIME_vert-")]]
    # ########################################### VERTICAL PROFILE LAYOUT END ##########################################
    # ##################################################################################################################

    # ############################################### MAP PROFILE LAYOUT ###############################################
    # ##################################################################################################################
    map_layout = [[Sg.Text("Timestep", pad=((40, 20), (30, 0))), Sg.Text("Pressure level", pad=((30, 20), (30, 0)))],
                  [Sg.InputCombo(values=[i for i in range(0, 60)], pad=((45, 0), (0, 20)), size=(7, 1), default_value="9", key="-TIME_map-"),
                   Sg.InputCombo(values=[i for i in range(0, 56)], pad=((40, 0), (0, 20)), size=(11, 1), default_value="7", key="-Pr_level-")],
                  [Sg.Checkbox("Add nightshade", default=True, tooltip="Adds night region on map", key="-NIGHT-")]]

    map2_latout = [[Sg.Text("Timestep", pad=((40, 20), (30, 0))), Sg.Text("Longitude", pad=((30, 20), (30, 0)))],
                   [Sg.InputCombo(values=[i for i in range(0, 60)], pad=((45, 0), (0, 20)), size=(7, 1), default_value="9", key="-TIME_map2-"),
                    Sg.InputCombo(values=lon_values, pad=((40, 0), (0, 20)), size=(11, 1), default_value="-57.5", key="-Lon_map2-")]]
    # ############################################# MAP PROFILE LAYOUT END #############################################
    # ##################################################################################################################

    # ################################################### MAIN LAYOUT ##################################################
    # ##################################################################################################################
    main_layout = [[Sg.Text("Menu", font=("Helvetica", 30))], [Sg.HSeparator()], [Sg.Text("Choose File")],
                   [Sg.InputCombo(("tiegcm2.0_res2.5_3years_sech_014_JH_QD_AllVars.nc",
                                   "tiegcm2.0_res2.5_3years_sech_015_JH_QD_AllVars.nc",
                                   "tiegcm2.0_res2.5_3years_sech_016_JH_QD_AllVars.nc"),
                    default_value="tiegcm2.0_res2.5_3years_sech_014_JH_QD_AllVars.nc", key="-FILE-")],
                   [Sg.Frame("Choose plots", templay4, pad=((0, 0), (40, 30)))],
                   [Sg.Text("Choose Profile")],
                   [Sg.TabGroup([[Sg.Tab("Vertical Profile", vert_layout), Sg.Tab("Map Profile (Lat-Lon)", map_layout),
                                  Sg.Tab("Map Profile (Lat-Alt)", map2_latout)]], key="-TABGROUP-")],
                   [button("Calculate Products", "Click to start calculation"), button("Calculate Error", "Click to start calculation"),
                    button("Help", "Click for info"),
                    button("Exit", "Click to exit program")]]
    # ################################################# MAIN LAYOUT END ################################################
    # ##################################################################################################################
    # ################################################# Create Window ##################################################
    window = Sg.Window("Joule Heating Statistical Error GUI", main_layout, element_justification="c", keep_on_top=False)
    prod_calculated = False
    help_text = "----------------------- FILE SELECTION ----------------- \n" \
                "From the main Menu use dropdown menu to select \n a TIEGCM-file." \
                "\n" \
                "\n" \
                "--------------------- ERRORS ------------------------- \n" \
                "Either choose the errors as percentage or click the \n button Real Error." \
                "\n" \
                "\n" \
                "----------------------- PROFILES -------------------------- \n" \
                "The vertical profile tab corresponds to height profile, \n " \
                "that is Latitude, Longitude and Timestep are \n " \
                "considered constant and the  height is in a range.\n" \
                " To choose Latitude, Longitude and Timestep use the \n " \
                "dropdown menus in the tab. Timestep values are not \n " \
                "time values, but correspond to TIEGCM file\u0027s \n " \
                "instances.\n " \
                "The map profile tab calculates values for a constant \n " \
                "height and Timestep  and for all the Latitude - \n" \
                " Longitude. Again the Timestep is not real time and \n " \
                "also neither the Pressure level, which corresponds \n " \
                "to different heights. \n " \
                "The nightshade option greys out the areas of the \n " \
                "map that are in night time." \
                "\n" \
                "\n" \
                "-------------------------- RUN --------------------------\n" \
                "To run the program after you have chosen a file, \n " \
                "the type of errors, the wanted profile and have \n " \
                "checked what plots are needed, then FIRST \n " \
                "click the button Calculate Products and after \n " \
                "that the button Calculate Error. For the plotting \n" \
                " process make sure the \"Choose Profile\" tab \n" \
                " and \"Choose Plots\" tab match before pressing the \n" \
                " calculate buttons.\n" \
                "The nightshade option greys out the areas of the \n " \
                "map that are in night time. \n" \
                "\n " \
                "-------------------- EXIT ------------------\n" \
                "To exit the program hit the Exit button."

    # #################### Event Loop ####################
    while True:
        event, values = window.read()
        if event in (Sg.WIN_CLOSED, "Exit"):
            break

        user_file_name = values["-FILE-"]

        if event == "Calculate Products" and values["-TABGROUP-"] == "Vertical Profile":
            user_lat = values["-LAT-"]
            user_lon = values["-LON-"]
            min_alt = values["-min_alt-"]
            max_alt = values["-max_alt-"]
            user_time = int(values["-TIME_vert-"])
            models_input(file_name=user_file_name, timer=user_time, lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon])
            Joule_Heating[:, :, :], Ohmic_Heating[:, :, :], Frictional_Heating[:, :, :], pedersen_con[:, :, :], hall_con[:, :, :], \
             parallel_con[:, :, :], nu_Op_sum[:, :, :], nu_O2p_sum[:, :, :], nu_NOp_sum[:, :, :], nu_e_sum[:, :, :], C_Op[:, :, :], C_O2p[:, :, :], \
             C_NOp[:, :, :], C_ion[:, :, :], J_pedersen[:, :, :], J_hall[:, :, :], J_ohmic[:, :, :], J_dens[:, :, :] = \
             calculate_products(Bx_in=Bx, By_in=By, Bz_in=Bz, Ex_in=Ex, Ey_in=Ey, Ez_in=Ez, Unx_in=Unx, Uny_in=Uny, Unz_in=Unz, Vix_in=Vi_vertx,
                                   Viy_in=Vi_verty, Viz_in=Vi_vertz, NO_in=NO, NO2_in=NO2, NN2_in=NN2, NOp_in=NOp, NO2p_in=NO2p, NNOp_in=NNOp,
                                   Ne_in=Ne, Ti_in=Ti, Te_in=Te, Tn_in=Tn, lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon])
            prod_calculated = True
            if values["-TABGROUP1-"] == "Vertical Profile Plots":
                if values["-COL-"]:
                    plot_collisions(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt, max_alt=max_alt)
                if values["-HR-"]:
                    plot_heating_rates(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt, max_alt=max_alt)
                if values["-CON-"]:
                    plot_conductivities(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt, max_alt=max_alt)
                if values["-CUR-"]:
                    plot_currents(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt, max_alt=max_alt)
                if values["-CR-"]:
                    plot_cross_sections(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt, max_alt=max_alt)
        if event == "Calculate Products" and values["-TABGROUP-"] == "Map Profile (Lat-Lon)":
            user_time = values["-TIME_map-"]
            user_lev = values["-Pr_level-"]
            models_input(file_name=user_file_name, timer=user_time, pressure_level=user_lev)
            Joule_Heating[:, :, :], Ohmic_Heating[:, :, :], Frictional_Heating[:, :, :], pedersen_con[:, :, :], hall_con[:, :, :], \
             parallel_con[:, :, :], nu_Op_sum[:, :, :], nu_O2p_sum[:, :, :], nu_NOp_sum[:, :, :], nu_e_sum[:, :, :], C_Op[:, :, :], C_O2p[:, :, :], \
             C_NOp[:, :, :], C_ion[:, :, :], J_pedersen[:, :, :], J_hall[:, :, :], J_ohmic[:, :, :], J_dens[:, :, :] = \
             calculate_products(Bx_in=Bx, By_in=By, Bz_in=Bz, Ex_in=Ex, Ey_in=Ey, Ez_in=Ez, Unx_in=Unx, Uny_in=Uny, Unz_in=Unz, Vix_in=Vi_vertx,
                                   Viy_in=Vi_verty, Viz_in=Vi_vertz, NO_in=NO, NO2_in=NO2, NN2_in=NN2, NOp_in=NOp, NO2p_in=NO2p, NNOp_in=NNOp,
                                   Ne_in=Ne, Ti_in=Ti, Te_in=Te, Tn_in=Tn, pressure_level=user_lev)
            prod_calculated = True
            if values["-TABGROUP1-"] == "Map Profile (Lat-Lon) Plots":
                user_lev = values["-Pr_level-"]
                night_shade = False
                if values["-NIGHT-"]:
                    night_shade = True
                if values["-HR_mapll-"]:
                    mapll_heating_rates_plot(pressure_level=user_lev, night_shade=night_shade)
                if values["-COL_mapll-"]:
                    mapll_collisions_plot(pressure_level=user_lev, night_shade=night_shade)
                if values["-CON_mapll-"]:
                    mapll_conductivities_plot(pressure_level=user_lev, night_shade=night_shade)
                if values["-CUR_mapll-"]:
                    mapll_currents_plot(pressure_level=user_lev, night_shade=night_shade)
                if values["-CR_mapll-"]:
                    mapll_csection_plot(pressure_level=user_lev, night_shade=night_shade)
        if event == "Calculate Products" and values["-TABGROUP-"] == "Map Profile (Lat-Alt)":
            user_lon = values["-Lon_map2-"]
            user_time = values["-TIME_map2-"]
            min_alt_la = values["-min_alt_la-"]
            max_alt_la = values["-max_alt_la-"]
            models_input(file_name=user_file_name, timer=user_time, lon_value=lon_dictionary[user_lon])
            Joule_Heating[:, :, :], Ohmic_Heating[:, :, :], Frictional_Heating[:, :, :], pedersen_con[:, :, :], hall_con[:, :, :], \
             parallel_con[:, :, :], nu_Op_sum[:, :, :], nu_O2p_sum[:, :, :], nu_NOp_sum[:, :, :], nu_e_sum[:, :, :], C_Op[:, :, :], C_O2p[:, :, :], \
             C_NOp[:, :, :], C_ion[:, :, :], J_pedersen[:, :, :], J_hall[:, :, :], J_ohmic[:, :, :], J_dens[:, :, :] = \
             calculate_products(Bx_in=Bx, By_in=By, Bz_in=Bz, Ex_in=Ex, Ey_in=Ey, Ez_in=Ez, Unx_in=Unx, Uny_in=Uny, Unz_in=Unz, Vix_in=Vi_vertx,
                                   Viy_in=Vi_verty, Viz_in=Vi_vertz, NO_in=NO, NO2_in=NO2, NN2_in=NN2, NOp_in=NOp, NO2p_in=NO2p, NNOp_in=NNOp,
                                   Ne_in=Ne, Ti_in=Ti, Te_in=Te, Tn_in=Tn, lon_value=lon_dictionary[user_lon])
            prod_calculated = True
            if values["-TABGROUP1-"] == "Map Profile (Lat-Alt) Plots":
                if values["-HR_mapla-"]:
                    mapla_heating_rates_plot(lon_value=lon_dictionary[user_lon], min_alt=min_alt_la, max_alt=max_alt_la)
                if values["-COL_mapla-"]:
                    mapla_collisions_plot(lon_value=lon_dictionary[user_lon], min_alt=min_alt_la, max_alt=max_alt_la)
                if values["-CON_mapla-"]:
                    mapla_conductivities_plot(lon_value=lon_dictionary[user_lon], min_alt=min_alt_la, max_alt=max_alt_la)
                if values["-CUR_mapla-"]:
                    mapla_currents_plot(lon_value=lon_dictionary[user_lon], min_alt=min_alt_la, max_alt=max_alt_la)
                if values["-CR_mapla-"]:
                    mapla_cross_section_plot(lon_value=lon_dictionary[user_lon], min_alt=min_alt_la, max_alt=max_alt_la)
        if event == "Calculate Error" and values["-TABGROUP-"] == 'Vertical Profile':
            if prod_calculated:
                user_lat = values["-LAT-"]
                user_lon = values["-LON-"]
                calculate_noise_2(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon])
                Joule_Heating_noisy[:, :, :], Ohmic_Heating_noisy[:, :, :], Frictional_Heating_noisy[:, :, :], pedersen_con_noisy[:, :, :], \
                 hall_con_noisy[:, :, :], parallel_con_noisy[:, :, :], nu_Op_sum_noisy[:, :, :], nu_O2p_sum_noisy[:, :, :], nu_NOp_sum_noisy[:, :, :],\
                 nu_e_sum_noisy[:, :, :], C_Op_noisy[:, :, :], C_O2p_noisy[:, :, :], C_NOp_noisy[:, :, :], C_ion_noisy[:, :, :], \
                 J_pedersen_noisy[:, :, :], J_hall_noisy[:, :, :], J_ohmic_noisy[:, :, :], J_dens_noisy[:, :, :] = \
                 calculate_products(Bx_in=Bx_noisy, By_in=By_noisy, Bz_in=Bz_noisy, Ex_in=Ex_noisy, Ey_in=Ey_noisy, Ez_in=Ez_noisy,
                                       Unx_in=Unx_noisy, Uny_in=Uny_noisy, Unz_in=Unz_noisy, Vix_in=Vi_vertx_noisy, Viy_in=Vi_verty_noisy,
                                       Viz_in=Vi_vertz_noisy, NO_in=NO_noisy, NO2_in=NO2_noisy, NN2_in=NN2_noisy, NOp_in=NOp_noisy,
                                       NO2p_in=NO2p_noisy, NNOp_in=NNOp_noisy, Ne_in=Ne_noisy, Ti_in=Ti_noisy, Te_in=Te_noisy, Tn_in=Tn_noisy,
                                       lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon])
                if values["-TABGROUP1-"] == "Vertical Profile Plots":
                    min_alt = values["-min_alt-"]
                    max_alt = values["-max_alt-"]
                    if values["-COL_abs-"]:
                        plot_collisions_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                              max_alt=max_alt)
                    if values["-COL_rel-"]:
                        plot_collisions_rel_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                                  max_alt=max_alt)
                    if values["-HR_abs-"]:
                        plot_heating_rates_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                                 max_alt=max_alt)
                    if values["-HR_rel-"]:
                        plot_heating_rates_rel_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                                     max_alt=max_alt)
                    if values["-CON_abs-"]:
                        plot_conductivities_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                                  max_alt=max_alt)
                    if values["-CON_rel-"]:
                        plot_conductivities_rel_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                                      max_alt=max_alt)
                    if values["-CUR_abs-"]:
                        plot_currents_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt, max_alt=max_alt)
                    if values["-CUR_rel-"]:
                        plot_currents_rel_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                                max_alt=max_alt)
                    if values["-CR_abs-"]:
                        plot_csections_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt, max_alt=max_alt)
                    if values["-CR_rel-"]:
                        plot_csections_rel_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                                 max_alt=max_alt)
            else:
                Sg.popup("First Products Must Be Calculated", title="Input Error", keep_on_top=True)
        if event == "Calculate Error" and values["-TABGROUP-"] == "Map Profile (Lat-Lon)":
            user_lev = values["-Pr_level-"]
            if prod_calculated:
                calculate_noise_2(pressure_level=user_lev)
                Joule_Heating_noisy[:, :, :], Ohmic_Heating_noisy[:, :, :], Frictional_Heating_noisy[:, :, :], pedersen_con_noisy[:, :, :], \
                 hall_con_noisy[:, :, :], parallel_con_noisy[:, :, :], nu_Op_sum_noisy[:, :, :], nu_O2p_sum_noisy[:, :, :], nu_NOp_sum_noisy[:, :, :], \
                 nu_e_sum_noisy[:, :, :], C_Op_noisy[:, :, :], C_O2p_noisy[:, :, :], C_NOp_noisy[:, :, :], C_ion_noisy[:, :, :], \
                 J_pedersen_noisy[:, :, :], J_hall_noisy[:, :, :], J_ohmic_noisy[:, :, :], J_dens_noisy[:, :, :] = \
                 calculate_products(Bx_in=Bx_noisy, By_in=By_noisy, Bz_in=Bz_noisy, Ex_in=Ex_noisy, Ey_in=Ey_noisy, Ez_in=Ez_noisy,
                                       Unx_in=Unx_noisy, Uny_in=Uny_noisy, Unz_in=Unz_noisy, Vix_in=Vi_vertx_noisy, Viy_in=Vi_verty_noisy,
                                       Viz_in=Vi_vertz_noisy, NO_in=NO_noisy, NO2_in=NO2_noisy, NN2_in=NN2_noisy, NOp_in=NOp_noisy,
                                       NO2p_in=NO2p_noisy, NNOp_in=NNOp_noisy, Ne_in=Ne_noisy, Ti_in=Ti_noisy, Te_in=Te_noisy, Tn_in=Tn_noisy,
                                       pressure_level=user_lev)
                if values["-TABGROUP1-"] == "Map Profile (Lat-Lon) Plots":
                    user_lev = values["-Pr_level-"]
                    night_shade = False
                    if values["-NIGHT-"]:
                        night_shade = True
                    if values["-HR_mapll_error-"]:
                        mapll_heating_rates_rel_error_plot(pressure_level=user_lev, night_shade=night_shade)
                    if values["-COL_mapll_error-"]:
                        mapll_collisions_rel_error_plot(pressure_level=user_lev, night_shade=night_shade)
                    if values["-CON_mapll_error-"]:
                        mapll_conductivities_rel_error_plot(pressure_level=user_lev, night_shade=night_shade)
                    if values["-CUR_mapll_error-"]:
                        mapll_currents_rel_error_plot(pressure_level=user_lev, night_shade=night_shade)
                    if values["-CR_mapll_error-"]:
                        mapll_csection_rel_error_plot(pressure_level=user_lev, night_shade=night_shade)
            else:
                Sg.popup("First Products Must Be Calculated", title="Input Error", keep_on_top=True)
        if event == "Calculate Error" and values["-TABGROUP-"] == "Map Profile (Lat-Alt)":
            user_lon = values["-Lon_map2-"]
            min_alt_la = values["-min_alt_la-"]
            max_alt_la = values["-max_alt_la-"]
            if prod_calculated:
                calculate_noise_2(lon_value=lon_dictionary[user_lon])
                Joule_Heating_noisy[:, :, :], Ohmic_Heating_noisy[:, :, :], Frictional_Heating_noisy[:, :, :], pedersen_con_noisy[:, :, :], \
                 hall_con_noisy[:, :, :], parallel_con_noisy[:, :, :], nu_Op_sum_noisy[:, :, :], nu_O2p_sum_noisy[:, :, :], nu_NOp_sum_noisy[:, :, :], \
                 nu_e_sum_noisy[:, :, :], C_Op_noisy[:, :, :], C_O2p_noisy[:, :, :], C_NOp_noisy[:, :, :], C_ion_noisy[:, :, :], \
                 J_pedersen_noisy[:, :, :], J_hall_noisy[:, :, :], J_ohmic_noisy[:, :, :], J_dens_noisy[:, :, :] = \
                 calculate_products(Bx_in=Bx_noisy, By_in=By_noisy, Bz_in=Bz_noisy, Ex_in=Ex_noisy, Ey_in=Ey_noisy, Ez_in=Ez_noisy,
                                       Unx_in=Unx_noisy, Uny_in=Uny_noisy, Unz_in=Unz_noisy, Vix_in=Vi_vertx_noisy, Viy_in=Vi_verty_noisy,
                                       Viz_in=Vi_vertz_noisy, NO_in=NO_noisy, NO2_in=NO2_noisy, NN2_in=NN2_noisy, NOp_in=NOp_noisy,
                                       NO2p_in=NO2p_noisy, NNOp_in=NNOp_noisy, Ne_in=Ne_noisy, Ti_in=Ti_noisy, Te_in=Te_noisy, Tn_in=Tn_noisy,
                                       lon_value=lon_dictionary[user_lon])
                if values["-TABGROUP1-"] == "Map Profile (Lat-Alt) Plots":
                    if values["-HR_mapla_error-"]:
                        mapla_heating_rates_rel_error_plot(lon_value=lon_dictionary[user_lon], min_alt=min_alt_la, max_alt=max_alt_la)
                    if values["-COL_mapla_error-"]:
                        mapla_collisions_rel_error_plot(lon_dictionary[user_lon], min_alt=min_alt_la, max_alt=max_alt_la)
                    if values["-CON_mapla_error-"]:
                        mapla_conductivities_rel_error_plot(lon_dictionary[user_lon], min_alt=min_alt_la, max_alt=max_alt_la)
                    if values["-CUR_mapla_error-"]:
                        mapla_currents_rel_error_plot(lon_dictionary[user_lon], min_alt=min_alt_la, max_alt=max_alt_la)
                    if values["-CR_mapla_error-"]:
                        mapla_cross_section_rel_error_plot(lon_dictionary[user_lon], min_alt=min_alt_la, max_alt=max_alt_la)
            else:
                Sg.popup("First Products Must Be Calculated", title="Input Error", keep_on_top=True)
        if event == "Help":
            Sg.popup(help_text, title="Help", keep_on_top=True)
    # ######################### Close Window #########################
    window.close()


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$ GUI CREATION END $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# Run gui
gui()
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END OF PROGRAM @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
