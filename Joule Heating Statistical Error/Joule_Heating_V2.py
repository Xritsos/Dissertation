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


def insert_noise(accuracy):

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

    # Create loading bar
    layout = [[Sg.Text('Getting data from models........')], [Sg.ProgressBar(210, orientation='h', size=(50, 20), key='progressbar')]]

    window = Sg.Window('Progress', keep_on_top=True).Layout(layout)
    progress_bar = window.FindElement('progressbar')
    event, values = window.Read(timeout=1)

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
    progress_bar.UpdateBar(0, 30)
    glev_in = tiegcm.variables['lev'][:]     # midpoint levels
    time_in = tiegcm.variables['time'][:]    # minutes since 2015-1-1 0:0:0
    O_in = tiegcm.variables['O_CM3'][:]      # atomic oxygen density (neutral) in cm^(-3)
    O2_in = tiegcm.variables['O2_CM3'][:]    # molecular oxygen density (neutral) in cm^(-3)
    N2_in = tiegcm.variables['N2_CM3'][:]    # molecular nitrogen density (neutral) in cm^(-3)
    progress_bar.UpdateBar(30, 80)
    Op_in = tiegcm.variables['OP'][:]        # atomic oxygen density (ion) in cm^(-3)
    O2p_in = tiegcm.variables['O2P'][:]      # molecular oxygen density (ion) in cm^(-3)
    NOp_in = tiegcm.variables['NOP_LAM'][:]  # nitric oxide density (ion) in cm^(-3)
    Te_in = tiegcm.variables['TE'][:]        # electron temperature in kelvin
    Tn_in = tiegcm.variables['TN'][:]        # neutral temperature in kelvin
    Ti_in = tiegcm.variables['TI'][:]        # ion temperature in kelvin
    progress_bar.UpdateBar(80, 140)
    Un_north_in = tiegcm.variables['VN_si'][:]  # neutral meridional wind (+north) in m/sec
    Un_east_in = tiegcm.variables['UN_si'][:]   # neutral zonal wind (+east) in m/sec
    Un_up_in = tiegcm.variables['WN_si'][:]     # neutral vertical wind (+up) in m/sec

    Ee_in = tiegcm.variables['EEX_si'][:]  # electric field (+east) in V/m
    En_in = tiegcm.variables['EEY_si'][:]  # electric field (+north) in V/m
    Eu_in = tiegcm.variables['EEZ_si'][:]  # electric field (+up) in V/m
    progress_bar.UpdateBar(140, 180)

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

                nu_Op_sum = nu_Op_N2 + nu_Op_O2 + nu_Op_O

                # nu-O2p = nu(O2p-N2) + nu(O2p-O) + nu(O2p-O2) (in Hz)
                # densities in cm^(-3)
                nu_O2p_N2 = 4.13 * NN2[lat, lon, lev] * 10 ** (-10)
                nu_O2p_O = 2.31 * NO[lat, lon, lev] * 10 ** (-10)
                nu_O2p_O2 = 2.59 * NO2[lat, lon, lev] * 10 ** (-11) * Tr ** (1 / 2) * (1 - 0.073 * np.log10(Tr)) ** 2

                nu_O2p_sum = nu_O2p_N2 + nu_O2p_O + nu_O2p_O2

                # nu-NOp = nu(NOp-N2) + nu(NOp-O) + nu(NOp-O2) (in Hz)
                # densities in cm^(-3)
                nu_NOp_N2 = 4.34 * NN2[lat, lon, lev] * 10 ** (-10)
                nu_NOp_O = 2.44 * NO[lat, lon, lev] * 10 ** (-10)
                nu_NOp_O2 = 4.27 * NO2[lat, lon, lev] * 10 ** (-10)

                nu_NOp_sum = nu_NOp_N2 + nu_NOp_O + nu_NOp_O2

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
                Vi_Op_star_x = (nu_Op_sum * omega_Op * Estar[0] + omega_Op ** 2 * EstarXbunit[0]) / \
                               (Bnorm * (nu_Op_sum ** 2 + omega_Op ** 2))
                Vi_Op_star_y = (nu_Op_sum * omega_Op * Estar[1] + omega_Op ** 2 * EstarXbunit[1]) / \
                               (Bnorm * (nu_Op_sum ** 2 + omega_Op ** 2))
                Vi_Op_star_z = (nu_Op_sum * omega_Op * Estar[2] + omega_Op ** 2 * EstarXbunit[2]) / \
                               (Bnorm * (nu_Op_sum ** 2 + omega_Op ** 2))
                # ############# O2+ #############
                Vi_O2p_star_x = (nu_O2p_sum * omega_O2p * Estar[0] + omega_O2p ** 2 * EstarXbunit[0]) / \
                                (Bnorm * (nu_O2p_sum ** 2 + omega_O2p ** 2))
                Vi_O2p_star_y = (nu_O2p_sum * omega_O2p * Estar[1] + omega_O2p ** 2 * EstarXbunit[1]) / \
                                (Bnorm * (nu_O2p_sum ** 2 + omega_O2p ** 2))
                Vi_O2p_star_z = (nu_O2p_sum * omega_O2p * Estar[2] + omega_O2p ** 2 * EstarXbunit[2]) / \
                                (Bnorm * (nu_O2p_sum ** 2 + omega_O2p ** 2))
                # ############ NO+ ############
                Vi_NOp_star_x = (nu_NOp_sum * omega_NOp * Estar[0] + omega_NOp ** 2 * EstarXbunit[0]) / \
                                (Bnorm * (nu_NOp_sum ** 2 + omega_NOp ** 2))
                Vi_NOp_star_y = (nu_NOp_sum * omega_NOp * Estar[1] + omega_NOp ** 2 * EstarXbunit[1]) / \
                                (Bnorm * (nu_NOp_sum ** 2 + omega_NOp ** 2))
                Vi_NOp_star_z = (nu_NOp_sum * omega_NOp * Estar[2] + omega_NOp ** 2 * EstarXbunit[2]) / \
                                (Bnorm * (nu_NOp_sum ** 2 + omega_NOp ** 2))

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
                print(Vi_vertx[lat, lon, lev], Vi_verty[lat, lon, lev], Vi_vertz[lat, lon, lev])

    progress_bar.UpdateBar(180, 210)
    time.sleep(3)
    window.close()
    Sg.popup("_"*40, "Data imported in : " + str(time.time() - start_time) + " sec!", "_"*40, title="Finished", keep_on_top=True)
    print('Data imported in: ', str(time.time() - start_time), ' sec!')
    print(' ')
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END OF INPUT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ GUI CREATION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
    col1 = [[Sg.Checkbox("Plot Heating Rates Absolute Error", default=False, tooltip="Plots heating rates absolute error", key="-HR_abs-")],
            [Sg.Checkbox("Plot Heating Rates Relative Error", default=False, tooltip="Plots heating rates relative error", key="-HR_rel-")],
            [Sg.Checkbox("Plot Collision Frequencies Absolute Error", default=False, tooltip="Plots collision frequencies absolute error",
                key="-COL_abs-")],
            [Sg.Checkbox("Plot Collision Frequencies Relative Error", default=False, tooltip="Plots collision frequencies relative error",
                key="-COL_rel-")],
            [Sg.Checkbox("Plot Conductivities Absolute Error", default=False, tooltip="Plots conductivities absolute error", key="-CON_abs-")],
            [Sg.Checkbox("Plot Conductivities Relative Error", default=False, tooltip="Plots conductivities relative error", key="-CON_rel-")],
            [Sg.Checkbox("Plot Currents Absolute Error", default=False, tooltip="Plots currents absolute error", key="-CUR_abs-")],
            [Sg.Checkbox("Plot Currents Relative Error", default=False, tooltip="Plots currents relative error", key="-CUR_rel-")],
            [Sg.Checkbox("Plot Cross Sections Absolute Error", default=False, tooltip="Plots cross sections absolute error", key="-CR_abs-")],
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

        if event=="Calculate Products" and values["-TABGROUP-"]=="Vertical Profile":
            user_lat = values["-LAT-"]
            user_lon = values["-LON-"]
            min_alt = values["-min_alt-"]
            max_alt = values["-max_alt-"]
            user_time = int(values["-TIME_vert-"])
            models_input(file_name=user_file_name, timer=user_time, lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon])
            prod_calculated = True
            if values["-TABGROUP1-"]=="Vertical Profile Plots":
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
        if event=="Calculate Products" and values["-TABGROUP-"]=="Map Profile (Lat-Lon)":
            user_time = values["-TIME_map-"]
            user_lev = values["-Pr_level-"]
            models_input(file_name=user_file_name, timer=user_time, pressure_level=user_lev)
            prod_calculated = True
            if values["-TABGROUP1-"]=="Map Profile (Lat-Lon) Plots":
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
        if event=="Calculate Products" and values["-TABGROUP-"]=="Map Profile (Lat-Alt)":
            user_lon = values["-Lon_map2-"]
            user_time = values["-TIME_map2-"]
            min_alt_la = values["-min_alt_la-"]
            max_alt_la = values["-max_alt_la-"]
            models_input(file_name=user_file_name, timer=user_time, lon_value=lon_dictionary[user_lon])
            prod_calculated = True
            if values["-TABGROUP1-"]=="Map Profile (Lat-Alt) Plots":
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
        if event=="Calculate Error" and values["-TABGROUP-"]=='Vertical Profile':
            if prod_calculated:
                user_lat = values["-LAT-"]
                user_lon = values["-LON-"]
                if values["-TABGROUP1-"]=="Vertical Profile Plots":
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
        if event=="Calculate Error" and values["-TABGROUP-"]=="Map Profile (Lat-Lon)":
            user_lev = values["-Pr_level-"]
            if values["-TABGROUP1-"]=="Map Profile (Lat-Lon) Plots":
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
        if event=="Calculate Error" and values["-TABGROUP-"]=="Map Profile (Lat-Alt)":
            user_lon = values["-Lon_map2-"]
            min_alt_la = values["-min_alt_la-"]
            max_alt_la = values["-max_alt_la-"]
            if values["-TABGROUP1-"]=="Map Profile (Lat-Alt) Plots":
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
        if event=="Help":
            Sg.popup(help_text, title="Help", keep_on_top=True)
    # ######################### Close Window #########################
    window.close()


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$ GUI CREATION END $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# Run gui
gui()
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END OF PROGRAM @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@