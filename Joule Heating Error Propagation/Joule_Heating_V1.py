"""
Joule Heating Error Propagation Code
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

# Electric field
Ex = np.zeros((72, 144, 57), order='F')
Ey = np.zeros((72, 144, 57), order='F')
Ez = np.zeros((72, 144, 57), order='F')

# Neutral wind
Unx = np.zeros((72, 144, 57), order='F')
Uny = np.zeros((72, 144, 57), order='F')
Unz = np.zeros((72, 144, 57), order='F')

# Densities
NO = np.zeros((72, 144, 57), order='F')
NO2 = np.zeros((72, 144, 57), order='F')
NN2 = np.zeros((72, 144, 57), order='F')
NOp = np.zeros((72, 144, 57), order='F')
NO2p = np.zeros((72, 144, 57), order='F')
NNOp = np.zeros((72, 144, 57), order='F')
Ne = np.zeros((72, 144, 57), order='F')

# Temperatures
Te = np.zeros((72, 144, 57), order='F')
Ti = np.zeros((72, 144, 57), order='F')
Tn = np.zeros((72, 144, 57), order='F')

# Collision frequencies
nu_Op_sum = np.zeros((72, 144, 57), order='F')
nu_O2p_sum = np.zeros((72, 144, 57), order='F')
nu_NOp_sum = np.zeros((72, 144, 57), order='F')
nu_e_sum = np.zeros((72, 144, 57), order='F')

# Conductivities
pedersen_con = np.zeros((72, 144, 57), order='F')
hall_con = np.zeros((72, 144, 57), order='F')
parallel_con = np.zeros((72, 144, 57), order='F')

# Ion velocities perpendicular to magnetic field
Vi_vertx = np.zeros((72, 144, 57), order='F')
Vi_verty = np.zeros((72, 144, 57), order='F')
Vi_vertz = np.zeros((72, 144, 57), order='F')

# Heating rates
Joule_Heating = np.zeros((72, 144, 57), order='F')
Frictional_Heating = np.zeros((72, 144, 57), order='F')
Ohmic_Heating = np.zeros((72, 144, 57), order='F')

# Cross sections
C_Op = np.zeros((72, 144, 57), order='F')
C_O2p = np.zeros((72, 144, 57), order='F')
C_NOp = np.zeros((72, 144, 57), order='F')
C_ion = np.zeros((72, 144, 57), order='F')

# Perpendicular currents
J_pedersen = np.zeros((72, 144, 57), order='F')
J_hall = np.zeros((72, 144, 57), order='F')
J_ohmic = np.zeros((72, 144, 57), order='F')
J_dens = np.zeros((72, 144, 57), order='F')

# Collision frequencies error
nuOp_error = np.zeros((72, 144, 57), order='F')
nuO2p_error = np.zeros((72, 144, 57), order='F')
nuNOp_error = np.zeros((72, 144, 57), order='F')
nue_error = np.zeros((72, 144, 57), order='F')

# Gyro-frequencies
Omega_ion = np.zeros((72, 144, 57), order='F')
Omega_e = np.zeros((72, 144, 57), order='F')

# Conductivities error
pedersen_con_error = np.zeros((72, 144, 57), order='F')
hall_con_error = np.zeros((72, 144, 57), order='F')
parallel_con_error = np.zeros((72, 144, 57), order='F')

# Heating rates error
Joule_Heating_error = np.zeros((72, 144, 57), order='F')
Frictional_Heating_error = np.zeros((72, 144, 57), order='F')
Ohmic_Heating_error = np.zeros((72, 144, 57), order='F')

# Cross sections error
C_Op_error = np.zeros((72, 144, 57), order='F')
C_O2p_error = np.zeros((72, 144, 57), order='F')
C_NOp_error = np.zeros((72, 144, 57), order='F')
C_ion_error = np.zeros((72, 144, 57), order='F')

# Perpendicular currents
J_pedersen_error = np.zeros((72, 144, 57), order='F')
J_hall_error = np.zeros((72, 144, 57), order='F')
J_ohmic_error = np.zeros((72, 144, 57), order='F')
J_dens_error = np.zeros((72, 144, 57), order='F')

# Collision frequencies contributions error
dnuion_Ti = np.zeros((72, 144, 57), order='F')
dnuion_Tn = np.zeros((72, 144, 57), order='F')
dnuion_Nneutral = np.zeros((72, 144, 57), order='F')
dnue_Te = np.zeros((72, 144, 57), order='F')
dnue_Nneutral = np.zeros((72, 144, 57), order='F')

# Pedersen conductivity contributions error
dsp_B = np.zeros((72, 144, 57), order='F')
dsp_Te = np.zeros((72, 144, 57), order='F')
dsp_Ti = np.zeros((72, 144, 57), order='F')
dsp_Tn = np.zeros((72, 144, 57), order='F')
dsp_Nion = np.zeros((72, 144, 57), order='F')
dsp_Ne = np.zeros((72, 144, 57), order='F')
dsp_Nneutral = np.zeros((72, 144, 57), order='F')

# Hall conductivity contributions error
dsh_B = np.zeros((72, 144, 57), order='F')
dsh_Te = np.zeros((72, 144, 57), order='F')
dsh_Ti = np.zeros((72, 144, 57), order='F')
dsh_Tn = np.zeros((72, 144, 57), order='F')
dsh_Nion = np.zeros((72, 144, 57), order='F')
dsh_Ne = np.zeros((72, 144, 57), order='F')
dsh_Nneutral = np.zeros((72, 144, 57), order='F')

# Joule Heating contributions error
dJH_B = np.zeros((72, 144, 57), order='F')
dJH_E = np.zeros((72, 144, 57), order='F')
dJH_Vi = np.zeros((72, 144, 57), order='F')
dJH_Un = np.zeros((72, 144, 57), order='F')
dJH_Ne = np.zeros((72, 144, 57), order='F')

# Ohmic Heating contributions error
dOH_B = np.zeros((72, 144, 57), order='F')
dOH_E = np.zeros((72, 144, 57), order='F')
dOH_Nneutral = np.zeros((72, 144, 57), order='F')
dOH_Nion = np.zeros((72, 144, 57), order='F')
dOH_Un = np.zeros((72, 144, 57), order='F')
dOH_Ne = np.zeros((72, 144, 57), order='F')
dOH_Te = np.zeros((72, 144, 57), order='F')
dOH_Tn = np.zeros((72, 144, 57), order='F')
dOH_Ti = np.zeros((72, 144, 57), order='F')
dOH_sp = np.zeros((72, 144, 57), order='F')

# Frictional Heating contributions error
dFH_B = np.zeros((72, 144, 57), order='F')
dFH_Nneutral = np.zeros((72, 144, 57), order='F')
dFH_Nion = np.zeros((72, 144, 57), order='F')
dFH_Un = np.zeros((72, 144, 57), order='F')
dFH_Tn = np.zeros((72, 144, 57), order='F')
dFH_Ti = np.zeros((72, 144, 57), order='F')
dFH_Vi = np.zeros((72, 144, 57), order='F')
dFH_nu = np.zeros((72, 144, 57), order='F')

# Ion cross section contribution error
dCion_Ti = np.zeros((72, 144, 57), order='F')
dCion_Tn = np.zeros((72, 144, 57), order='F')
dCion_nu = np.zeros((72, 144, 57), order='F')
dCion_Nneutral = np.zeros((72, 144, 57), order='F')

# Perpendicular currents contributions error
dJohm_B = np.zeros((72, 144, 57), order='F')
dJohm_E = np.zeros((72, 144, 57), order='F')
dJohm_Un = np.zeros((72, 144, 57), order='F')
dJohm_sp = np.zeros((72, 144, 57), order='F')
dJohm_sh = np.zeros((72, 144, 57), order='F')
dJohm_Ti = np.zeros((72, 144, 57), order='F')
dJohm_Tn = np.zeros((72, 144, 57), order='F')
dJohm_Te = np.zeros((72, 144, 57), order='F')
dJohm_Ne = np.zeros((72, 144, 57), order='F')
dJohm_Nneutral = np.zeros((72, 144, 57), order='F')
dJohm_Nion = np.zeros((72, 144, 57), order='F')

dJd_B = np.zeros((72, 144, 57), order='F')
dJd_E = np.zeros((72, 144, 57), order='F')
dJd_Vi = np.zeros((72, 144, 57), order='F')
dJd_Un = np.zeros((72, 144, 57), order='F')
dJd_Ne = np.zeros((72, 144, 57), order='F')


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

    progress_bar.UpdateBar(180, 210)
    time.sleep(3)
    window.close()
    Sg.popup("_"*40, "Data imported in : " + str(time.time() - start_time) + " sec!", "_"*40, title="Finished", keep_on_top=True)
    print('Data imported in: ', str(time.time() - start_time), ' sec!')
    print(' ')
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END OF INPUT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ CALCULATE PRODUCTS USING MODEL INPUTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def products(lat_value=-1, lon_value=-1, pressure_level=-1):
    start_time = time.time()
    print('Calculating Products.....')
    print(' ')
    print(title)
    # distinguish Map from Vertical profile
    lev_range = 0
    lat_range = 0
    lon_range = 0

    lev_start = 0
    lat_start = 0
    lon_start = 0
    max_bar = 0

    # Lat - Alt map profile
    if lat_value == -1 and pressure_level == -1:
        lev_range = len(glev_in) - 1
        lat_range = len(glat_in)
        lon_start = lon_value
        lon_range = lon_start + 1
        max_bar = lev_range

    # Lat - Lon map profile
    if lat_value == -1 and lon_value == -1:
        lev_start = pressure_level
        lev_range = lev_start + 1
        lat_range = len(glat_in)
        lon_range = len(glon_in)
        max_bar = lat_range

    # Vertical profile
    if lat_value != -1 and lon_value != -1:
        lat_start = lat_value
        lon_start = lon_value
        lat_range = lat_start + 1
        lon_range = lon_start + 1
        lev_range = len(glev_in) - 1
        max_bar = lev_range

    # Create loading bar
    layout = [[Sg.Text('Calculating Products.....')], [Sg.ProgressBar(max_bar, orientation='h', size=(50, 20), key='progressbar')]]

    window = Sg.Window(title=title, keep_on_top=True).Layout(layout)
    progress_bar = window.FindElement('progressbar')
    event, values = window.Read(timeout=1)

    for lev in range(lev_start, lev_range):
        for lat in range(lat_start, lat_range):
            for lon in range(lon_start, lon_range):
                # progress bar index
                if max_bar == lev_range:
                    i = lev
                else:
                    i = lat
                # ########################## COLLISION FREQUENCIES ##########################
                # nu-Op = nu(Op-N2) + nu(Op-O2) + nu(Op-O) (in Hz)
                # densities in cm^(-3)
                nu_Op_N2 = 6.82 * NN2[lat, lon, lev] * 10 ** (-10)
                nu_Op_O2 = 6.64 * NO2[lat, lon, lev] * 10 ** (-10)

                Tr = (Ti[lat, lon, lev] + Tn[lat, lon, lev]) / 2  # in kelvin
                nu_Op_O = fb * (3.67 * NO[lat, lon, lev] * 10 ** (-11) * Tr ** (1 / 2) * (1 - 0.064 * np.log10(Tr)) ** 2)

                nu_Op_sum[lat, lon, lev] = nu_Op_N2 + nu_Op_O2 + nu_Op_O

                # nu-O2p = nu(O2p-N2) + nu(O2p-O) + nu(O2p-O2) (in Hz)
                # densities in cm^(-3)
                nu_O2p_N2 = 4.13 * NN2[lat, lon, lev] * 10 ** (-10)
                nu_O2p_O = 2.31 * NO[lat, lon, lev] * 10 ** (-10)
                nu_O2p_O2 = 2.59 * NO2[lat, lon, lev] * 10 ** (-11) * Tr ** (1 / 2) * (1 - 0.073 * np.log10(Tr)) ** 2

                nu_O2p_sum[lat, lon, lev] = nu_O2p_N2 + nu_O2p_O + nu_O2p_O2

                # nu-NOp = nu(NOp-N2) + nu(NOp-O) + nu(NOp-O2) (in Hz)
                # densities in cm^(-3)
                nu_NOp_N2 = 4.34 * NN2[lat, lon, lev] * 10 ** (-10)
                nu_NOp_O = 2.44 * NO[lat, lon, lev] * 10 ** (-10)
                nu_NOp_O2 = 4.27 * NO2[lat, lon, lev] * 10 ** (-10)

                nu_NOp_sum[lat, lon, lev] = nu_NOp_N2 + nu_NOp_O + nu_NOp_O2

                # nu-e = nu(e-N2) + nu(e-O) + nu(e-O2) (in Hz)
                # densities in cm^(-3)
                nu_e_N2 = 2.33 * 10 ** (-11) * NN2[lat, lon, lev] * Te[lat, lon, lev] * (1 - 1.21 * 10 ** (-4) * Te[lat, lon, lev])
                nu_e_O2 = 1.82 * 10 ** (-10) * NO2[lat, lon, lev] * Te[lat, lon, lev] ** (1 / 2) * (1 +
                                                                                                    3.6 * 10 ** (-2) * Te[lat, lon, lev] ** (1 / 2))
                nu_e_O = 8.9 * 10 ** (-11) * NO[lat, lon, lev] * Te[lat, lon, lev] ** (1 / 2) * (1 + 5.7 * 10 ** (-4) * Te[lat, lon, lev])

                nu_e_sum[lat, lon, lev] = nu_e_N2 + nu_e_O2 + nu_e_O
                # ################ GYRO-FREQUENCIES(OMEGAS) ################
                # Magnetic field vector (in tesla)
                B = [Bx[lat, lon, lev], By[lat, lon, lev], Bz[lat, lon, lev]]
                # Magnetic field magnitude (in tesla)
                Bnorm = np.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)
                # Magnetic field unit vector
                bunit = [B[0] / Bnorm, B[1] / Bnorm, B[2] / Bnorm]

                # qe(in coulomb), mk(masses in kg), omegas(in Hz)
                omega_Op = (qe * Bnorm) / mkO
                omega_O2p = (qe * Bnorm) / mkO2
                omega_NOp = (qe * Bnorm) / mkNO
                omega_e = (qe * Bnorm) / me

                Omega_ion[lat, lon, lev] = (omega_Op + omega_O2p + omega_NOp) / 3
                Omega_e[lat, lon, lev] = omega_e
                # ################## RATIOS ##################
                # dimensionless
                r_Op = nu_Op_sum[lat, lon, lev] / omega_Op
                r_O2p = nu_O2p_sum[lat, lon, lev] / omega_O2p
                r_NOp = nu_NOp_sum[lat, lon, lev] / omega_NOp
                r_e = nu_e_sum[lat, lon, lev] / omega_e
                # ############################# CONDUCTIVITIES #############################
                # Pedersen conductivity (in siemens/meter)
                # qe(in coulomb), B_norm(in tesla), N(densities in m^(-3)), ratios(dimensionless)
                term_a_ped = (Ne[lat, lon, lev] * ccm) * (r_e / (1 + r_e ** 2))
                term_b_ped = (NOp[lat, lon, lev] * ccm) * (r_Op / (1 + r_Op ** 2))
                term_c_ped = (NO2p[lat, lon, lev] * ccm) * (r_O2p / (1 + r_O2p ** 2))
                term_d_ped = (NNOp[lat, lon, lev] * ccm) * (r_NOp / (1 + r_NOp ** 2))

                pedersen_con[lat, lon, lev] = (qe / Bnorm) * (term_a_ped + term_b_ped + term_c_ped + term_d_ped)

                # Hall conductivity (in siemens/meter)
                # qe(in coulomb), B_norm(in tesla), N(densities in m^(-3)), ratios(dimensionless)
                term_a_hall = (Ne[lat, lon, lev] * ccm) / (1 + r_e ** 2)
                term_b_hall = (NOp[lat, lon, lev] * ccm) / (1 + r_Op ** 2)
                term_c_hall = (NO2p[lat, lon, lev] * ccm) / (1 + r_O2p ** 2)
                term_d_hall = (NNOp[lat, lon, lev] * ccm) / (1 + r_NOp ** 2)

                hall_con[lat, lon, lev] = (qe / Bnorm) * (term_a_hall - term_b_hall - term_c_hall - term_d_hall)

                # Parallel conductivity (in siemens/meter)
                # qe(in coulomb), me(mass in tesla), N(density) (in m^(-3)), collision frequency(in Hz)
                parallel_con[lat, lon, lev] = (Ne[lat, lon, lev] * ccm * qe ** 2) / (me * nu_e_sum[lat, lon, lev])

                # ################################ HEATING RATES ################################
                # Electric field vector(in volt/meter)
                E = [Ex[lat, lon, lev], Ey[lat, lon, lev], Ez[lat, lon, lev]]

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
                Vi_Op_star_x = (nu_Op_sum[lat, lon, lev] * omega_Op * Estar[0] + omega_Op ** 2 * EstarXbunit[0]) / \
                               (Bnorm * (nu_Op_sum[lat, lon, lev] ** 2 + omega_Op ** 2))
                Vi_Op_star_y = (nu_Op_sum[lat, lon, lev] * omega_Op * Estar[1] + omega_Op ** 2 * EstarXbunit[1]) / \
                               (Bnorm * (nu_Op_sum[lat, lon, lev] ** 2 + omega_Op ** 2))
                Vi_Op_star_z = (nu_Op_sum[lat, lon, lev] * omega_Op * Estar[2] + omega_Op ** 2 * EstarXbunit[2]) / \
                               (Bnorm * (nu_Op_sum[lat, lon, lev] ** 2 + omega_Op ** 2))
                # ############# O2+ #############
                Vi_O2p_star_x = (nu_O2p_sum[lat, lon, lev] * omega_O2p * Estar[0] + omega_O2p ** 2 * EstarXbunit[0]) / \
                                (Bnorm * (nu_O2p_sum[lat, lon, lev] ** 2 + omega_O2p ** 2))
                Vi_O2p_star_y = (nu_O2p_sum[lat, lon, lev] * omega_O2p * Estar[1] + omega_O2p ** 2 * EstarXbunit[1]) / \
                                (Bnorm * (nu_O2p_sum[lat, lon, lev] ** 2 + omega_O2p ** 2))
                Vi_O2p_star_z = (nu_O2p_sum[lat, lon, lev] * omega_O2p * Estar[2] + omega_O2p ** 2 * EstarXbunit[2]) / \
                                (Bnorm * (nu_O2p_sum[lat, lon, lev] ** 2 + omega_O2p ** 2))
                # ############ NO+ ############
                Vi_NOp_star_x = (nu_NOp_sum[lat, lon, lev] * omega_NOp * Estar[0] + omega_NOp ** 2 * EstarXbunit[0]) / \
                                (Bnorm * (nu_NOp_sum[lat, lon, lev] ** 2 + omega_NOp ** 2))
                Vi_NOp_star_y = (nu_NOp_sum[lat, lon, lev] * omega_NOp * Estar[1] + omega_NOp ** 2 * EstarXbunit[1]) / \
                                (Bnorm * (nu_NOp_sum[lat, lon, lev] ** 2 + omega_NOp ** 2))
                Vi_NOp_star_z = (nu_NOp_sum[lat, lon, lev] * omega_NOp * Estar[2] + omega_NOp ** 2 * EstarXbunit[2]) / \
                                (Bnorm * (nu_NOp_sum[lat, lon, lev] ** 2 + omega_NOp ** 2))

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

                # Measurements for each ion cannot be made separately so we take the average
                # no coordinate system assigned to satellite, assumed as point measurement
                # perpendicular ion velocity corresponds to measured ion velocity
                # in contrast to neutral wind, which we get as model input

                # Ion velocity perpendicular to magnetic field
                Vi_vertx[lat, lon, lev] = (Vi_Op_x + Vi_O2p_x + Vi_NOp_x) / 3
                Vi_verty[lat, lon, lev] = (Vi_Op_y + Vi_O2p_y + Vi_NOp_y) / 3
                Vi_vertz[lat, lon, lev] = (Vi_Op_z + Vi_O2p_z + Vi_NOp_z) / 3

                # Vi perpendicular to magnetic field vector
                Vi_vert = [Vi_vertx[lat, lon, lev], Vi_verty[lat, lon, lev], Vi_vertz[lat, lon, lev]]

                # ################################## JOULE HEATING ##################################
                # ###################################################################################
                # Joule Heating = qeNe(Vi_vert - Un_vert)dot(E_vert + Un_vert cross B)(in watt/m^3)
                # qe(in coulomb), B(in tesla), E(in volt/meter), Vi,Un(in meter/sec)
                Joule_Heating[lat, lon, lev] = qe * (Ne[lat, lon, lev] * ccm) * \
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
                Ohmic_Heating[lat, lon, lev] = pedersen_con[lat, lon, lev] * (term_ohm_x + term_ohm_y + term_ohm_z)
                # ############################### FRICTIONAL HEATING ################################
                # ###################################################################################
                # Frictional Heating = m_ion * nu_ion * N_ion * |Vi_vert - Un_vert|^2 (in watt/m^3)
                # m_ion(in kg), nu_ion(in Hz), N_ion(in m^(-3)), Vi,Un(in meter/sec)
                term_fric_x = (Vi_vert[0] - Un_vert[0]) ** 2
                term_fric_y = (Vi_vert[1] - Un_vert[1]) ** 2
                term_fric_z = (Vi_vert[2] - Un_vert[2]) ** 2
                term_Op = mkO * nu_Op_sum[lat, lon, lev] * (NOp[lat, lon, lev] * ccm)
                term_O2p = mkO2 * nu_O2p_sum[lat, lon, lev] * (NO2p[lat, lon, lev] * ccm)
                term_NOp = mkNO * nu_NOp_sum[lat, lon, lev] * (NNOp[lat, lon, lev] * ccm)
                Frictional_Heating[lat, lon, lev] = (term_Op + term_O2p + term_NOp) * (term_fric_x + term_fric_y + term_fric_z)
                # ############################ CROSS SECTIONS ############################
                # C = (nu_ion / N_neutral) / (sqrt(2 * boltzmann * T_i / m_ion)) (in m^2)
                # nu(in Hz), N(in m^(-3)), T(in kelvin), mass(in kg)
                N_neutral = NO[lat, lon, lev] + NO2[lat, lon, lev] + NN2[lat, lon, lev]
                N_neutral = N_neutral * ccm
                # ####### O+ #######
                C_Op[lat, lon, lev] = (nu_Op_sum[lat, lon, lev] / N_neutral) / (np.sqrt(2 * boltzmann * Ti[lat, lon, lev] / mkO))
                # ####### O2+ #######
                C_O2p[lat, lon, lev] = (nu_O2p_sum[lat, lon, lev] / N_neutral) / (np.sqrt(2 * boltzmann * Ti[lat, lon, lev] / mkO2))
                # ####### NO+ #######
                C_NOp[lat, lon, lev] = (nu_NOp_sum[lat, lon, lev] / N_neutral) / (np.sqrt(2 * boltzmann * Ti[lat, lon, lev] / mkNO))
                # ####### ION #######
                nu_ion = nu_Op_sum[lat, lon, lev] + nu_O2p_sum[lat, lon, lev] + nu_NOp_sum[lat, lon, lev]
                # Average collision frequency
                nu_ion = nu_ion / 3  # in Hz
                m_ion = mkO + mkO2 + mkNO
                # Average mass
                m_ion = m_ion / 3  # in kg
                # Because measurements for each species cannot be made we assume an average ion cross section
                C_ion[lat, lon, lev] = (nu_ion / N_neutral) / (np.sqrt(2 * boltzmann * Ti[lat, lon, lev] / m_ion))
                # ################################# PERPENDICULAR CURRENTS ##############################
                # ###############################  1st Methodology - Ohms law ###########################
                # Pedersen current = sigmaPedersen * E_star = (Evert + Unvert cross B) (in ampere/meter^2)
                # sigmaPedersen(in siemens/meter), E(in volt/meter)
                J_px = pedersen_con[lat, lon, lev] * Estar[0]
                J_py = pedersen_con[lat, lon, lev] * Estar[1]
                J_pz = pedersen_con[lat, lon, lev] * Estar[2]

                J_pedersen[lat, lon, lev] = np.sqrt(J_px ** 2 + J_py ** 2 + J_pz ** 2)

                # b_unit cross E_star
                x1 = bunit[1] * Estar[2] - bunit[2] * Estar[1]
                y1 = bunit[2] * Estar[0] - bunit[0] * Estar[2]
                z1 = bunit[0] * Estar[1] - bunit[1] * Estar[0]

                # Hall current = sigmaHall * (b_unit cross E_star) (in ampere/meter^2)
                # sigmaHall(in siemens/meter), E(in volt/meter)
                J_hx = hall_con[lat, lon, lev] * x1
                J_hy = hall_con[lat, lon, lev] * y1
                J_hz = hall_con[lat, lon, lev] * z1

                J_hall[lat, lon, lev] = np.sqrt(J_hx ** 2 + J_hy ** 2 + J_hz ** 2)

                # Ohmic current = |JPedersen + JHall|
                J_ohmic[lat, lon, lev] = np.sqrt((J_px + J_hx) ** 2 + (J_py + J_hy) ** 2 + (J_pz + J_hz) ** 2)

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
                J_denx = qe * (Ne[lat, lon, lev] * ccm) * (Vi_vert[0] - Ve_vert[0])
                J_deny = qe * (Ne[lat, lon, lev] * ccm) * (Vi_vert[1] - Ve_vert[1])
                J_denz = qe * (Ne[lat, lon, lev] * ccm) * (Vi_vert[2] - Ve_vert[2])
                J_dens[lat, lon, lev] = np.sqrt(J_denx ** 2 + J_deny ** 2 + J_denz ** 2)

                progress_bar.UpdateBar(i)

    time.sleep(3)
    window.close()
    Sg.popup("_" * 50, "Products calculated in : " + str(time.time() - start_time) + " sec!", "_" * 50, title="Finished", keep_on_top=True)
    print('Products calculated in: ', time.time() - start_time, ' sec!')
    print(' ')
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END OF PRODUCT CALCULATION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ CALCULATE ERROR $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def error(error_flag, B_error, E_error, NO_error, NO2_error, NN2_error, NOp_error, NO2p_error, NNOp_error, Ne_error, Te_error, Ti_error, Tn_error,
          Un_error, Vi_error, lat_value=-1, lon_value=-1, pressure_level=-1):

    start_time = time.time()
    print('Calculating Error.....')
    print(' ')

    # Distinguish Map from Vertical profile
    lev_range = 0
    lat_range = 0
    lon_range = 0

    lev_start = 0
    lat_start = 0
    lon_start = 0
    max_bar = 0

    # Lat - Alt map profile
    if lat_value == -1 and pressure_level == -1:
        lev_range = len(glev_in) - 1
        lat_range = len(glat_in)
        lon_start = lon_value
        lon_range = lon_start + 1
        max_bar = lev_range

    # Lat - Lon map profile
    if lat_value == -1 and lon_value == -1:
        lev_start = pressure_level
        lev_range = lev_start + 1
        lat_range = len(glat_in)
        lon_range = len(glon_in)
        max_bar = lat_range

    # Vertical profile
    if lat_value != -1 and lon_value != -1:
        lat_start = lat_value
        lon_start = lon_value
        lat_range = lat_start + 1
        lon_range = lon_start + 1
        lev_range = len(glev_in) - 1
        max_bar = lev_range

    # Create loading bar
    layout = [[Sg.Text('Calculating Error.....')], [Sg.ProgressBar(max_bar, orientation='h', size=(50, 20), key='progressbar')]]

    window = Sg.Window(title=title, keep_on_top=True).Layout(layout)
    progress_bar = window.FindElement('progressbar')
    event, values = window.Read(timeout=1)

    for lev in range(lev_start, lev_range):
        for lat in range(lat_start, lat_range):
            for lon in range(lon_start, lon_range):
                # progress bar index
                if max_bar == 56:
                    i = lev
                else:
                    i = lat

                # Magnetic field vector(in tesla)
                B = [Bx[lat, lon, lev], By[lat, lon, lev], Bz[lat, lon, lev]]
                # Magnetic field norm and unit vector
                Bnorm = np.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)
                bunit = [B[0] / Bnorm, B[1] / Bnorm, B[2] / Bnorm]
                # Electric field vector(in volt/meter)
                E = [Ex[lat, lon, lev], Ey[lat, lon, lev], Ez[lat, lon, lev]]
                # Neutral wind vector(in meter/sec)
                Un = [Unx[lat, lon, lev], Uny[lat, lon, lev], Unz[lat, lon, lev]]
                # Ion velocity vector(in meter/sec)
                # Vi_vert ≡ Vi as mentioned in products function above
                Vi_vert = [Vi_vertx[lat, lon, lev], Vi_verty[lat, lon, lev], Vi_vertz[lat, lon, lev]]

                # ############################### ASSIGNING ERRORS ###############################
                # distinguish percentage from real errors
                if not error_flag:
                    # percentage errors, squared
                    dBx = ((B_error / 100) * B[0]) ** 2
                    dBy = ((B_error / 100) * B[1]) ** 2
                    dBz = ((B_error / 100) * B[2]) ** 2
                    # |theta_|B| / theta_Bx|^2
                    thB_Bx = (B[0] / Bnorm) ** 2
                    # |theta_|B| / theta_By|^2
                    thB_By = (B[1] / Bnorm) ** 2
                    # |theta_|B| / theta_Bz|^2
                    thB_Bz = (B[2] / Bnorm) ** 2
                    # d|B|^2
                    dB = thB_Bx * dBx + thB_By * dBy + thB_Bz * dBz
                    dEx = ((E_error / 100) * E[0]) ** 2
                    dEy = ((E_error / 100) * E[1]) ** 2
                    dEz = ((E_error / 100) * E[2]) ** 2
                    # ################### N in cm^(-3) ####################
                    dNO = ((NO_error / 100) * NO[lat, lon, lev]) ** 2
                    dNO2 = ((NO2_error / 100) * NO2[lat, lon, lev]) ** 2
                    dNN2 = ((NN2_error / 100) * NN2[lat, lon, lev]) ** 2
                    dNOp = ((NOp_error / 100) * NOp[lat, lon, lev]) ** 2
                    dNO2p = ((NO2p_error / 100) * NO2p[lat, lon, lev]) ** 2
                    dNNOp = ((NNOp_error / 100) * NNOp[lat, lon, lev]) ** 2
                    dNe = ((Ne_error / 100) * Ne[lat, lon, lev]) ** 2
                    # #####################################################
                    dTe = ((Te_error / 100) * Te[lat, lon, lev]) ** 2
                    dTi = ((Ti_error / 100) * Ti[lat, lon, lev]) ** 2
                    dTn = ((Tn_error / 100) * Tn[lat, lon, lev]) ** 2
                    dUnx = ((Un_error / 100) * Un[0]) ** 2
                    dUny = ((Un_error / 100) * Un[1]) ** 2
                    dUnz = ((Un_error / 100) * Un[2]) ** 2
                    dVix = ((Vi_error / 100) * Vi_vert[0]) ** 2
                    dViy = ((Vi_error / 100) * Vi_vert[1]) ** 2
                    dViz = ((Vi_error / 100) * Vi_vert[2]) ** 2
                elif error_flag:
                    # real errors, squared
                    dBx = () ** 2
                    dBy = () ** 2
                    dBz = () ** 2
                    # |theta_|B| / theta_Bx|^2
                    thB_Bx = (B[0] / Bnorm) ** 2
                    # |theta_|B| / theta_By|^2
                    thB_By = (B[1] / Bnorm) ** 2
                    # |theta_|B| / theta_Bz|^2
                    thB_Bz = (B[2] / Bnorm) ** 2
                    # d|B|^2
                    dB = thB_Bx * dBx + thB_By * dBy + thB_Bz * dBz

                    dEx = () ** 2
                    dEy = () ** 2
                    dEz = () ** 2
                    # ################### N in cm^(-3) ####################
                    dNO = () ** 2
                    dNO2 = () ** 2
                    dNN2 = () ** 2
                    dNOp = () ** 2
                    dNO2p = () ** 2
                    dNNOp = () ** 2
                    dNe = () ** 2
                    # #####################################################
                    dTe = () ** 2
                    dTi = () ** 2
                    dTn = () ** 2
                    dUnx = () ** 2
                    dUny = () ** 2
                    dUnz = () ** 2
                    dVix = () ** 2
                    dViy = () ** 2
                    dViz = () ** 2

                # ################# COLLISION FREQUENCIES ERROR #################
                # ######################### O+ #########################
                # |dnuOp|^2 = |dnuOp-O|^2 + |dnuOp-O2|^2 + |dnuOp-N2|^2
                Tr = (Ti[lat, lon, lev] + Tn[lat, lon, lev]) / 2
                # |theta_nuOp-O / theta_NO|^2
                thOp_O_NO = (fb * 3.67 * 10 ** (-11) * Tr ** (1 / 2) * (1 - 0.064 * np.log10(Tr)) ** 2) ** 2
                # |theta_nuOp-O / theta_Ti|^2
                thOp_O_Ti = ((9.175 * 10 ** (-12) * fb * NO[lat, lon, lev] * (1 - 0.064 * (np.log(Tr) / np.log(10))) ** 2) / Tr ** (1 / 2) -
                             (2.3488 * 10 ** (-12) * fb * NO[lat, lon, lev] * (1 - 0.064 * (np.log(Tr) / np.log(10)))) /
                             (np.log(10) * Tr ** (1 / 2))) ** 2
                # |theta_nuOp-O / theta_Tn|^2
                thOp_O_Tn = thOp_O_Ti
                # |dnuOp_O|^2
                dnuOp_O = thOp_O_NO * dNO + thOp_O_Ti * dTi + thOp_O_Tn * dTn

                # |theta_nuOp-O2 / theta_NO2|^2
                thOp_O2_NO2 = (6.64 * 10 ** (-10)) ** 2
                # |dnuOp_O2|^2
                dnuOp_O2 = thOp_O2_NO2 * dNO2

                # |theta_nuOp-N2 / theta_NN2|^2
                thOp_N2_NN2 = (6.82 * 10 ** (-10)) ** 2
                # |dnuOp_N2|^2
                dnuOp_N2 = thOp_N2_NN2 * dNN2

                # |dnuOp|
                nuOp_error[lat, lon, lev] = np.sqrt(dnuOp_O + dnuOp_O2 + dnuOp_N2)
                # ############################# O2+ #############################
                # |dnuO2p|^2 = |dnuO2p-O2|^2 + |dnuOp-O|^2 + |dnuOp-N2|^2
                # |theta_nuO2p-O2 / theta_NO2|^2
                thO2p_O2_NO2 = (2.59 * 10 ** (-11) * Tr ** (1 / 2) * (1 - 0.073 * np.log10(Tr)) ** 2) ** 2
                # |theta_nuO2p-O2 / theta_Ti|^2
                thO2p_O2_Ti = ((6.475 * 10 ** (-12) * NO2[lat, lon, lev] * (1 - 0.073 * (np.log(Tr) / np.log(10))) ** 2) / Tr ** (1 / 2) -
                               (1.8907 * 10 ** (-12) * NO2[lat, lon, lev] * (1 - 0.073 * (np.log(Tr) / np.log(10)))) /
                               (np.log(10) * Tr ** (1 / 2))) ** 2
                # |theta_nuO2p-O2 / theta_Tn|^2
                thO2p_O2_Tn = thO2p_O2_Ti

                # |dnuO2p_O2|^2
                dnuO2p_O2 = thO2p_O2_NO2 * dNO2 + thO2p_O2_Ti * dTi + thO2p_O2_Tn * dTn

                # |theta_nuO2p-O / theta_NO|^2
                thO2p_O_NO = (2.31 * 10 ** (-10)) ** 2
                # |dnuO2p_O|^2
                dnuO2p_O = thO2p_O_NO * dNO
                # |theta_nuO2p-N2 / theta_NN2|^2
                thO2p_N2_NN2 = (4.13 * 10 ** (-10)) ** 2
                # |dnuO2p_N2|^2
                dnuO2p_N2 = thO2p_N2_NN2 * dNN2

                # |dnuO2p|
                nuO2p_error[lat, lon, lev] = np.sqrt(dnuO2p_O2 + dnuO2p_O + dnuO2p_N2)
                # ########################## NO+ ###########################
                # |dnuNOp|^2 = |dnuNOp-O|^2 + |dnuNOp-O2|^2 + |dnuNOp-N2|^2
                # |theta_nuNOp-O / theta_NO|^2
                thNOp_O_NO = (2.44 * 10 ** (-10)) ** 2
                # |dnuNOp_O|^2
                dnuNOp_O = thNOp_O_NO * dNO
                # |theta_nuNOp_O2 / theta_NO2|^2
                thNOp_O2_NO2 = (4.27 * 10 ** (-10)) ** 2
                # |dnuNOp_O2|^2
                dnuNOp_O2 = thNOp_O2_NO2 * dNO2
                # |theta_nuNOp-N2 / theta_NN2|^2
                thNOp_N2_NN2 = (4.34 * 10 ** (-10)) ** 2
                # |dnuNOp_N2|^2
                dnuNOp_N2 = thNOp_N2_NN2 * dNN2

                # |dnuNOp|
                nuNOp_error[lat, lon, lev] = np.sqrt(dnuNOp_O + dnuNOp_N2 + dnuNOp_O2)
                # ############################ Ion collision frequency contributions error ############################
                # squared
                dnuion_Ti[lat, lon, lev] = thOp_O_Ti * dTi + thO2p_O2_Ti * dTi
                dnuion_Tn[lat, lon, lev] = thOp_O_Tn * dTn + thO2p_O2_Tn * dTn
                dnuion_Nneutral[lat, lon, lev] = ((thOp_O_NO + thO2p_O_NO + thNOp_O_NO) / 9) * dNO + (
                                                  (thOp_O2_NO2 + thO2p_O2_NO2 + thNOp_O2_NO2) / 9) * dNO2 + (
                                                  (thOp_N2_NN2 + thO2p_N2_NN2 + thNOp_N2_NN2) / 9) * dNN2
                # ######################################################################################################
                # ######################## e #######################
                # |dnue|^2 = |dnue-O|^2 + |dnue-O2|^2 + |dnue-N2|^2
                # |theta_nue-O / theta_Te|^2
                the_O_Te = (4.45 * 10 ** (-11) * NO[lat, lon, lev] * Te[lat, lon, lev] ** (-1 / 2) +
                            7.6095 * 10 ** (-14) * NO[lat, lon, lev] * Te[lat, lon, lev] ** (1 / 2)) ** 2
                # |theta_nue-O / theta_NO|^2
                the_O_NO = (8.9 * 10 ** (-11) * Te[lat, lon, lev] ** (1 / 2) * (1 + 5.7 * 10 ** (-4) * Te[lat, lon, lev])) ** 2
                # |dne_O|^2
                dnue_O = the_O_Te * dTe + the_O_NO * dNO
                # |theta_nue-O2 / theta_Te|^2
                the_O2_Te = (9.1 * 10 ** (-11) * NO2[lat, lon, lev] * Te[lat, lon, lev] ** (-1 / 2) +
                             6.552 * 10 ** (-12) * NO2[lat, lon, lev]) ** 2
                # |theta_nue-O2 / theta_NO2|^2
                the_O2_NO2 = (1.82 * 10 ** (-10) * Te[lat, lon, lev] ** (1 / 2) * (1 + 3.6 * 10 ** (-2) * Te[lat, lon, lev] ** (1 / 2))) ** 2
                # |dnue_O2|^2
                dnue_O2 = the_O2_Te * dTe + the_O2_NO2 * dNO2
                # |theta_nue-N2 / theta_Te|^2
                the_N2_Te = (2.33 * 10 ** (-11) * NN2[lat, lon, lev] - 5.6386 * 10 ** (-15) * NN2[lat, lon, lev] * Te[lat, lon, lev]) ** 2
                # |theta_nue-N2 / theta_NN2|^2
                the_N2_NN2 = (2.33 * 10 ** (-11) * Te[lat, lon, lev] * (1 - 1.21 * 10 ** (-4) * Te[lat, lon, lev])) ** 2
                # |dnue_N2|^2
                dnue_N2 = the_N2_Te * dTe + the_N2_NN2 * dNN2

                # |dnue|
                nue_error[lat, lon, lev] = np.sqrt(dnue_O + dnue_O2 + dnue_N2)
                # ############################### e collision frequency contributions error ##############################
                # squared
                dnue_Te[lat, lon, lev] = the_O_Te * dTe + the_O2_Te * dTe + the_N2_Te * dTe
                dnue_Nneutral[lat, lon, lev] = the_O_NO * dNO + the_O2_NO2 * dNO2 + the_N2_NN2 * dNN2
                # ########################################################################################################
                # ###################### OMEGAS ERROR ######################

                # qe(in coulomb), mk(masses in kg), omegas(in Hz)
                omega_Op = (qe * Bnorm) / mkO
                omega_O2p = (qe * Bnorm) / mkO2
                omega_NOp = (qe * Bnorm) / mkNO
                omega_e = (qe * Bnorm) / me

                # ############## O+ ##############
                # |theta_omegaOp / theta_B|^2
                th_omOp_B = (qe / mkO) ** 2
                # |domegaOp|^2
                domegaOp = th_omOp_B * dB
                # ############# O2+ ##############
                # |theta_omegaO2p / theta_B|^2
                th_omO2p_B = (qe / mkO2) ** 2
                # |domegaO2p|^2
                domegaO2p = th_omO2p_B * dB
                # ############# NO+ ##############
                # |theta_omegaNOp / theta_B|^2
                th_omNOp_B = (qe / mkNO) ** 2
                # |domegaOp|^2
                domegaNOp = th_omNOp_B * dB
                # ############## e ###############
                # |theta_omega_e / theta_B|^2
                th_ome_B = (qe / me) ** 2
                # |domega_e|^2
                domegae = th_ome_B * dB
                # ####################### RATIOS ERROR #######################
                r_Op = nu_Op_sum[lat, lon, lev] / omega_Op
                r_O2p = nu_O2p_sum[lat, lon, lev] / omega_O2p
                r_NOp = nu_NOp_sum[lat, lon, lev] / omega_NOp
                r_e = nu_e_sum[lat, lon, lev] / omega_e
                # ############ O+ ###########
                # |theta_rOp / theta_nuOp|^2
                th_rOp_nuOp = (1 / omega_Op) ** 2
                # |theta_rOp / theta_omegaOp|^2
                th_rOp_omOp = (nu_Op_sum[lat, lon, lev] / omega_Op ** 2) ** 2
                # |drOp|^2
                drOp = th_rOp_nuOp * nuOp_error[lat, lon, lev] ** 2 + th_rOp_omOp * domegaOp
                # ########### O2+ ###########
                # |theta_rO2p / theta_nuO2p|^2
                th_rO2p_nuO2p = (1 / omega_O2p) ** 2
                # |theta_rO2p / theta_omegaO2p|^2
                th_rO2p_omO2p = (nu_O2p_sum[lat, lon, lev] / omega_O2p ** 2) ** 2
                # |drO2p|^2
                drO2p = th_rO2p_nuO2p * nuO2p_error[lat, lon, lev] ** 2 + th_rO2p_omO2p * domegaO2p
                # ########### NO+ ###########
                # |theta_rNOp / theta_nuNOp|^2
                th_rNOp_nuNOp = (1 / omega_NOp) ** 2
                # |theta_rNOp / theta_omegaNOp|^2
                th_rNOp_omNOp = (nu_NOp_sum[lat, lon, lev] / omega_NOp ** 2) ** 2
                # |drNOp|^2
                drNOp = th_rNOp_nuNOp * nuNOp_error[lat, lon, lev] ** 2 + th_rNOp_omNOp * domegaNOp
                # ########## e ##########
                # |theta_re / theta_nue|^2
                th_re_nue = (1 / omega_e) ** 2
                # |theta_re / theta_omegae|^2
                th_re_ome = (nu_e_sum[lat, lon, lev] / omega_e ** 2) ** 2
                # |dre|^2
                dre = th_re_nue * nue_error[lat, lon, lev] ** 2 + th_re_ome * domegae
                # ################################ CONDUCTIVITIES ERROR ################################
                # ############################# PEDERSEN CONDUCTIVITY ERROR ############################
                # densities error in m^(-3)
                # |theta_sigmaP/ theta_B|^2
                term_e = (Ne[lat, lon, lev] * ccm * r_e) / (1 + r_e ** 2)
                term_Op = (NOp[lat, lon, lev] * ccm * r_Op) / (1 + r_Op ** 2)
                term_O2p = (NO2p[lat, lon, lev] * ccm * r_O2p) / (1 + r_O2p ** 2)
                term_NOp = (NNOp[lat, lon, lev] * ccm * r_NOp) / (1 + r_NOp ** 2)
                thP_B = (qe * (term_e + term_Op + term_O2p + term_NOp) / Bnorm ** 2) ** 2
                # |theta_sigmaP / theta_Ne|^2
                thP_Ne = (qe * r_e / (Bnorm * (1 + r_e ** 2))) ** 2
                # |theta_sigmaP / theta_re|^2
                thP_re = (qe * Ne[lat, lon, lev] * ccm * (1 - r_e ** 2) / (Bnorm * (1 + r_e ** 2) ** 2)) ** 2
                # |theta_sigmaP / theta_NOp|^2
                thP_NOp = (qe * r_Op / (Bnorm * (1 + r_Op ** 2))) ** 2
                # |theta_sigmaP / theta_rOp|^2
                thP_rOp = (qe * NOp[lat, lon, lev] * ccm * (1 - r_Op ** 2) / (Bnorm * (1 + r_Op ** 2) ** 2)) ** 2
                # |theta_sigmaP / theta_NO2p|^2
                thP_NO2p = (qe * r_O2p / (Bnorm * (1 + r_O2p ** 2))) ** 2
                # |theta_sigmaP / theta_rO2p|^2
                thP_rO2p = (qe * NO2p[lat, lon, lev] * ccm * (1 - r_O2p ** 2) / (Bnorm * (1 + r_O2p ** 2) ** 2)) ** 2
                # |theta_sigmaP / theta_NNOp|^2
                thP_NNOp = (qe * r_NOp / (Bnorm * (1 + r_NOp ** 2))) ** 2
                # |theta_sigmaP / theta_rOp|^2
                thP_rNOp = (qe * NNOp[lat, lon, lev] * ccm * (1 - r_NOp ** 2) / (Bnorm * (1 + r_NOp ** 2) ** 2)) ** 2
                # |dsigmaP|
                pedersen_con_error[lat, lon, lev] = np.sqrt(thP_B * dB + thP_Ne * dNe * ccm ** 2 + thP_re * dre +
                                                            thP_NOp * dNOp * ccm ** 2 + thP_rOp * drOp +
                                                            thP_NO2p * dNO2p * ccm ** 2 + thP_rO2p * drO2p +
                                                            thP_NNOp * dNNOp * ccm ** 2 + thP_rNOp * drNOp)
                # ############################ Pedersen conductivity contributions error ############################
                # squared
                dsp_B[lat, lon, lev] = thP_B * dB + thP_re * th_re_ome * th_ome_B * dB + thP_rOp * th_rOp_omOp * th_omOp_B * dB + \
                                       thP_rO2p * th_rO2p_omO2p * th_omO2p_B * dB + thP_rNOp * th_rNOp_omNOp * th_omNOp_B * dB
                dsp_Ti[lat, lon, lev] = thP_rOp * th_rOp_nuOp * thOp_O_Ti * dTi + thP_rO2p * th_rO2p_nuO2p * thO2p_O2_Ti * dTi
                dsp_Te[lat, lon, lev] = thP_re * th_re_nue * the_O_Te * dTe + thP_re * th_re_nue * the_O2_Te * dTe + \
                                        thP_re * th_re_nue * the_N2_Te * dTe
                dsp_Tn[lat, lon, lev] = thP_rOp * th_rOp_nuOp * thOp_O_Tn * dTn + thP_rO2p * th_rO2p_nuO2p * thO2p_O2_Tn * dTn
                dsp_Nion[lat, lon, lev] = thP_NOp * dNOp * ccm ** 2 + thP_NO2p * dNO2p * ccm ** 2 + thP_NNOp * dNNOp * ccm ** 2
                dsp_Ne[lat, lon, lev] = thP_Ne * dNe * ccm ** 2
                dsp_Nneutral[lat, lon, lev] = thP_rOp * th_rOp_nuOp * thOp_O_NO * dNO + thP_rOp * th_rOp_nuOp * thOp_O2_NO2 * dNO2 + \
                                              thP_rOp * th_rOp_nuOp * thOp_N2_NN2 * dNN2 + thP_rO2p * th_rO2p_nuO2p * thO2p_O2_NO2 * dNO2 + \
                                              thP_rO2p * th_rO2p_nuO2p * thO2p_O_NO * dNO + thP_rO2p * th_rO2p_nuO2p * thO2p_N2_NN2 * dNN2 + \
                                              thP_rNOp * th_rNOp_nuNOp * thNOp_O_NO * dNO + thP_rNOp * th_rNOp_nuNOp * thNOp_O2_NO2 * dNO2 + \
                                              thP_rNOp * th_rNOp_nuNOp * thNOp_N2_NN2 * dNN2 + thP_re * th_re_nue * the_O_NO * dNO + \
                                              thP_re * th_re_nue * the_O2_NO2 * dNO2 + thP_re * th_re_nue * the_N2_NN2 * dNN2
                # ####################################################################################################
                # ####################### HALL CONDUCTIVITY ERROR #######################
                # densities error in m^(-3)
                # |theta_sigmaH / theta_B|^2
                thH_B = (qe * (term_e / r_e - term_Op / r_Op - term_O2p / r_O2p - term_NOp / r_NOp) / Bnorm ** 2) ** 2
                # |theta_sigmaH / theta_Ne|^2
                thH_Ne = (qe / (Bnorm * (1 + r_e ** 2))) ** 2
                # |theta_sigmaH / theta_re|^2
                thH_re = (qe * Ne[lat, lon, lev] * ccm * 2 * r_e / (Bnorm * (1 + r_e ** 2) ** 2)) ** 2
                # |theta_sigmaH / theta_NOp|^2
                thH_NOp = (qe / (Bnorm * (1 + r_Op ** 2))) ** 2
                # |theta_sigmaH / theta_rOp|^2
                thH_rOp = (qe * NOp[lat, lon, lev] * ccm * 2 * r_Op / (Bnorm * (1 + r_Op ** 2) ** 2)) ** 2
                # |theta_sigmaH / theta_NO2p|^2
                thH_NO2p = (qe / (Bnorm * (1 + r_O2p ** 2))) ** 2
                # |theta_sigmaH / theta_rO2p|^2
                thH_rO2p = (qe * NO2p[lat, lon, lev] * ccm * 2 * r_O2p / (Bnorm * (1 + r_O2p ** 2) ** 2)) ** 2
                # |theta_sigmaH / theta_NNOp|^2
                thH_NNOp = (qe / (Bnorm * (1 + r_NOp ** 2))) ** 2
                # |theta_sigmaH / theta_rNOp|^2
                thH_rNOp = (qe * NNOp[lat, lon, lev] * ccm * 2 * r_NOp / (Bnorm * (1 + r_NOp ** 2) ** 2)) ** 2
                # |dsigmaH|
                hall_con_error[lat, lon, lev] = np.sqrt(thH_B * dB + thH_Ne * dNe * ccm ** 2 + thH_re * dre +
                                                        thH_NOp * dNOp * ccm ** 2 + thH_rOp * drOp +
                                                        thH_NO2p * dNO2p * ccm ** 2 + thH_rO2p * drO2p +
                                                        thH_NNOp * dNNOp * ccm ** 2 + thH_rNOp * drNOp)
                # ###################################### Hall conductivity contributions error ######################################
                dsh_B[lat, lon, lev] = thH_B * dB + thH_re * th_re_ome * th_ome_B * dB + thH_rOp * th_rOp_omOp * th_omOp_B * dB + \
                                       thH_rO2p * th_rO2p_omO2p * th_omO2p_B * dB + thH_rNOp * th_rNOp_omNOp * th_omNOp_B * dB
                dsh_Ti[lat, lon, lev] = thH_rOp * th_rOp_nuOp * thOp_O_Ti * dTi + thH_rO2p * th_rO2p_nuO2p * thO2p_O2_Ti * dTi
                dsh_Te[lat, lon, lev] = thH_re * th_re_nue * the_O_Te * dTe + thH_re * th_re_nue * the_O2_Te * dTe + \
                                        thH_re * th_re_nue * the_N2_Te * dTe
                dsh_Tn[lat, lon, lev] = thH_rOp * th_rOp_nuOp * thOp_O_Tn * dTn + thH_rO2p * th_rO2p_nuO2p * thO2p_O2_Tn * dTn
                dsh_Nion[lat, lon, lev] = thH_NOp * dNOp * ccm ** 2 + thH_NO2p * dNO2p * ccm ** 2 + thH_NNOp * dNNOp * ccm ** 2
                dsh_Ne[lat, lon, lev] = thH_Ne * dNe * ccm ** 2
                dsh_Nneutral[lat, lon, lev] = thH_rOp * th_rOp_nuOp * thOp_O_NO * dNO + thH_rOp * th_rOp_nuOp * thOp_O2_NO2 * dNO2 + \
                                              thH_rOp * th_rOp_nuOp * thOp_N2_NN2 * dNN2 + thH_rO2p * th_rO2p_nuO2p * thO2p_O2_NO2 * dNO2 + \
                                              thH_rO2p * th_rO2p_nuO2p * thO2p_O_NO * dNO + thH_rO2p * th_rO2p_nuO2p * thO2p_N2_NN2 * dNN2 + \
                                              thH_rNOp * th_rNOp_nuNOp * thNOp_O_NO * dNO + thH_rNOp * th_rNOp_nuNOp * thNOp_O2_NO2 * dNO2 + \
                                              thH_rNOp * th_rNOp_nuNOp * thNOp_N2_NN2 * dNN2 + thH_re * th_re_nue * the_O_NO * dNO + \
                                              thH_re * th_re_nue * the_O2_NO2 * dNO2 + thH_re * th_re_nue * the_N2_NN2 * dNN2
                # #####################################################################################################################
                # ############################# PARALLEL CONDUCTIVITY ERROR #############################
                # densities error in m^(-3)
                # |theta_sigma_paral / theta_Ne|^2
                thparal_Ne = (qe ** 2 / (me * nu_e_sum[lat, lon, lev])) ** 2
                # |theta_sigma_paral / theta_nue|^2
                thparal_nue = (Ne[lat, lon, lev] * ccm * qe ** 2 / (me * nu_e_sum[lat, lon, lev] ** 2)) ** 2
                # |dsigma_paral|
                parallel_con_error[lat, lon, lev] = np.sqrt(thparal_Ne * dNe * ccm ** 2 + thparal_nue * nue_error[lat, lon, lev] ** 2)
                # ########################### HEATING RATES ERROR ###########################
                # ################## bunit error ##################
                # ######### bx #########
                # |theta_bx / theta_Bx|^2
                thbx_Bx = (1 / Bnorm) ** 2
                # |theta_bx / theta_B|^2
                thbx_B = (B[0] / Bnorm ** 2) ** 2
                # |dbx|^2
                dbx = thbx_Bx * dBx + thbx_B * dB
                # ######### by #########
                # |theta_by / theta_By|^2
                thby_By = (1 / Bnorm) ** 2
                # |theta_by / theta_B|^2
                thby_B = (B[1] / Bnorm ** 2) ** 2
                # |dby|^2
                dby = thby_By * dBy + thby_B * dB
                # ######### bz #########
                # |theta_bz / theta_Bz|^2
                thbz_Bz = (1 / Bnorm) ** 2
                # |theta_bz / theta_B|^2
                thbz_B = (B[2] / Bnorm ** 2) ** 2
                # |dbz|^2
                dbz = thbz_Bz * dBz + thbz_B * dB
                # ############### Un vertical error ###############
                # ########## Un_vertx ##########
                # |theta_Unvertx / theta_Uny|^2
                thUn_vertx_Uny = (bunit[2]) ** 2
                # |theta_Unvertx / theta_Unz|^2
                thUn_vertx_Unz = (bunit[1]) ** 2
                # |theta_Unvertx / theta_bz|^2
                thUn_vertx_bz = (Un[1]) ** 2
                # |theta_Unvertx / theta_by|^2
                thUn_vertx_by = (Un[2]) ** 2
                # |dUn_vertx|^2
                dUn_vertx = thUn_vertx_Uny * dUny + thUn_vertx_Unz * dUnz + thUn_vertx_bz * dbz + thUn_vertx_by * dby
                # ########## Un_verty ##########
                # |theta_Unverty / theta_Unz|^2
                thUn_verty_Unz = (bunit[0]) ** 2
                # |theta_Unverty / theta_Unx|^2
                thUn_verty_Unx = (bunit[2]) ** 2
                # |theta_Unverty / theta_bx|^2
                thUn_verty_bx = (Un[2]) ** 2
                # |theta_Unverty / theta_bz|^2
                thUn_verty_bz = (Un[0]) ** 2
                # |dUn_verty|^2
                dUn_verty = thUn_verty_Unz * dUnz + thUn_verty_Unx * dUnx + thUn_verty_bx * dbx + thUn_verty_bz * dbz
                # ########## Un_vertz ##########
                # |theta_Unvertz / theta_Unx|^2
                thUn_vertz_Unx = (bunit[1]) ** 2
                # |theta_Unvertz / theta_Uny|^2
                thUn_vertz_Uny = (bunit[0]) ** 2
                # |theta_Unvertz / theta_by|^2
                thUn_vertz_by = (Un[0]) ** 2
                # |theta_Unvertz / theta_bx|^2
                thUn_vertz_bx = (Un[1]) ** 2
                # |dUn_vertz|^2
                dUn_vertz = thUn_vertz_Unx * dUnx + thUn_vertz_Uny * dUny + thUn_vertz_by * dby + thUn_vertz_bx * dbx
                # ###################### E vertical error ######################
                # ######### E_vertx #########
                # |theta_Evertx / theta_Ey|^2
                thEvertx_Ey = (bunit[2]) ** 2
                # |theta_Evertx / theta_Ez|^2
                thEvertx_Ez = (bunit[1]) ** 2
                # |theta_Evertx / theta_bz|^2
                thEvertx_bz = (E[1]) ** 2
                # |theta_Evertx / theta_by|^2
                thEvertx_by = (E[2]) ** 2
                # |dEvertx|^2
                dEvertx = thEvertx_Ey * dEy + thEvertx_Ez * dEz + thEvertx_bz * dbz + thEvertx_by * dby
                # ######### E_verty #########
                # |theta_Everty / theta_Ez|^2
                thEverty_Ez = (bunit[0]) ** 2
                # |theta_Everty / theta_Ex|^2
                thEverty_Ex = (bunit[2]) ** 2
                # |theta_Everty / theta_bx|^2
                thEverty_bx = (E[2]) ** 2
                # |theta_Everty / theta_bz|^2
                thEverty_bz = (E[0]) ** 2
                # |dEverty|^2
                dEverty = thEverty_Ez * dEz + thEverty_Ex * dEx + thEverty_bx * dbx + thEverty_bz * dbz
                # ######### E_vertz #########
                # |theta_Evertz / theta_Ex|^2
                thEvertz_Ex = (bunit[1]) ** 2
                # |theta_Evertz / theta_Ey|^2
                thEvertz_Ey = (bunit[0]) ** 2
                # |theta_Evertz / theta_by|^2
                thEvertz_by = (E[0]) ** 2
                # |theta_Evertz / theta_bx|^2
                thEvertz_bx = (E[1]) ** 2
                # |dEvertz|^2
                dEvertz = thEvertz_Ex * dEx + thEvertz_Ey * dEy + thEvertz_by * dby + thEvertz_bx * dbx
                # ############################ JOULE HEATING ERROR ############################
                # #############################################################################
                # densities error in m^(-3)
                # Electric field perpendicular to magnetic field
                # Evert = E cross bunit
                Evertx = E[1] * bunit[2] - E[2] * bunit[1]
                Everty = E[2] * bunit[0] - E[0] * bunit[2]
                Evertz = E[0] * bunit[1] - E[1] * bunit[0]

                # E vertical vector
                Evert = [Evertx, Everty, Evertz]

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

                # Estar: perpendicular electric field in the neutral frame, Estar = Evert + Unvert cross B
                # vector addition
                Estar_x = Evert[0] + UnvXB[0]
                Estar_y = Evert[1] + UnvXB[1]
                Estar_z = Evert[2] + UnvXB[2]

                Estar = [Estar_x, Estar_y, Estar_z]

                # |theta_JH / theta_Ne|^2
                thJH_Ne = (Joule_Heating[lat, lon, lev] / (Ne[lat, lon, lev] * ccm)) ** 2
                # |theta_JH / theta_Vi_vertx|^2
                thJH_Vi_vertx = (qe * Ne[lat, lon, lev] * ccm * Estar[0]) ** 2
                # |theta_JH / theta_Vi_verty|^2
                thJH_Vi_verty = (qe * Ne[lat, lon, lev] * ccm * Estar[1]) ** 2
                # |theta_JH / theta_Vi_vertz|^2
                thJH_Vi_vertz = (qe * Ne[lat, lon, lev] * ccm * Estar[2]) ** 2
                # |theta_JH / theta_Un_vertx|^2
                thJH_Un_vertx = (qe * Ne[lat, lon, lev] * ccm * (Vi_vert[2] * B[1] - Vi_vert[1] * B[2] - Evert[0])) ** 2
                # |theta_JH / theta_Un_verty|^2
                thJH_Un_verty = (qe * Ne[lat, lon, lev] * ccm * (Vi_vert[0] * B[2] - Vi_vert[2] * B[0] - Evert[1])) ** 2
                # |theta_JH / theta_Un_vertz|^2
                thJH_Un_vertz = (qe * Ne[lat, lon, lev] * ccm * (Vi_vert[1] * B[0] - Vi_vert[0] * B[1] - Evert[2])) ** 2
                # |theta_JH / theta_Evertx|^2
                thJH_Evertx = (qe * Ne[lat, lon, lev] * ccm * (Vi_vert[0] - Un_vert[0])) ** 2
                # |theta_JH / theta_Everty|^2
                thJH_Everty = (qe * Ne[lat, lon, lev] * ccm * (Vi_vert[1] - Un_vert[1])) ** 2
                # |theta_JH / theta_Evertz|^2
                thJH_Evertz = (qe * Ne[lat, lon, lev] * ccm * (Vi_vert[2] - Un_vert[2])) ** 2
                # |theta_JH / theta_Bx|^2
                thJH_Bx = (qe * Ne[lat, lon, lev] * ccm * (Vi_vert[1] * Un_vert[2] - Vi_vert[2] * Un_vert[1])) ** 2
                # |theta_JH / theta_By|^2
                thJH_By = (qe * Ne[lat, lon, lev] * ccm * (Vi_vert[2] * Un_vert[0] - Vi_vert[0] * Un_vert[2])) ** 2
                # |theta_JH / theta_Bz|^2
                thJH_Bz = (qe * Ne[lat, lon, lev] * ccm * (Vi_vert[0] * Un_vert[1] - Vi_vert[1] * Un_vert[0])) ** 2

                # |dJH|
                Joule_Heating_error[lat, lon, lev] = np.sqrt(thJH_Ne * dNe * ccm ** 2 + thJH_Vi_vertx * dVix + thJH_Vi_verty * dViy +
                                                             thJH_Vi_vertz * dViz + thJH_Un_vertx * dUn_vertx + thJH_Un_verty * dUn_verty +
                                                             thJH_Un_vertz * dUn_vertz + thJH_Evertx * dEvertx + thJH_Everty * dEverty +
                                                             thJH_Evertz * dEvertz + thJH_Bx * dBx + thJH_By * dBy + thJH_Bz * dBz)
                # ######################################### Joule Heating contributions error #########################################
                # squared
                dJH_B[lat, lon, lev] = thJH_Bx * dBx + thJH_Un_verty * thUn_verty_bx * dbx + thJH_Un_vertz * thUn_vertz_bx * dbx + \
                                       thJH_Everty * thEverty_bx * dbx + thJH_Evertz * thEvertz_bx * dbx + thJH_By * dBy + \
                                       thJH_Un_vertx * thUn_vertx_by * dby + thJH_Un_vertz * thUn_vertz_by * dby + thJH_Evertx * thEvertx_by * dby + \
                                       thJH_Evertz * thEvertz_by * dby + thJH_Bz * dBz + thJH_Un_vertx * thUn_vertx_bz * dbz + \
                                       thJH_Un_verty * thUn_verty_bz * dbz + thJH_Evertx * thEvertx_bz * dbz + thJH_Everty * thEverty_bz * dbz
                dJH_E[lat, lon, lev] = thJH_Evertx * thEvertx_Ey * dEy + thJH_Evertx * thEvertx_Ez * dEz + thJH_Everty * thEverty_Ex * dEx + \
                                       thJH_Everty * thEverty_Ez * dEz + thJH_Evertz * thEvertz_Ex * dEx + thJH_Evertz * thEvertz_Ey * dEy
                dJH_Vi[lat, lon, lev] = thJH_Vi_vertx * dVix + thJH_Vi_verty * dViy + thJH_Vi_vertz * dViz
                dJH_Un[lat, lon, lev] = thJH_Un_vertx * thUn_vertx_Uny * dUny + thJH_Un_vertx * thUn_vertx_Unz * dUnz + \
                                        thJH_Un_verty * thUn_verty_Unx * dUnx + thJH_Un_verty * thUn_verty_Unz * dUnz + \
                                        thJH_Un_vertz * thUn_vertz_Unx * dUnx + thJH_Un_vertz * thUn_vertz_Unx * dUnx
                dJH_Ne[lat, lon, lev] = thJH_Ne * dNe * ccm ** 2
                # ######################################################################################################################
                # ############################ OHMIC HEATING ERROR ############################
                # #############################################################################
                # |theta_OH / theta_sigmaP|^2
                thOH_sigmaP = (Estar[0] ** 2 + Estar[1] ** 2 + Estar[2] ** 2) ** 2
                # |theta_OH /theta_Evertx|^2
                thOH_Evertx = (2 * pedersen_con[lat, lon, lev] * Estar[0]) ** 2
                # |theta_OH /theta_Everty|^2
                thOH_Everty = (2 * pedersen_con[lat, lon, lev] * Estar[1]) ** 2
                # |theta_OH /theta_Evertz|^2
                thOH_Evertz = (2 * pedersen_con[lat, lon, lev] * Estar[2]) ** 2
                # |theta_OH / thetaUn_vertx|^2
                thOH_Un_vertx = (2 * pedersen_con[lat, lon, lev] * (B[1] * Estar[2] - B[2] * Estar[1])) ** 2
                # |theta_OH / thetaUn_verty|^2
                thOH_Un_verty = (2 * pedersen_con[lat, lon, lev] * (B[2] * Estar[0] - B[0] * Estar[2])) ** 2
                # |theta_OH / thetaUn_vertz|^2
                thOH_Un_vertz = (2 * pedersen_con[lat, lon, lev] * (B[0] * Estar[1] - B[1] * Estar[0])) ** 2
                # |theta_OH / theta_Bx|^2
                thOH_Bx = (2 * pedersen_con[lat, lon, lev] * (Un[2] * Estar[1] - Un[1] * Estar[2])) ** 2
                # |theta_OH / theta_By|^2
                thOH_By = (2 * pedersen_con[lat, lon, lev] * (Un[0] * Estar[2] - Un[2] * Estar[0])) ** 2
                # |theta_OH / theta_Bz|^2
                thOH_Bz = (2 * pedersen_con[lat, lon, lev] * (Un[1] * Estar[0] - Un[0] * Estar[1])) ** 2

                # |dOH|
                Ohmic_Heating_error[lat, lon, lev] = np.sqrt(thOH_sigmaP * pedersen_con_error[lat, lon, lev] ** 2 + thOH_Evertx * dEvertx +
                                                             thOH_Everty * dEverty + thOH_Evertz * dEvertz + thOH_Un_vertx * dUn_vertx +
                                                             thOH_Un_verty * dUn_verty + thOH_Un_vertz * dUn_vertz + thOH_Bx * dBx + thOH_By * dBy +
                                                             thOH_Bz * dBz)
                # ######################################### Ohmic Heating contributions error #########################################
                # squared
                dOH_B[lat, lon, lev] = thOH_Bx * dBx + thOH_By * dBy + thOH_Bz * dBz + thOH_sigmaP * dsp_B[lat, lon, lev] + \
                                       thOH_Evertx * thEvertx_by * dby + thOH_Evertx * thEvertx_bz * dbz + thOH_Everty * thEverty_bx * dbx + \
                                       thOH_Everty * thEverty_bz * dbz + thOH_Evertz * thEvertz_bx * dbx + thOH_Evertz * thEvertz_by * dby + \
                                       thOH_Un_vertx * thUn_vertx_by * dby + thOH_Un_vertx * thUn_vertx_bz * dbz + \
                                       thOH_Un_verty * thUn_verty_bx * dbx + thOH_Un_verty * thUn_verty_bz * dbz + \
                                       thOH_Un_vertz * thUn_vertz_bx * dbx + thOH_Un_vertz * thUn_vertz_by * dby
                dOH_E[lat, lon, lev] = thOH_Evertx * thEvertx_Ey * dEy + thOH_Evertx * thEvertx_Ez * dEz + thOH_Everty * thEverty_Ex * dEx + \
                                       thOH_Everty * thEverty_Ez * dEz + thOH_Evertz * thEvertz_Ex * dEx + thOH_Evertz * thEvertz_Ey * dEy
                dOH_Un[lat, lon, lev] = thOH_Un_vertx * thUn_vertx_Uny * dUny + thOH_Un_vertx * thUn_vertx_Unz * dUnz + \
                                        thOH_Un_verty * thUn_verty_Unx * dUnx + thOH_Un_verty * thUn_verty_Unz * dUnz + \
                                        thOH_Un_vertz * thUn_vertz_Unx * dUnx + thOH_Un_vertz * thUn_vertz_Uny * dUny
                dOH_Nion[lat, lon, lev] = thOH_sigmaP * dsp_Nion[lat, lon, lev]
                dOH_Nneutral[lat, lon, lev] = thOH_sigmaP * dsp_Nneutral[lat, lon, lev]
                dOH_Ne[lat, lon, lev] = thOH_sigmaP * dsp_Ne[lat, lon, lev]
                dOH_Te[lat, lon, lev] = thOH_sigmaP * dsp_Te[lat, lon, lev]
                dOH_Ti[lat, lon, lev] = thOH_sigmaP * dsp_Ti[lat, lon, lev]
                dOH_Tn[lat, lon, lev] = thOH_sigmaP * dsp_Tn[lat, lon, lev]
                dOH_sp[lat, lon, lev] = thOH_sigmaP * pedersen_con_error[lat, lon, lev] ** 2
                # #####################################################################################################################
                # ############################ FRICTIONAL HEATING ERROR ############################
                # ##################################################################################
                # densities error in m^(-3)
                termFH_Op = mkO * nu_Op_sum[lat, lon, lev] * NOp[lat, lon, lev] * ccm
                termFH_O2p = mkO2 * nu_O2p_sum[lat, lon, lev] * NO2p[lat, lon, lev] * ccm
                termFH_NOp = mkNO * nu_NOp_sum[lat, lon, lev] * NNOp[lat, lon, lev] * ccm

                # ion velocity - neutral wind difference vector
                DV = [Vi_vert[0] - Un_vert[0], Vi_vert[1] - Un_vert[1], Vi_vert[2] - Un_vert[2]]

                # |theta_FH / theta_Vi_vertx|^2
                thFH_Vi_vertx = (2 * (DV[0]) * (termFH_Op + termFH_O2p + termFH_NOp)) ** 2
                # |theta_FH / theta_Vi_verty|^2
                thFH_Vi_verty = (2 * (DV[1]) * (termFH_Op + termFH_O2p + termFH_NOp)) ** 2
                # |theta_FH / theta_Vi_vertz|^2
                thFH_Vi_vertz = (2 * (DV[2]) * (termFH_Op + termFH_O2p + termFH_NOp)) ** 2
                # |theta_FH / theta_Un_vertx|^2
                thFH_Un_vertx = thFH_Vi_vertx
                # |theta_FH / theta_Un_verty|^2
                thFH_Un_verty = thFH_Vi_verty
                # |theta_FH / theta_Un_vertz|^2
                thFH_Un_vertz = thFH_Vi_vertz
                # |theta_FH / theta_NOp|^2
                thFH_NOp = (mkO * nu_Op_sum[lat, lon, lev] * (DV[0] ** 2 + DV[1] ** 2 + DV[2] ** 2)) ** 2
                # |theta_FH / theta_NO2p|^2
                thFH_NO2p = (mkO2 * nu_O2p_sum[lat, lon, lev] * (DV[0] ** 2 + DV[1] ** 2 + DV[2] ** 2)) ** 2
                # |theta_FH / theta_NNOp|^2
                thFH_NNOp = (mkNO * nu_NOp_sum[lat, lon, lev] * (DV[0] ** 2 + DV[1] ** 2 + DV[2] ** 2)) ** 2
                # |theta_FH / theta_nuOp|^2
                thFH_nuOp = (mkO * NOp[lat, lon, lev] * ccm * (DV[0] ** 2 + DV[1] ** 2 + DV[2] ** 2)) ** 2
                # |theta_FH / theta_nuO2p|^2
                thFH_nuO2p = (mkO2 * NO2p[lat, lon, lev] * ccm * (DV[0] ** 2 + DV[1] ** 2 + DV[2] ** 2)) ** 2
                # |theta_FH / theta_nuNOp|^2
                thFH_nuNOp = (mkNO * NNOp[lat, lon, lev] * ccm * (DV[0] ** 2 + DV[1] ** 2 + DV[2] ** 2)) ** 2

                # |dFH|
                Frictional_Heating_error[lat, lon, lev] = np.sqrt(thFH_Vi_vertx * dVix + thFH_Vi_verty * dViy + thFH_Vi_vertz * dViz +
                                                                  thFH_Un_vertx * dUn_vertx + thFH_Un_verty * dUn_verty + thFH_Un_vertz * dUn_vertz +
                                                                  thFH_NOp * dNOp * ccm ** 2 + thFH_NO2p * dNO2p * ccm ** 2 +
                                                                  thFH_NNOp * dNNOp * ccm ** 2 + thFH_nuOp * nuOp_error[lat, lon, lev] ** 2 +
                                                                  thFH_nuO2p * nuO2p_error[lat, lon, lev] ** 2 +
                                                                  thFH_nuNOp * nuNOp_error[lat, lon, lev] ** 2)
                # ############################# Frictional Heating contributions error #############################
                # squared
                dFH_B[lat, lon, lev] = thFH_Un_vertx * thUn_vertx_by * dby + thFH_Un_vertx * thUn_vertx_bz * dbz + \
                                       thFH_Un_verty * thUn_verty_bx * dbx + thFH_Un_verty * thUn_verty_bz * dbz + \
                                       thFH_Un_vertz * thUn_vertz_bx * dbx + thFH_Un_vertz * thUn_vertz_by * dby
                dFH_Un[lat, lon, lev] = thFH_Un_vertx * thUn_vertx_Uny * dUny + thFH_Un_vertx * thUn_vertx_Unz * dUnz + \
                                        thFH_Un_verty * thUn_verty_Unx * dUnx + thFH_Un_verty * thUn_verty_Unz * dUnz + \
                                        thFH_Un_vertz * thUn_vertz_Unx * dUnx + thFH_Un_vertz * thUn_vertz_Uny * dUny
                dFH_Vi[lat, lon, lev] = thFH_Vi_vertx * dVix + thFH_Vi_verty * dViy + thFH_Vi_vertz * dViz
                dFH_Nion[lat, lon, lev] = thFH_NOp * dNOp * ccm ** 2 + thFH_NO2p * dNO2p * ccm ** 2 + thFH_NNOp * dNNOp * ccm ** 2
                dFH_Nneutral[lat, lon, lev] = thFH_nuOp * (thOp_O_NO * dNO + thOp_O2_NO2 * dNO2 + thOp_N2_NN2 * dNN2) + \
                                              thFH_nuO2p * (thO2p_O_NO * dNO + thO2p_O2_NO2 * dNO2 + thO2p_N2_NN2 * dNN2) + \
                                              thFH_nuNOp * (thNOp_O_NO * dNO + thNOp_O2_NO2 * dNO2 + thNOp_N2_NN2 * dNN2)
                dFH_Ti[lat, lon, lev] = thFH_nuOp * thOp_O_Ti * dTi + thFH_nuO2p * thO2p_O2_Ti * dTi
                dFH_Tn[lat, lon, lev] = thFH_nuOp * thOp_O_Tn * dTn + thFH_nuO2p * thO2p_O2_Tn * dTn
                dFH_nu[lat, lon, lev] = thFH_nuOp * nuOp_error[lat, lon, lev] ** 2 + thFH_nuO2p * nuO2p_error[lat, lon, lev] ** 2 + \
                                        thFH_nuNOp * nuNOp_error[lat, lon, lev] ** 2
                # ####################################################################################################
                # ########################## CROSS SECTIONS ERROR ##########################
                # nu(in Hz), N(in m^(-3)), T(in kelvin), mass(in kg)
                N_neutral = NO[lat, lon, lev] + NO2[lat, lon, lev] + NN2[lat, lon, lev]
                N_neutral = N_neutral * ccm
                nu_ion = nu_Op_sum[lat, lon, lev] + nu_O2p_sum[lat, lon, lev] + nu_NOp_sum[lat, lon, lev]
                nu_ion = nu_ion / 3
                m_ion = mkO + mkO2 + mkNO
                m_ion = m_ion / 3
                # |dnu_ion|^2
                dnu_ion = (nuOp_error[lat, lon, lev] ** 2 + nuO2p_error[lat, lon, lev] ** 2 + nuNOp_error[lat, lon, lev] ** 2) / 3
                # |dN_neutral|^2
                dN_neutral = dNO * ccm ** 2 + dNO2 * ccm ** 2 + dNN2 * ccm ** 2
                # ############ O+ ###########
                # |theta_COp / theta_nuOp|^2
                thCOp_nuOp = (np.sqrt(mkO / (2 * boltzmann * Ti[lat, lon, lev])) / N_neutral) ** 2
                # |theta_COp / theta_N_neutral|^2
                thCOp_N_neutral = (nu_Op_sum[lat, lon, lev] * np.sqrt(mkO / (2 * boltzmann * Ti[lat, lon, lev])) / N_neutral ** 2) ** 2
                # |theta_COp / theta_Ti|^2
                thCOp_Ti = (nu_Op_sum[lat, lon, lev] * np.sqrt(mkO / (2 * boltzmann * Ti[lat, lon, lev])) / (2 * N_neutral * Ti[lat, lon, lev])) ** 2
                # |dCOp|
                C_Op_error[lat, lon, lev] = np.sqrt(thCOp_nuOp * nuOp_error[lat, lon, lev] ** 2 + thCOp_N_neutral * dN_neutral + thCOp_Ti * dTi)
                # ############ O2+ ############
                # |theta_CO2p / theta_nuO2p|^2
                thCO2p_nuO2p = (np.sqrt(mkO2 / (2 * boltzmann * Ti[lat, lon, lev])) / N_neutral) ** 2
                # |theta_CO2p / theta_N_neutral|^2
                thCO2p_N_neutral = (nu_O2p_sum[lat, lon, lev] * np.sqrt(mkO2 / (2 * boltzmann * Ti[lat, lon, lev])) / N_neutral ** 2) ** 2
                # |theta_CO2p / theta_Ti|^2
                thCO2p_Ti = (nu_O2p_sum[lat, lon, lev] * np.sqrt(mkO2 / (2 * boltzmann * Ti[lat, lon, lev])) /
                             (2 * N_neutral * Ti[lat, lon, lev])) ** 2
                # |dCOp|
                C_O2p_error[lat, lon, lev] = np.sqrt(thCO2p_nuO2p * nuO2p_error[lat, lon, lev] ** 2 + thCO2p_N_neutral * dN_neutral + thCO2p_Ti * dTi)
                # ############ NO+ ############
                # |theta_CNOp / theta_nuNOp|^2
                thCNOp_nuNOp = (np.sqrt(mkNO / (2 * boltzmann * Ti[lat, lon, lev])) / N_neutral) ** 2
                # |theta_CNOp / theta_N_neutral|^2
                thCNOp_N_neutral = (nu_NOp_sum[lat, lon, lev] * np.sqrt(mkNO / (2 * boltzmann * Ti[lat, lon, lev])) / N_neutral ** 2) ** 2
                # |theta_CNOp / theta_Ti|^2
                thCNOp_Ti = (nu_NOp_sum[lat, lon, lev] * np.sqrt(mkNO / (2 * boltzmann * Ti[lat, lon, lev])) /
                             (2 * N_neutral * Ti[lat, lon, lev])) ** 2
                # |dCNOp|
                C_NOp_error[lat, lon, lev] = np.sqrt(thCNOp_nuNOp * nuNOp_error[lat, lon, lev] ** 2 + thCNOp_N_neutral * dN_neutral + thCNOp_Ti * dTi)
                # ############ ion #############
                # |theta_Cion / theta_nu_ion|^2
                thCion_nuion = (np.sqrt(m_ion / (2 * boltzmann * Ti[lat, lon, lev])) / N_neutral) ** 2
                # |theta_Cion / theta_N_neutral|^2
                thCion_N_neutral = (nu_ion * np.sqrt(m_ion / (2 * boltzmann * Ti[lat, lon, lev])) / N_neutral ** 2) ** 2
                # |theta_Cion / theta_Ti|^2
                thCion_Ti = (nu_ion * np.sqrt(m_ion / (2 * boltzmann * Ti[lat, lon, lev])) / (2 * N_neutral * Ti[lat, lon, lev])) ** 2
                # |dCion|
                C_ion_error[lat, lon, lev] = np.sqrt(thCion_nuion * dnu_ion + thCion_N_neutral * dN_neutral + thCion_Ti * dTi)
                # ###################################### Cross section contributions error ######################################
                # squared
                dCion_Ti[lat, lon, lev] = thCion_Ti * dTi + thCion_nuion * (thOp_O_Ti * dTi + thO2p_O2_Ti * dTi) / 3
                dCion_Tn[lat, lon, lev] = thCion_nuion * (thOp_O_Tn * dTn + thO2p_O2_Tn * dTn) / 3
                dCion_nu[lat, lon, lev] = thCion_nuion * dnu_ion
                dCion_Nneutral[lat, lon, lev] = thCion_N_neutral * dN_neutral + thCion_nuion * (
                                                thOp_O_NO * dNO + thOp_O2_NO2 * dNO2 + thOp_N2_NN2 * dNN2 + thO2p_O_NO * dNO + thO2p_O2_NO2 * dNO2 +
                                                thO2p_N2_NN2 * dNN2 + thNOp_O_NO * dNO + thNOp_O2_NO2 * dNO2 + thNOp_N2_NN2 * dNN2) / 3
                # ###############################################################################################################
                # ########################## CURRENTS ERROR PROPAGATION ##########################
                # #################### 1st Methodology - Ohms law ####################
                # ########### Pedersen current error ###########
                mag = np.sqrt(Estar[0] ** 2 + Estar[1] ** 2 + Estar[2] ** 2)

                # |theta_JP / theta_sigmaP|^2
                thJp_sigmaP = mag ** 2
                # |theta_JP / theta_Evertx|^2
                thJp_Evertx = (pedersen_con[lat, lon, lev] * Estar[0] / mag) ** 2
                # |theta_JP / theta_Everty|^2
                thJp_Everty = (pedersen_con[lat, lon, lev] * Estar[1] / mag) ** 2
                # |theta_JP / theta_Evertz|^2
                thJp_Evertz = (pedersen_con[lat, lon, lev] * Estar[2] / mag) ** 2
                # |theta_JP / theta_Un_vertx|^2
                thJp_Un_vertx = (pedersen_con[lat, lon, lev] * (B[1] * Estar[2] - B[2] * Estar[1]) / mag) ** 2
                # |theta_JP / theta_Un_verty|^2
                thJp_Un_verty = (pedersen_con[lat, lon, lev] * (B[2] * Estar[0] - B[0] * Estar[2]) / mag) ** 2
                # |theta_JP / theta_Un_vertz|^2
                thJp_Un_vertz = (pedersen_con[lat, lon, lev] * (B[0] * Estar[1] - B[1] * Estar[0]) / mag) ** 2
                # |theta_JP / theta_Bx|^2
                thJp_Bx = (pedersen_con[lat, lon, lev] * (Un_vert[2] * Estar[1] - Un_vert[1] * Estar[2]) / mag) ** 2
                # |theta_JP / theta_By|^2
                thJp_By = (pedersen_con[lat, lon, lev] * (Un_vert[0] * Estar[2] - Un_vert[2] * Estar[0]) / mag) ** 2
                # |theta_JP / theta_Bz|^2
                thJp_Bz = (pedersen_con[lat, lon, lev] * (Un_vert[1] * Estar[0] - Un_vert[0] * Estar[1]) / mag) ** 2

                # |dJP|
                J_pedersen_error[lat, lon, lev] = np.sqrt(thJp_sigmaP * pedersen_con_error[lat, lon, lev] ** 2 + thJp_Evertx * dEvertx +
                                                          thJp_Everty * dEverty + thJp_Evertz * dEvertz + thJp_Un_vertx * dUn_vertx +
                                                          thJp_Un_verty * dUn_verty + thJp_Un_vertz * dUn_vertz + thJp_Bx * dBx + thJp_By * dBy +
                                                          thJp_Bz * dBz)
                # ############# Hall current error #############
                x = bunit[1] * Estar[2] - bunit[2] * Estar[1]
                y = bunit[2] * Estar[0] - bunit[0] * Estar[2]
                z = bunit[0] * Estar[1] - bunit[1] * Estar[0]
                mag1 = np.sqrt(x ** 2 + y ** 2 + z ** 2)

                # |theta_JH / theta_sigmaH|^2
                thJh_sigmaH = mag1 ** 2
                # |theta_JH / theta_Evertx|^2
                thJh_Evertx = (hall_con[lat, lon, lev] * (bunit[2] * y - bunit[1] * z) / mag1) ** 2
                # |theta_JH / theta_Everty|^2
                thJh_Everty = (hall_con[lat, lon, lev] * (bunit[0] * z - bunit[2] * x) / mag1) ** 2
                # |theta_JH / theta_Evertz|^2
                thJh_Evertz = (hall_con[lat, lon, lev] * (bunit[1] * x - bunit[0] * y) / mag1) ** 2
                # |theta_JH / theta_bx|^2
                thJh_bx = (hall_con[lat, lon, lev] * (Estar[1] * z - Estar[2] * y) / mag1) ** 2
                # |theta_JH / theta_by|^2
                thJh_by = (hall_con[lat, lon, lev] * (Estar[2] * x - Estar[0] * z) / mag1) ** 2
                # |theta_JH / theta_bz|^2
                thJh_bz = (hall_con[lat, lon, lev] * (Estar[0] * y - Estar[1] * x) / mag1) ** 2
                # |theta_JH / theta_Bx|^2
                thJh_Bx = (hall_con[lat, lon, lev] * (- x * (bunit[1] * Un_vert[1] + bunit[2] * Un_vert[2]) +
                                                      bunit[0] * (Un_vert[1] * y + Un_vert[2] * z)) / mag1) ** 2
                # |theta_JH / theta_By|^2
                thJh_By = (hall_con[lat, lon, lev] * (- y * (bunit[0] * Un_vert[0] + bunit[2] * Un_vert[2]) +
                                                      bunit[1] * (Un_vert[0] * x + Un_vert[2] * z)) / mag1) ** 2
                # |theta_JH / theta_Bz|^2
                thJh_Bz = (hall_con[lat, lon, lev] * (- z * (bunit[0] * Un_vert[0] + bunit[1] * Un_vert[1]) +
                                                      bunit[2] * (Un_vert[0] * x + Un_vert[1] * y)) / mag1) ** 2
                # |theta_JH / theta_Un_vertx|^2
                thJh_Un_vertx = (hall_con[lat, lon, lev] * (x * (bunit[1] * B[1] + bunit[2] * B[2]) - bunit[0] * (B[1] * y + B[2] * z)) / mag1) ** 2
                # |theta_JH / theta_Un_verty|^2
                thJh_Un_verty = (hall_con[lat, lon, lev] * (y * (bunit[0] * B[0] + bunit[2] * B[2]) - bunit[1] * (B[0] * x + B[2] * z)) / mag1) ** 2
                # |theta_JH / theta_Un_vertz|^2
                thJh_Un_vertz = (hall_con[lat, lon, lev] * (z * (bunit[0] * B[0] + bunit[1] * B[1]) - bunit[2] * (B[0] * x + B[1] * y)) / mag1) ** 2

                # |dJh|
                J_hall_error[lat, lon, lev] = np.sqrt(thJh_sigmaH * hall_con_error[lat, lon, lev] ** 2 + thJh_Evertx * dEvertx +
                                                      thJh_Everty * dEverty + thJh_Evertz * dEvertz + thJh_bx * dbx + thJh_by * dby + thJh_bz * dbz +
                                                      thJh_Bx * dBx + thJh_By * dBy + thJh_Bz * dBz + thJh_Un_vertx * dUn_vertx +
                                                      thJh_Un_verty * dUn_verty + thJh_Un_vertz * dUn_vertz)

                # ############################## TOTAL CURRENT ERROR J_OHMIC ##############################
                x1 = pedersen_con[lat, lon, lev] * Estar[0] + hall_con[lat, lon, lev] * x
                y1 = pedersen_con[lat, lon, lev] * Estar[1] + hall_con[lat, lon, lev] * y
                z1 = pedersen_con[lat, lon, lev] * Estar[2] + hall_con[lat, lon, lev] * z
                mag2 = np.sqrt(x1 ** 2 + y1 ** 2 + z1 ** 2)

                # |theta_Johm / theta_sigmaP|^2
                thJohm_sigmaP = ((x1 * Estar[0] + y1 * Estar[1] + z1 * Estar[2]) / mag2) ** 2
                # |theta_Johm / theta_sigmaH|^2
                thJohm_sigmaH = ((x1 * x + y1 * y + z1 * z) / mag2) ** 2
                # |theta_Johm / theta_Evertx|^2
                thJohm_Evertx = ((x1 * pedersen_con[lat, lon, lev] + hall_con[lat, lon, lev] * (y1 * bunit[2] - z1 * bunit[1])) / mag2) ** 2
                # |theta_Johm / theta_Everty|^2
                thJohm_Everty = ((y1 * pedersen_con[lat, lon, lev] + hall_con[lat, lon, lev] * (z1 * bunit[0] - x1 * bunit[2])) / mag2) ** 2
                # |theta_Johm / theta_Evertz|^2
                thJohm_Evertz = ((z1 * pedersen_con[lat, lon, lev] + hall_con[lat, lon, lev] * (x1 * bunit[1] - y1 * bunit[0])) / mag2) ** 2
                # |theta_Johm / theta_Un_vertx|^2
                thJohm_Un_vertx = ((x1 * hall_con[lat, lon, lev] * (B[1] * bunit[1] + B[2] * bunit[2]) -
                                    y1 * (pedersen_con[lat, lon, lev] * B[2] + hall_con[lat, lon, lev] * B[1] * bunit[0]) +
                                    z1 * (pedersen_con[lat, lon, lev] * B[1] - hall_con[lat, lon, lev] * B[2] * bunit[0])) / mag2) ** 2
                # |theta_Johm / theta_Un_verty|^2
                thJohm_Un_verty = ((y1 * hall_con[lat, lon, lev] * (B[2] * bunit[2] + B[0] * bunit[0]) -
                                    z1 * (pedersen_con[lat, lon, lev] * B[0] + hall_con[lat, lon, lev] * B[2] * bunit[1]) +
                                    x1 * (pedersen_con[lat, lon, lev] * B[2] - hall_con[lat, lon, lev] * B[0] * bunit[1])) / mag2) ** 2
                # |theta_Johm / theta_Un_vertz|^2
                thJohm_Un_vertz = ((z1 * hall_con[lat, lon, lev] * (B[0] * bunit[0] + B[1] * bunit[1]) -
                                    x1 * (pedersen_con[lat, lon, lev] * B[1] + hall_con[lat, lon, lev] * B[1] * bunit[2]) +
                                    y1 * (pedersen_con[lat, lon, lev] * B[0] - hall_con[lat, lon, lev] * B[1] * bunit[2])) / mag2) ** 2
                # |theta_Johm / theta_bx|^2
                thJohm_bx = ((hall_con[lat, lon, lev] * (z1 * Estar[1] - y1 * Estar[2])) / mag2) ** 2
                # |theta_Johm / theta_by|^2
                thJohm_by = ((hall_con[lat, lon, lev] * (x1 * Estar[2] - z1 * Estar[0])) / mag2) ** 2
                # |theta_Johm / theta_bz|^2
                thJohm_bz = ((hall_con[lat, lon, lev] * (y1 * Estar[0] - x1 * Estar[1])) / mag2) ** 2
                # |theta_Johm / theta_Bx|^2
                thJohm_Bx = ((- x1 * hall_con[lat, lon, lev] * (Un_vert[1] * bunit[1] + Un_vert[2] * bunit[2]) +
                              y1 * (pedersen_con[lat, lon, lev] * Un_vert[2] + hall_con[lat, lon, lev] * Un_vert[1] * bunit[0]) +
                              z1 * (hall_con[lat, lon, lev] * Un_vert[2] * bunit[0] - pedersen_con[lat, lon, lev] * Un_vert[1])) / mag2) ** 2
                # |theta_Johm / theta_By|^2
                thJohm_By = ((- y1 * hall_con[lat, lon, lev] * (Un_vert[2] * bunit[2] + Un_vert[0] * bunit[0]) +
                              z1 * (pedersen_con[lat, lon, lev] * Un_vert[0] + hall_con[lat, lon, lev] * Un_vert[2] * bunit[1]) +
                              x1 * (hall_con[lat, lon, lev] * Un_vert[0] * bunit[1] - pedersen_con[lat, lon, lev] * Un_vert[2])) / mag2) ** 2
                # |theta_Johm / theta_Bz|^2
                thJohm_Bz = ((- z1 * hall_con[lat, lon, lev] * (Un_vert[0] * bunit[0] + Un_vert[1] * bunit[1]) +
                              x1 * (pedersen_con[lat, lon, lev] * Un_vert[1] + hall_con[lat, lon, lev] * Un_vert[0] * bunit[2]) +
                              y1 * (hall_con[lat, lon, lev] * Un_vert[1] * bunit[2] - pedersen_con[lat, lon, lev] * Un_vert[0])) / mag2) ** 2

                # |dJohmic|
                J_ohmic_error[lat, lon, lev] = np.sqrt(thJohm_sigmaP * pedersen_con_error[lat, lon, lev] ** 2 +
                                                       thJohm_sigmaH * hall_con_error[lat, lon, lev] ** 2 + thJohm_Evertx * dEvertx +
                                                       thJohm_Everty * dEverty + thJohm_Evertz * dEvertz + thJohm_Un_vertx * dUn_vertx +
                                                       thJohm_Un_verty * dUn_verty + thJohm_Un_vertz * dUn_vertz + thJohm_bx * dbx + thJohm_by * dby +
                                                       thJohm_bz * dbz + thJohm_Bx * dBx + thJohm_By * dBy + thJohm_Bz * dBz)
                # ######################################### Johmic contributions error #########################################
                # squared
                dJohm_B[lat, lon, lev] = thJohm_bx * dbx + thJohm_by * dby + thJohm_bz * dbz + thJohm_Bx * dBx + thJohm_By * dBy + thJohm_Bz * dBz + \
                                         thJohm_sigmaP * dsp_B[lat, lon, lev] + thJohm_sigmaH * dsh_B[lat, lon, lev] + \
                                         thJohm_Evertx * thEvertx_by * dby + thJohm_Evertx * thEvertx_bz * dbz + thJohm_Everty * thEverty_bx * dbx + \
                                         thJohm_Everty * thEverty_bz * dbz + thJohm_Evertz * thEvertz_bx * dbx + thJohm_Evertz * thEvertz_by * dby + \
                                         thJohm_Un_vertx * thUn_vertx_by * dby + thJohm_Un_vertx * thUn_vertx_bz * dbz + \
                                         thJohm_Un_verty * thUn_verty_bx * dbx + thJohm_Un_verty * thUn_verty_bz * dbz + \
                                         thJohm_Un_vertz * thUn_vertz_bx * dbx + thJohm_Un_vertz * thUn_vertz_by * dby
                dJohm_E[lat, lon, lev] = thJohm_Evertx * thEvertx_Ey * dEy + thJohm_Evertx * thEvertx_Ez * dEz + thJohm_Everty * thEverty_Ex * dEx + \
                                         thJohm_Everty * thEverty_Ez * dEz + thJohm_Evertz * thEvertz_Ex * dEx + thJohm_Evertz * thEvertz_Ey * dEy
                dJohm_Un[lat, lon, lev] = thJohm_Un_vertx * thUn_vertx_Uny * dUny + thJohm_Un_vertx * thUn_vertx_Unz * dUnz + \
                                          thJohm_Un_verty * thUn_verty_Unx * dUnx + thJohm_Un_verty * thUn_verty_Unz * dUnz + \
                                          thJohm_Un_vertz * thUn_vertz_Unx * dUnx + thJohm_Un_vertz * thUn_vertz_Uny * dUny
                dJohm_sp[lat, lon, lev] = thJohm_sigmaP * pedersen_con_error[lat, lon, lev] ** 2
                dJohm_sh[lat, lon, lev] = thJohm_sigmaH * hall_con_error[lat, lon, lev] ** 2
                dJohm_Ti[lat, lon, lev] = thJohm_sigmaP * dsp_Ti[lat, lon, lev] + thJohm_sigmaH * dsh_Ti[lat, lon, lev]
                dJohm_Tn[lat, lon, lev] = thJohm_sigmaP * dsp_Tn[lat, lon, lev] + thJohm_sigmaH * dsh_Tn[lat, lon, lev]
                dJohm_Te[lat, lon, lev] = thJohm_sigmaP * dsp_Te[lat, lon, lev] + thJohm_sigmaH * dsh_Te[lat, lon, lev]
                dJohm_Ne[lat, lon, lev] = thJohm_sigmaP * dsp_Ne[lat, lon, lev] + thJohm_sigmaH * dsh_Ne[lat, lon, lev]
                dJohm_Nneutral[lat, lon, lev] = thJohm_sigmaP * dsp_Nneutral[lat, lon, lev] + thJohm_sigmaH * dsh_Nneutral[lat, lon, lev]
                dJohm_Nion[lat, lon, lev] = thJohm_sigmaP * dsp_Nion[lat, lon, lev] + thJohm_sigmaH * dsh_Nion[lat, lon, lev]
                # ##############################################################################################################
                # #################### 2nd Methodology - Densities (current definition) ####################
                x2 = Vi_vert[0] - Un_vert[0] + (Estar[2] * B[1] - Estar[1] * B[2]) / Bnorm ** 2
                y2 = Vi_vert[1] - Un_vert[1] + (Estar[0] * B[2] - Estar[2] * B[0]) / Bnorm ** 2
                z2 = Vi_vert[2] - Un_vert[2] + (Estar[1] * B[0] - Estar[0] * B[1]) / Bnorm ** 2
                mag3 = np.sqrt(x2 ** 2 + y2 ** 2 + z2 ** 2)

                # |theta_Jd / theta_Vi_vertx|^2
                thJd_Vi_vertx = (qe * Ne[lat, lon, lev] * ccm * x2 / mag3) ** 2
                # |theta_Jd / theta_Vi_verty|^2
                thJd_Vi_verty = (qe * Ne[lat, lon, lev] * ccm * y2 / mag3) ** 2
                # |theta_Jd / theta_Vi_vertz|^2
                thJd_Vi_vertz = (qe * Ne[lat, lon, lev] * ccm * z2 / mag3) ** 2
                # |theta_Jd / theta_Un_vertx|^2
                thJd_Un_vertx = (qe * Ne[lat, lon, lev] * ccm * (x2 * ((B[2] ** 2 + B[1] ** 2) / Bnorm ** 2 - 1) - y2 * B[0] * B[1] / Bnorm ** 2 -
                                                                 z2 * B[0] * B[2] / Bnorm ** 2) / mag3) ** 2
                # |theta_Jd / theta_Un_verty|^2
                thJd_Un_verty = (qe * Ne[lat, lon, lev] * ccm * (y2 * ((B[0] ** 2 + B[2] ** 2) / Bnorm ** 2 - 1) - x2 * B[0] * B[1] / Bnorm ** 2 -
                                                                 z2 * B[1] * B[2] / Bnorm ** 2) / mag3) ** 2
                # |theta_Jd / theta_Un_vertz|^2
                thJd_Un_vertz = (qe * Ne[lat, lon, lev] * ccm * (z2 * ((B[0] ** 2 + B[1] ** 2) / Bnorm ** 2 - 1) - x2 * B[0] * B[2] / Bnorm ** 2 -
                                                                 y2 * B[1] * B[2] / Bnorm ** 2) / mag3) ** 2
                # |theta_Jd / theta_Evertx|^2
                thJd_Evertx = (qe * Ne[lat, lon, lev] * ccm * (y2 * B[2] / Bnorm ** 2 - z2 * B[1] / Bnorm ** 2) / mag3) ** 2
                # |theta_Jd / theta_Everty|^2
                thJd_Everty = (qe * Ne[lat, lon, lev] * ccm * (z2 * B[0] / Bnorm ** 2 - x2 * B[2] / Bnorm ** 2) / mag3) ** 2
                # |theta_Jd / theta_Evertz|^2
                thJd_Evertz = (qe * Ne[lat, lon, lev] * ccm * (x2 * B[1] / Bnorm ** 2 - y2 * B[0] / Bnorm ** 2) / mag3) ** 2
                # |theta_Jd / theta_Bx|^2
                thJd_Bx = (qe * Ne[lat, lon, lev] * ccm * (- x2 * (Un_vert[2] * B[2] + Un_vert[1] * B[1]) / Bnorm ** 2 +
                                                             y2 * (Un_vert[1] * B[0] - Estar[2]) / Bnorm ** 2 +
                                                             z2 * (Un_vert[2] * B[0] + Estar[1])) / mag3) ** 2
                # |theta_Jd / theta_By|^2
                thJd_By = (qe * Ne[lat, lon, lev] * ccm * (- y2 * (Un_vert[0] * B[0] + Un_vert[2] * B[2]) / Bnorm ** 2 +
                                                             z2 * (Un_vert[2] * B[1] - Estar[0]) / Bnorm ** 2 +
                                                             x2 * (Un_vert[0] * B[1] + Estar[2])) / mag3) ** 2
                # |theta_Jd / theta_Bz|^2
                thJd_Bz = (qe * Ne[lat, lon, lev] * ccm * (- z2 * (Un_vert[1] * B[1] + Un_vert[0] * B[0]) / Bnorm ** 2 +
                                                             x2 * (Un_vert[0] * B[1] - Estar[1]) / Bnorm ** 2 +
                                                             y2 * (Un_vert[1] * B[2] + Estar[0])) / mag3) ** 2
                # |theta_Jd / theta_B|^2
                thJd_B = (qe * Ne[lat, lon, lev] * ccm * (2 * x2 * (Estar[1] * B[2] - Estar[2] * B[1]) / Bnorm ** 3 +
                                                          2 * y2 * (Estar[2] * B[0] - Estar[0] * B[2]) / Bnorm ** 3 +
                                                          2 * z2 * (Estar[0] * B[1] - Estar[1] * B[0]) / Bnorm ** 3) / mag3) ** 2
                # |theta_Jd / theta_Ne|^2
                thJd_Ne = (qe * mag3) ** 2

                # |dJd|
                J_dens_error[lat, lon, lev] = np.sqrt(thJd_Vi_vertx * dVix + thJd_Vi_verty * dViy + thJd_Vi_vertz * dViz + thJd_Un_vertx * dUn_vertx +
                                                      thJd_Un_verty * dUn_verty + thJd_Un_vertz * dUn_vertz + thJd_Evertx * dEvertx +
                                                      thJd_Everty * dEverty + thJd_Evertz * dEvertz + thJd_Bx * dBx + thJd_By * dBy + thJd_Bz * dBz +
                                                      thJd_B * dB + thJd_Ne * dNe * ccm ** 2)
                # ###################################### J(densities) contributions error ######################################
                # squared
                dJd_B[lat, lon, lev] = thJd_Bx * dBx + thJd_By * dBy + thJd_Bz * dBz + thJd_B * dB + thJd_Un_vertx * thUn_vertx_by * dby + \
                                       thJd_Un_vertx * thUn_vertx_bz * dbz + thJd_Un_verty * thUn_verty_bx * dbx + \
                                       thJd_Un_verty * thUn_verty_bz * dbz + thJd_Un_vertz * thUn_vertz_bx * dbx + \
                                       thJd_Un_vertz * thUn_vertz_by * dby + thJd_Evertx * thEvertx_by * dby + thJd_Evertx * thEvertx_bz * dbz + \
                                       thJd_Everty * thEverty_bx * dbx + thJd_Everty * thEverty_bz * dbz + thJd_Evertz * thEvertz_bx * dbx + \
                                       thJd_Evertz * thEvertz_by * dby
                dJd_E[lat, lon, lev] = thJd_Evertx * thEvertx_Ey * dEy + thJd_Evertx * thEvertx_Ez * dEz + thJd_Everty * thEverty_Ex * dEx + \
                                       thJd_Everty * thEverty_Ez * dEz + thJd_Evertz * thEvertz_Ex * dEx + thJd_Evertz * thEvertz_Ey * dEy
                dJd_Vi[lat, lon, lev] = thJd_Vi_vertx * dVix + thJd_Vi_verty * dViy + thJd_Vi_vertz * dViz
                dJd_Un[lat, lon, lev] = thJd_Un_vertx * thUn_vertx_Uny * dUny + thJd_Un_vertx * thUn_vertx_Unz * dUnz + \
                                        thJd_Un_verty * thUn_verty_Unx * dUnx + thJd_Un_verty * thUn_verty_Unz * dUnz + \
                                        thJd_Un_vertz * thUn_vertz_Unx * dUnx + thJd_Un_vertz * thUn_vertz_Uny * dUny
                dJd_Ne[lat, lon, lev] = thJd_Ne * dNe * ccm ** 2
                # ##############################################################################################################
                progress_bar.UpdateBar(i)

    time.sleep(3)
    window.close()
    Sg.popup("_" * 50, "Error calculated in : " + str(time.time() - start_time) + " sec!", "_" * 50, title="Finished", keep_on_top=True)
    print('Calculated Errors in: ', time.time() - start_time, ' sec !')
    print(' ')
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END OF ERROR CALCULATION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ PLOTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ############################ Vertical Profile Plots ############################
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
    fig1.add_trace(go.Scatter(x=Omega_ion[lat, lon, :-1], y=heights[lat, lon, :-1], name="$Ω_{i}$", mode='lines',
                              line=dict(shape='spline', color='brown')))
    fig1.add_trace(go.Scatter(x=Omega_e[lat, lon, :-1], y=heights[lat, lon, :-1], name="$Ω_{e}$", mode='lines',
                              line=dict(shape='spline', color='black')))

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


def plot_collisions_error(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=nuOp_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="$ν_{O^{+}}$ error", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=nuO2p_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="$ν_{O2^{+}}$ error", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=nuNOp_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="$ν_{NO^{+}}$ error", mode='lines',
                             line=dict(shape='spline', color='green')))
    fig.add_trace(go.Scatter(x=nue_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="$ν_{e}$ error", mode='lines',
                             line=dict(shape='spline', color='purple')))

    # updating the layout of the figure
    fig.update_layout(xaxis_type="log", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                      xaxis_title="$Frequency \ (Hz)$", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Collision Frequencies Absolute Error' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

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

    nuOp_rel = nuOp_error[lat, lon, :-1] / nu_Op_sum[lat, lon, :-1]
    nuO2p_rel = nuO2p_error[lat, lon, :-1] / nu_O2p_sum[lat, lon, :-1]
    nuNOp_rel = nuNOp_error[lat, lon, :-1] / nu_NOp_sum[lat, lon, :-1]
    nu_ion = (nu_Op_sum[lat, lon, :-1] + nu_O2p_sum[lat, lon, :-1] + nu_NOp_sum[lat, lon, :-1]) / 3
    nuion_error = (nuOp_error[lat, lon, :-1] + nuO2p_error[lat, lon, :-1] + nuNOp_error[lat, lon, :-1]) / 3
    nuion_rel = nuion_error / nu_ion
    nue_rel = nue_error[lat, lon, :-1] / nu_e_sum[lat, lon, :-1]

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


def plot_collisions_contr(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    nu_ion = (nu_Op_sum[lat, lon, :-1] + nu_O2p_sum[lat, lon, :-1] + nu_NOp_sum[lat, lon, :-1]) / 3
    nuion_error = (nuOp_error[lat, lon, :-1] + nuO2p_error[lat, lon, :-1] + nuNOp_error[lat, lon, :-1]) / 3
    nuion_rel = nuion_error / nu_ion
    nue_rel = nue_error[lat, lon, :-1] / nu_e_sum[lat, lon, :-1]

    dnui_Tn = dnuion_Tn[lat, lon, :-1] ** (1/2) / nu_ion
    dnui_Ti = dnuion_Ti[lat, lon, :-1] ** (1/2) / nu_ion
    dnui_Nn = dnuion_Nneutral[lat, lon, :-1] ** (1 / 2) / nu_ion

    dnu_e_Te = dnue_Te[lat, lon, :-1] ** (1/2) / nu_e_sum[lat, lon, :-1]
    dnu_e_Nneutral = dnue_Nneutral[lat, lon, :-1] ** (1/2) / nu_e_sum[lat, lon, :-1]

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=nuion_rel, y=heights[lat, lon, :-1], name="νion error", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=dnui_Tn, y=heights[lat, lon, :-1], name="dTn(i)", mode='lines',
                             line=dict(shape='spline', dash="dot", color='red')))
    fig.add_trace(go.Scatter(x=dnui_Ti, y=heights[lat, lon, :-1], name="dTi(i)", mode='lines',
                             line=dict(shape='spline', dash="dash", color='red')))
    fig.add_trace(go.Scatter(x=dnui_Nn, y=heights[lat, lon, :-1], name="dNn(i)", mode='lines',
                             line=dict(shape='spline', dash="dot", color='brown')))

    fig.add_trace(go.Scatter(x=nue_rel, y=heights[lat, lon, :-1], name="νe error", mode='lines', line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=dnu_e_Te, y=heights[lat, lon, :-1], name="dTe(e)", mode='lines', line=dict(shape='spline', dash="dot", color='blue')))
    fig.add_trace(go.Scatter(x=dnu_e_Nneutral, y=heights[lat, lon, :-1], name="dNn(e)", mode='lines',
                             line=dict(shape='spline', dash="dash", color='blue')))

    # updating the layout of the figure
    fig.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)), xaxis_title="",
                      yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Collision Frequencies Relative Error Contributions' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center',
                             'yanchor': 'top'})

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
    fig.add_trace(go.Scatter(x=Ohmic_Heating_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="Ohmic Heating error", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=Frictional_Heating_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="Frictional Heating error", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=Joule_Heating_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="Joule Heating error", mode='lines',
                             line=dict(shape='spline', color='green')))

    # updating the layout of the figure
    fig.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                      xaxis_title="$(W/m^{3})$", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Heating Rates Absolute Error' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

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

    Ohmic_rel = Ohmic_Heating_error[lat, lon, :-1] / Ohmic_Heating[lat, lon, :-1]
    Frict_rel = Frictional_Heating_error[lat, lon, :-1] / Frictional_Heating[lat, lon, :-1]
    Joule_rel = Joule_Heating_error[lat, lon, :-1] / Joule_Heating[lat, lon, :-1]

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


def plot_heating_rates_contr(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    # Ohmic Heating
    Ohmic_rel = Ohmic_Heating_error[lat, lon, :-1] / Ohmic_Heating[lat, lon, :-1]
    dOhm_dB = dOH_B[lat, lon, :-1] ** (1/2) / Ohmic_Heating[lat, lon, :-1]
    dOhm_dE = dOH_E[lat, lon, :-1] ** (1/2) / Ohmic_Heating[lat, lon, :-1]
    dOhm_dNneutral = dOH_Nneutral[lat, lon, :-1] ** (1/2) / Ohmic_Heating[lat, lon, :-1]
    dOhm_dNion = dOH_Nion[lat, lon, :-1] ** (1/2) / Ohmic_Heating[lat, lon, :-1]
    dOhm_dUn = dOH_Un[lat, lon, :-1] ** (1/2) / Ohmic_Heating[lat, lon, :-1]
    dOhm_dNe = dOH_Ne[lat, lon, :-1] ** (1/2) / Ohmic_Heating[lat, lon, :-1]
    dOhm_dTe = dOH_Te[lat, lon, :-1] ** (1/2) / Ohmic_Heating[lat, lon, :-1]
    dOhm_dTi = dOH_Ti[lat, lon, :-1] ** (1/2) / Ohmic_Heating[lat, lon, :-1]
    dOhm_dTn = dOH_Tn[lat, lon, :-1] ** (1/2) / Ohmic_Heating[lat, lon, :-1]
    dOhm_dsp = dOH_sp[lat, lon, :-1] ** (1/2) / Ohmic_Heating[lat, lon, :-1]

    fig1 = go.Figure()

    # adding the various plots
    fig1.add_trace(go.Scatter(x=Ohmic_rel, y=heights[lat, lon, :-1], name="Ohmic Heating error", mode='lines',
                              line=dict(shape='spline', color='red')))
    fig1.add_trace(go.Scatter(x=dOhm_dB, y=heights[lat, lon, :-1], name="dB", mode='lines',
                              line=dict(shape='spline', dash="dot", color='red')))
    fig1.add_trace(go.Scatter(x=dOhm_dE, y=heights[lat, lon, :-1], name="dE", mode='lines',
                              line=dict(shape='spline', dash="dash", color='red')))
    fig1.add_trace(go.Scatter(x=dOhm_dNneutral, y=heights[lat, lon, :-1], name="dNn", mode='lines',
                              line=dict(shape='spline', dash="dash", color='lightcoral')))
    fig1.add_trace(go.Scatter(x=dOhm_dNion, y=heights[lat, lon, :-1], name="dNion", mode='lines',
                              line=dict(shape='spline', dash="dot", color='red')))
    fig1.add_trace(go.Scatter(x=dOhm_dUn, y=heights[lat, lon, :-1], name="dUn", mode='lines',
                              line=dict(shape='spline', dash="dot", color='maroon')))
    fig1.add_trace(go.Scatter(x=dOhm_dNe, y=heights[lat, lon, :-1], name="dNe", mode='lines',
                              line=dict(shape='spline', dash="dash", color='maroon')))
    fig1.add_trace(go.Scatter(x=dOhm_dTe, y=heights[lat, lon, :-1], name="dTe", mode='lines',
                              line=dict(shape='spline', dash="dot", color='sienna')))
    fig1.add_trace(go.Scatter(x=dOhm_dTi, y=heights[lat, lon, :-1], name="dTi", mode='lines',
                              line=dict(shape='spline', dash="dash", color='sienna')))
    fig1.add_trace(go.Scatter(x=dOhm_dTn, y=heights[lat, lon, :-1], name="dTn", mode='lines',
                              line=dict(shape='spline', dash="dot", color='rosybrown')))
    fig1.add_trace(go.Scatter(x=dOhm_dsp, y=heights[lat, lon, :-1], name="dσPedersen", mode='lines',
                              line=dict(shape='spline', dash="dash", color='rosybrown')))

    # updating the layout of the figure
    fig1.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                       tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                       xaxis_title="$(W/m^{3})$", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                       title={'text': 'Ohmic Heating Error Contributions' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig1.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig1.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig1.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig1.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig1.show()

    # Frictional Heating
    Frict_rel = Frictional_Heating_error[lat, lon, :-1] / Frictional_Heating[lat, lon, :-1]
    dFric_dB = dFH_B[lat, lon, :-1] ** (1/2) / Frictional_Heating[lat, lon, :-1]
    dFric_dNneutral = dFH_Nneutral[lat, lon, :-1] ** (1/2) / Frictional_Heating[lat, lon, :-1]
    dFric_dNion = dFH_Nion[lat, lon, :-1] ** (1/2) / Frictional_Heating[lat, lon, :-1]
    dFric_dUn = dFH_Un[lat, lon, :-1] ** (1/2) / Frictional_Heating[lat, lon, :-1]
    dFric_dVi = dFH_Vi[lat, lon, :-1] ** (1/2) / Frictional_Heating[lat, lon, :-1]
    dFric_dTn = dFH_Tn[lat, lon, :-1] ** (1/2) / Frictional_Heating[lat, lon, :-1]
    dFric_dTi = dFH_Ti[lat, lon, :-1] ** (1/2) / Frictional_Heating[lat, lon, :-1]
    dFric_dnu = dFH_nu[lat, lon, :-1] ** (1/2) / Frictional_Heating[lat, lon, :-1]

    fig2 = go.Figure()

    fig2.add_trace(go.Scatter(x=Frict_rel, y=heights[lat, lon, :-1], name="Frictional Heating error", mode='lines',
                              line=dict(shape='spline', color='blue')))
    fig2.add_trace(go.Scatter(x=dFric_dB, y=heights[lat, lon, :-1], name="dB", mode='lines',
                              line=dict(shape='spline', dash="dot", color='blue')))
    fig2.add_trace(go.Scatter(x=dFric_dNneutral, y=heights[lat, lon, :-1], name="dNn", mode='lines',
                              line=dict(shape='spline', dash="dash", color='blue')))
    fig2.add_trace(go.Scatter(x=dFric_dNion, y=heights[lat, lon, :-1], name="dNion", mode='lines',
                              line=dict(shape='spline', dash="dot", color='turquoise')))
    fig2.add_trace(go.Scatter(x=dFric_dUn, y=heights[lat, lon, :-1], name="dUn", mode='lines',
                              line=dict(shape='spline', dash="dash", color='turquoise')))
    fig2.add_trace(go.Scatter(x=dFric_dVi, y=heights[lat, lon, :-1], name="dVi", mode='lines',
                              line=dict(shape='spline', dash="dot", color='darkslategrey')))
    fig2.add_trace(go.Scatter(x=dFric_dTn, y=heights[lat, lon, :-1], name="dTn", mode='lines',
                              line=dict(shape='spline', dash="dash", color='darkslategrey')))
    fig2.add_trace(go.Scatter(x=dFric_dTi, y=heights[lat, lon, :-1], name="dTi", mode='lines',
                              line=dict(shape='spline', dash="dot", color='steelblue')))
    fig2.add_trace(go.Scatter(x=dFric_dnu, y=heights[lat, lon, :-1], name="dν", mode='lines',
                              line=dict(shape='spline', dash="dash", color='steelblue')))
    # updating the layout of the figure
    fig2.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                       tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                       xaxis_title="", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                       title={'text': 'Frictional Heating Error Contributions' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig2.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig2.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig2.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig2.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig2.show()

    # Joule Heating
    Joule_rel = Joule_Heating_error[lat, lon, :-1] / Joule_Heating[lat, lon, :-1]
    dJoule_dB = dJH_B[lat, lon, :-1] ** (1/2) / Joule_Heating[lat, lon, :-1]
    dJoule_dE = dJH_E[lat, lon, :-1] ** (1/2) / Joule_Heating[lat, lon, :-1]
    dJoule_dVi = dJH_Vi[lat, lon, :-1] ** (1/2) / Joule_Heating[lat, lon, :-1]
    dJoule_dUn = dJH_Un[lat, lon, :-1] ** (1/2) / Joule_Heating[lat, lon, :-1]
    dJoule_dNe = dJH_Ne[lat, lon, :-1] ** (1/2) / Joule_Heating[lat, lon, :-1]

    fig3 = go.Figure()

    fig3.add_trace(go.Scatter(x=Joule_rel, y=heights[lat, lon, :-1], name="Joule Heating error", mode='lines',
                              line=dict(shape='spline', color='green')))
    fig3.add_trace(go.Scatter(x=dJoule_dB, y=heights[lat, lon, :-1], name="dB", mode='lines',
                              line=dict(shape='spline', dash="dot", color='green')))
    fig3.add_trace(go.Scatter(x=dJoule_dE, y=heights[lat, lon, :-1], name="dE", mode='lines',
                              line=dict(shape='spline', dash="dash", color='green')))
    fig3.add_trace(go.Scatter(x=dJoule_dVi, y=heights[lat, lon, :-1], name="dVi", mode='lines',
                              line=dict(shape='spline', dash="dot", color='lime')))
    fig3.add_trace(go.Scatter(x=dJoule_dUn, y=heights[lat, lon, :-1], name="dUn", mode='lines',
                              line=dict(shape='spline', dash="dash", color='lime')))
    fig3.add_trace(go.Scatter(x=dJoule_dNe, y=heights[lat, lon, :-1], name="dNe", mode='lines',
                              line=dict(shape='spline', dash="dot", color='springgreen')))

    # updating the layout of the figure
    fig3.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', xaxis=dict(range=[0, 0.4]),
                       yaxis=dict(range=[min_alt, max_alt], tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)), xaxis_title="",
                       yaxis_title="$Altitude \ (km)$", width=800, height=650,
                       title={'text': 'Joule Heating Error Contributions' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig3.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig3.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig3.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig3.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig3.show()


def plot_conductivities_error(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=pedersen_con_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="σPedersen error", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=hall_con_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="σHall error", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=parallel_con_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="σParallel error",
                             mode='lines', line=dict(shape='spline', color='green'), visible="legendonly"))
    # updating the layout of the figure
    fig.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                      xaxis_title="$(S/m)$", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Conductivities Absolute Error' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

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

    pedersen_rel = pedersen_con_error[lat, lon, :-1] / pedersen_con[lat, lon, :-1]
    hall_rel = hall_con_error[lat, lon, :-1] / hall_con[lat, lon, :-1]
    parallel_rel = parallel_con_error[lat, lon, :-1] / parallel_con[lat, lon, :-1]

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


def plot_conductivities_contr(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    # Pedersen Conductivity
    pedersen_rel = pedersen_con_error[lat, lon, :-1] / pedersen_con[lat, lon, :-1]
    dped_dB = dsp_B[lat, lon, :-1] ** (1/2) / pedersen_con[lat, lon, :-1]
    dped_dTe = dsp_Te[lat, lon, :-1] ** (1/2) / pedersen_con[lat, lon, :-1]
    dped_dTi = dsp_Ti[lat, lon, :-1] ** (1/2) / pedersen_con[lat, lon, :-1]
    dped_dTn = dsp_Tn[lat, lon, :-1] ** (1/2) / pedersen_con[lat, lon, :-1]
    dped_dNion = dsp_Nion[lat, lon, :-1] ** (1/2) / pedersen_con[lat, lon, :-1]
    dped_dNneutral = dsp_Nneutral[lat, lon, :-1] ** (1/2) / pedersen_con[lat, lon, :-1]
    dped_dNe = dsp_Ne[lat, lon, :-1] ** (1/2) / pedersen_con[lat, lon, :-1]

    fig1 = go.Figure()

    # adding the various plots
    fig1.add_trace(go.Scatter(x=pedersen_rel, y=heights[lat, lon, :-1], name="σPedersen error", mode='lines',
                              line=dict(shape='spline', color='red')))
    fig1.add_trace(go.Scatter(x=dped_dB, y=heights[lat, lon, :-1], name="dB", mode='lines',
                              line=dict(shape='spline', dash='dot', color='red')))
    fig1.add_trace(go.Scatter(x=dped_dTe, y=heights[lat, lon, :-1], name="dTe", mode='lines',
                              line=dict(shape='spline', dash='dot', color='gold')))
    fig1.add_trace(go.Scatter(x=dped_dTn, y=heights[lat, lon, :-1], name="dTn", mode='lines',
                              line=dict(shape='spline', dash='dash', color='coral')))
    fig1.add_trace(go.Scatter(x=dped_dTi, y=heights[lat, lon, :-1], name="dTi", mode='lines',
                              line=dict(shape='spline', dash='dot', color='sienna')))
    fig1.add_trace(go.Scatter(x=dped_dNion, y=heights[lat, lon, :-1], name="dNion", mode='lines',
                              line=dict(shape='spline', dash='dash', color='brown')))
    fig1.add_trace(go.Scatter(x=dped_dNneutral, y=heights[lat, lon, :-1], name="dNn", mode='lines',
                              line=dict(shape='spline', dash='dot', color='tan')))
    fig1.add_trace(go.Scatter(x=dped_dNe, y=heights[lat, lon, :-1], name="dNe", mode='lines',
                              line=dict(shape='spline', dash='dash', color='peru')))

    # updating the layout of the figure
    fig1.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                       tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                       xaxis_title="", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                       title={'text': 'Pedersen Conductivity Error Contributions' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig1.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig1.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig1.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig1.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig1.show()

    # Hall Conductivity
    hall_rel = hall_con_error[lat, lon, :-1] / hall_con[lat, lon, :-1]
    dhall_dB = dsh_B[lat, lon, :-1] ** (1/2) / hall_con[lat, lon, :-1]
    dhall_dTe = dsh_Te[lat, lon, :-1] ** (1/2) / hall_con[lat, lon, :-1]
    dhall_dTi = dsh_Ti[lat, lon, :-1] ** (1/2) / hall_con[lat, lon, :-1]
    dhall_dTn = dsh_Tn[lat, lon, :-1] ** (1/2) / hall_con[lat, lon, :-1]
    dhall_Nion = dsh_Nion[lat, lon, :-1] ** (1/2) / hall_con[lat, lon, :-1]
    dhall_dNe = dsh_Ne[lat, lon, :-1] ** (1/2) / hall_con[lat, lon, :-1]
    dhall_dNneutral = dsh_Nneutral[lat, lon, :-1] ** (1/2) / hall_con[lat, lon, :-1]

    fig2 = go.Figure()

    fig2.add_trace(go.Scatter(x=hall_rel, y=heights[lat, lon, :-1], name="σHall error", mode='lines', line=dict(shape='spline', color='blue')))
    fig2.add_trace(go.Scatter(x=dhall_dB, y=heights[lat, lon, :-1], name="dB", mode='lines',
                              line=dict(shape='spline', dash="dot", color='blue')))
    fig2.add_trace(go.Scatter(x=dhall_dTe, y=heights[lat, lon, :-1], name="dTe", mode='lines',
                              line=dict(shape='spline', dash="dash", color='blue')))
    fig2.add_trace(go.Scatter(x=dhall_dTi, y=heights[lat, lon, :-1], name="dTi", mode='lines',
                              line=dict(shape='spline', dash="dot", color='dodgerblue')))
    fig2.add_trace(go.Scatter(x=dhall_dTn, y=heights[lat, lon, :-1], name="dTn", mode='lines',
                              line=dict(shape='spline', dash="dash", color='dodgerblue')))
    fig2.add_trace(go.Scatter(x=dhall_Nion, y=heights[lat, lon, :-1], name="dNion", mode='lines',
                              line=dict(shape='spline', dash="dot", color='deepskyblue')))
    fig2.add_trace(go.Scatter(x=dhall_dNe, y=heights[lat, lon, :-1], name="dNe", mode='lines',
                              line=dict(shape='spline', dash="dash", color='deepskyblue')))
    fig2.add_trace(go.Scatter(x=dhall_dNneutral, y=heights[lat, lon, :-1], name="dNn", mode='lines',
                              line=dict(shape='spline', dash="dash", color='indigo')))

    # updating the layout of the figure
    fig2.update_layout(xaxis_type="log", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                       tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                       xaxis_title="", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                       title={'text': 'Hall Conductivity Error Contributions' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig2.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig2.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig2.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig2.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig2.show()


def plot_currents_error(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=J_pedersen_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="Pedersen Current error", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=J_hall_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="Hall Current error", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=J_ohmic_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="Ohmic current error", mode='lines',
                             line=dict(shape='spline', color='green')))
    fig.add_trace(go.Scatter(x=J_dens_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="Densities current error", mode='lines',
                             line=dict(shape='spline', color='black')))

    x_range = max(J_dens_error[lat, lon, :-1])

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

    Jp_rel = J_pedersen_error[lat, lon, :-1] / J_pedersen[lat, lon, :-1]
    Jh_rel = J_hall_error[lat, lon, :-1] / J_hall[lat, lon, :-1]
    Johmic_rel = J_ohmic_error[lat, lon, :-1] / J_ohmic[lat, lon, :-1]
    Jdens_rel = J_dens_error[lat, lon, :-1] / J_dens[lat, lon, :-1]

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


def plot_currents_contr(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    # JOhmic
    Johmic_rel = J_ohmic_error[lat, lon, :-1] / J_ohmic[lat, lon, :-1]
    dJohm_dB = dJohm_B[lat, lon, :-1] ** (1/2) / J_ohmic[lat, lon, :-1]
    dJohm_dE = dJohm_E[lat, lon, :-1] ** (1/2) / J_ohmic[lat, lon, :-1]
    dJohm_dUn = dJohm_Un[lat, lon, :-1] ** (1/2) / J_ohmic[lat, lon, :-1]
    dJohm_dsp = dJohm_sp[lat, lon, :-1] ** (1/2) / J_ohmic[lat, lon, :-1]
    dJohm_dsh = dJohm_sh[lat, lon, :-1] ** (1/2) / J_ohmic[lat, lon, :-1]
    dJohm_dTi = dJohm_Ti[lat, lon, :-1] ** (1/2) / J_ohmic[lat, lon, :-1]
    dJohm_dTn = dJohm_Tn[lat, lon, :-1] ** (1/2) / J_ohmic[lat, lon, :-1]
    dJohm_dTe = dJohm_Te[lat, lon, :-1] ** (1/2) / J_ohmic[lat, lon, :-1]
    dJohm_dNe = dJohm_Ne[lat, lon, :-1] ** (1/2) / J_ohmic[lat, lon, :-1]
    dJohm_dNn = dJohm_Nneutral[lat, lon, :-1] ** (1/2) / J_ohmic[lat, lon, :-1]
    dJohm_dNion = dJohm_Nion[lat, lon, :-1] ** (1/2) / J_ohmic[lat, lon, :-1]

    fig1 = go.Figure()

    # adding the various plots
    fig1.add_trace(go.Scatter(x=Johmic_rel, y=heights[lat, lon, :-1], name="Ohmic current error", mode='lines',
                              line=dict(shape='spline', color='green')))
    fig1.add_trace(go.Scatter(x=dJohm_dB, y=heights[lat, lon, :-1], name="dB", mode='lines',
                              line=dict(shape='spline', dash="dot", color='green')))
    fig1.add_trace(go.Scatter(x=dJohm_dE, y=heights[lat, lon, :-1], name="dE", mode='lines',
                              line=dict(shape='spline', dash="dot", color='coral')))
    fig1.add_trace(go.Scatter(x=dJohm_dUn, y=heights[lat, lon, :-1], name="dUn", mode='lines',
                              line=dict(shape='spline', dash="dash", color='coral')))
    fig1.add_trace(go.Scatter(x=dJohm_dsp, y=heights[lat, lon, :-1], name="dσP", mode='lines',
                              line=dict(shape='spline', dash="dash", color='sienna')))
    fig1.add_trace(go.Scatter(x=dJohm_dsh, y=heights[lat, lon, :-1], name="dσH", mode='lines',
                              line=dict(shape='spline', dash="dot", color='gold')))
    fig1.add_trace(go.Scatter(x=dJohm_dTi, y=heights[lat, lon, :-1], name="dTi", mode='lines',
                              line=dict(shape='spline', dash="dot", color='orange')))
    fig1.add_trace(go.Scatter(x=dJohm_dTn, y=heights[lat, lon, :-1], name="dTn", mode='lines',
                              line=dict(shape='spline', dash="dash", color='orange')))
    fig1.add_trace(go.Scatter(x=dJohm_dTe, y=heights[lat, lon, :-1], name="dTe", mode='lines',
                              line=dict(shape='spline', dash="dash", color='crimson')))
    fig1.add_trace(go.Scatter(x=dJohm_dNe, y=heights[lat, lon, :-1], name="dNe", mode='lines',
                              line=dict(shape='spline', dash="dot", color='crimson')))
    fig1.add_trace(go.Scatter(x=dJohm_dNn, y=heights[lat, lon, :-1], name="dNn", mode='lines',
                              line=dict(shape='spline', dash="dot", color='mediumblue')))
    fig1.add_trace(go.Scatter(x=dJohm_dNion, y=heights[lat, lon, :-1], name="dNion", mode='lines',
                              line=dict(shape='spline', dash="dash", color='mediumblue')))

    # updating the layout of the figure
    fig1.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                       tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                       xaxis_title="", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                       title={'text': 'Ohmic Current Contributions' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig1.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig1.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig1.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig1.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig1.show()

    # Jdensities
    Jdens_rel = J_dens_error[lat, lon, :-1] / J_dens[lat, lon, :-1]
    Jdens_dB = dJd_B[lat, lon, :-1] ** (1/2) / J_dens[lat, lon, :-1]
    Jdens_dE = dJd_E[lat, lon, :-1] ** (1/2) / J_dens[lat, lon, :-1]
    Jdens_dVi = dJd_Vi[lat, lon, :-1] ** (1/2) / J_dens[lat, lon, :-1]
    Jdens_dUn = dJd_Un[lat, lon, :-1] ** (1/2) / J_dens[lat, lon, :-1]
    Jdens_dNe = dJd_Ne[lat, lon, :-1] ** (1/2) / J_dens[lat, lon, :-1]

    fig2 = go.Figure()

    fig2.add_trace(go.Scatter(x=Jdens_rel, y=heights[lat, lon, :-1], name="Densities current error", mode='lines',
                              line=dict(shape='spline', color='black')))
    fig2.add_trace(go.Scatter(x=Jdens_dB, y=heights[lat, lon, :-1], name="dB", mode='lines',
                              line=dict(shape='spline', dash="dot", color='gold')))
    fig2.add_trace(go.Scatter(x=Jdens_dE, y=heights[lat, lon, :-1], name="dE", mode='lines',
                              line=dict(shape='spline', dash="dash", color='gold')))
    fig2.add_trace(go.Scatter(x=Jdens_dVi, y=heights[lat, lon, :-1], name="dVi", mode='lines',
                              line=dict(shape='spline', dash="dot", color='red')))
    fig2.add_trace(go.Scatter(x=Jdens_dUn, y=heights[lat, lon, :-1], name="dUn", mode='lines',
                              line=dict(shape='spline', dash="dash", color='red')))
    fig2.add_trace(go.Scatter(x=Jdens_dNe, y=heights[lat, lon, :-1], name="dNe", mode='lines',
                              line=dict(shape='spline', dash="dash", color='coral')))

    # updating the layout of the figure
    fig2.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', xaxis=dict(range=[0, 0.4]),
                       yaxis=dict(range=[min_alt, max_alt],
                       tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                       xaxis_title="", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                       title={'text': 'Densities Current Contributions' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig2.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig2.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig2.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig2.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig2.show()


def plot_csections_error(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=C_Op_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="$O^{+}$", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=C_O2p_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="$O_{2}^{+}$", mode='lines',
                             line=dict(shape='spline', color='blue')))
    fig.add_trace(go.Scatter(x=C_NOp_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="$NO^{+}$", mode='lines',
                             line=dict(shape='spline', color='green')))
    fig.add_trace(go.Scatter(x=C_ion_error[lat, lon, :-1], y=heights[lat, lon, :-1], name="$Avg$", mode='lines',
                             line=dict(shape='spline', color='black')))

    # updating the layout of the figure
    fig.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)), xaxis=dict(range=[0, max(C_Op_error[lat, lon, :-1])]),
                      xaxis_title="$(m^{2})$", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Cross Sections Absolute Error' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

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

    COp_rel = C_Op_error[lat, lon, :-1] / C_Op[lat, lon, :-1]
    CO2p_rel = C_O2p_error[lat, lon, :-1] / C_O2p[lat, lon, :-1]
    CNOp_rel = C_NOp_error[lat, lon, :-1] / C_NOp[lat, lon, :-1]
    Cion_rel = C_ion_error[lat, lon, :-1] / C_ion[lat, lon, :-1]

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


def plot_csections_contr(lat_value, lon_value, min_alt, max_alt):
    print("Plotting.....")
    Sg.popup("_" * 50, "Plotting.....", "_" * 50, title=title, auto_close=True, keep_on_top=True)

    lat = lat_value
    lon = lon_value
    min_alt = min_alt
    max_alt = max_alt

    Cion_rel = C_ion_error[lat, lon, :-1] / C_ion[lat, lon, :-1]
    Cion_dTi = dCion_Ti[lat, lon, :-1] ** (1/2) / C_ion[lat, lon, :-1]
    Cion_dTn = dCion_Tn[lat, lon, :-1] ** (1/2) / C_ion[lat, lon, :-1]
    Cion_dnu = dCion_nu[lat, lon, :-1] ** (1/2) / C_ion[lat, lon, :-1]
    Cion_dNn = dCion_Nneutral[lat, lon, :-1] ** (1/2) / C_ion[lat, lon, :-1]

    fig = go.Figure()

    # adding the various plots
    fig.add_trace(go.Scatter(x=Cion_rel, y=heights[lat, lon, :-1], name="Avg Ion Csection error", mode='lines',
                             line=dict(shape='spline', color='red')))
    fig.add_trace(go.Scatter(x=Cion_dTi, y=heights[lat, lon, :-1], name="dTi", mode='lines',
                             line=dict(shape='spline', dash="dot", color='blue')))
    fig.add_trace(go.Scatter(x=Cion_dTn, y=heights[lat, lon, :-1], name="dTn", mode='lines',
                             line=dict(shape='spline', dash="dash", color='green')))
    fig.add_trace(go.Scatter(x=Cion_dnu, y=heights[lat, lon, :-1], name="dν", mode='lines',
                             line=dict(shape='spline', dash="dot", color='green')))
    fig.add_trace(go.Scatter(x=Cion_dNn, y=heights[lat, lon, :-1], name="dNn", mode='lines',
                             line=dict(shape='spline', dash="dash", color='blue')))

    # updating the layout of the figure
    fig.update_layout(xaxis_type="linear", xaxis_showexponent='all', xaxis_exponentformat='power', yaxis=dict(range=[min_alt, max_alt],
                      tickmode='array', tickvals=np.arange(min_alt, max_alt + 5, 5)),
                      xaxis_title="", yaxis_title="$Altitude \ (km)$", width=800, height=650,
                      title={'text': 'Average Ion Csection Error Contributions' + title, 'y': 0.9, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'})

    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='grey')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig.show()


# ############################### Lat-Lon Map Profile Plots ###############################
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

    sc1 = m1.imshow(Joule_Heating_error[:, :, lev] / Joule_Heating[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

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

    sc2 = m2.imshow(Ohmic_Heating_error[:, :, lev] / Ohmic_Heating[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

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

    sc3 = m3.imshow(Frictional_Heating_error[:, :, lev] / Frictional_Heating[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

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

    sc1 = m1.imshow(((nuOp_error[:, :, lev] + nuO2p_error[:, :, lev] + nuNOp_error[:, :, lev]) / 3) /
                    ((nu_Op_sum[:, :, lev] + nu_O2p_sum[:, :, lev] + nu_NOp_sum[:, :, lev]) / 3),
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

    sc2 = m2.imshow(nue_error[:, :, lev] / nu_e_sum[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

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

    sc1 = m1.imshow(pedersen_con_error[:, :, lev] / pedersen_con[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

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

    sc2 = m2.imshow(hall_con_error[:, :, lev] / hall_con[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

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

    sc3 = m3.imshow(parallel_con_error[:, :, lev] / parallel_con[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

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

    sc1 = m1.imshow(J_ohmic_error[:, :, lev] / J_ohmic[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

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

    sc2 = m2.imshow(J_dens_error[:, :, lev] / J_dens[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

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

    sc1 = m1.imshow(C_ion_error[:, :, lev] / C_ion[:, :, lev], cmap=cm.batlow, interpolation='bicubic')

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


# ####################################### Lat - Alt Profile Plots #######################################
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
    cp1 = plt.contourf(heights_la[:-1], glat_in[:], Joule_Heating_error[:, lon, :-1] / Joule_Heating[:, lon, :-1], locator=ticker.LogLocator(),
                       cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Joule Heating Relative Error' + title)
    cbar = plt.colorbar(cp1)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    # Ohmic Heating
    plt.figure(figsize=(12, 12))
    cp2 = plt.contourf(heights_la[:-1], glat_in[:], Ohmic_Heating_error[:, lon, :-1] / Ohmic_Heating[:, lon, :-1], locator=ticker.LogLocator(),
                       cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Ohmic Heating Relative Error' + title)
    cbar = plt.colorbar(cp2)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    # Frictional Heating
    plt.figure(figsize=(12, 12))
    cp3 = plt.contourf(heights_la[:-1], glat_in[:], Frictional_Heating_error[:, lon, :-1] / Frictional_Heating[:, lon, :-1],
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
    cp1 = plt.contourf(heights_la[:-1], glat_in[:], ((nuOp_error[:, lon, :-1] + nuO2p_error[:, lon, :-1] + nuNOp_error[:, lon, :-1]) / 3) /
                                                    ((nu_Op_sum[:, lon, :-1] + nu_O2p_sum[:, lon, :-1] + nu_NOp_sum[:, lon, :-1]) / 3),
                       cmap=cm.batlow, interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Average Ion Collision Frequency Relative Error' + title)
    cbar = plt.colorbar(cp1)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    # Electron Collision Frequency
    plt.figure(figsize=(12, 12))
    cp2 = plt.contourf(heights_la[:-1], glat_in[:], nue_error[:, lon, :-1] / nu_e_sum[:, lon, :-1], cmap=cm.batlow, interpolation='bicubic')
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
    cp1 = plt.contourf(heights_la[:-1], glat_in[:], pedersen_con_error[:, lon, :-1] / pedersen_con[:, lon, :-1], cmap=cm.batlow,
                       interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Pedersen Conductivity Relative Error' + title)
    cbar = plt.colorbar(cp1)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    # Hall Conductivity
    plt.figure(figsize=(12, 12))
    cp2 = plt.contourf(heights_la[:-1], glat_in[:], hall_con_error[:, lon, :-1] / hall_con[:, lon, :-1], locator=ticker.LogLocator(), cmap=cm.batlow,
                       interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Hall Conductivity Relative Error' + title)
    cbar = plt.colorbar(cp2)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    # Parallel Conductivity
    plt.figure(figsize=(12, 12))
    cp3 = plt.contourf(heights_la[:-1], glat_in[:], parallel_con_error[:, lon, :-1] / parallel_con[:, lon, :-1], cmap=cm.batlow,
                       interpolation='bicubic')
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
    cp1 = plt.contourf(heights_la[:-1], glat_in[:], J_ohmic_error[:, lon, :-1] / J_ohmic[:, lon, :-1], locator=ticker.LogLocator(), cmap=cm.batlow,
                       interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Ohmic Current Relative Error' + title)
    cbar = plt.colorbar(cp1)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    # Densities Current
    plt.figure(figsize=(12, 12))
    cp2 = plt.contourf(heights_la[:-1], glat_in[:], J_dens_error[:, lon, :-1] / J_dens[:, lon, :-1], locator=ticker.LogLocator(), cmap=cm.batlow,
                       interpolation='bicubic')
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
    cp1 = plt.contourf(heights_la[:-1], glat_in[:], C_ion_error[:, lon, :-1] / C_ion[:, lon, :-1], locator=ticker.LogLocator(), cmap=cm.batlow,
                       interpolation='bicubic')
    plt.xlim(min_alt, max_alt)
    plt.xlabel('~$Altitude \ (km)$')
    plt.ylabel('$Latitude \ (deg)$')

    plt.title('Average Ion Cross Section Relative Error' + title)
    cbar = plt.colorbar(cp1)
    cbar.set_label(label='', size='large', weight='bold', rotation=270, labelpad=30)

    plt.show()
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END OF PLOTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


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

    # #################### Define Layouts ####################
    # layout for densities frame
    densities_layout = [[Sg.Text("O\u207a", pad=((20, 0), (10, 0)), tooltip="Atomic Oxygen Ion"), Sg.Text("NO\u207a", pad=((20, 0), (10, 0)),
                         tooltip="Nitric Oxide Ion"), Sg.Text("O\u2082\u207a", pad=((20, 0), (10, 0)), tooltip="Molecular Oxygen Ion"),
                         Sg.Text("e\u207b", pad=((25, 0), (10, 0)), tooltip="Electron"), Sg.Text("O", pad=((30, 0), (10, 0)),
                         tooltip="Atomic Oxygen Neutral"), Sg.Text("O\u2082", pad=((30, 0), (10, 0)), tooltip="Molecular Oxygen Neutral"),
                         Sg.Text("N\u2082", pad=((30, 0), (10, 0)), tooltip="Molecular Nitrogen Neutral")],
                        [Sg.Spin(size=(3, 3), values=[i for i in range(0, 101, 2)], initial_value=4, key="-O+-"),
                         Sg.Spin(size=(3, 3), values=[i for i in range(0, 101, 2)], initial_value=4, key="-NO+-"),
                         Sg.Spin(size=(3, 3), values=[i for i in range(0, 101, 2)], initial_value=4, key="-O2+-"),
                         Sg.Spin(size=(3, 3), values=[i for i in range(0, 101, 2)], initial_value=4, key="-e-"),
                         Sg.Spin(size=(3, 3), values=[i for i in range(0, 101, 2)], initial_value=4, key="-O-"),
                         Sg.Spin(size=(3, 3), values=[i for i in range(0, 101, 2)], initial_value=4, key="-O2-"),
                         Sg.Spin(size=(3, 3), values=[i for i in range(0, 101, 2)], initial_value=4, key="-N2-")]]

    # layout for temperatures frame
    temp_layout = [[Sg.Text("T\u1d62", pad=((40, 0), (10, 0)), tooltip="Ion Temperature"), Sg.Text("T\u2099", pad=((85, 0), (10, 0)),
                    tooltip="Neutral Temperature"), Sg.Text("T\u2091", pad=((85, 0), (10, 0)), tooltip="Electron Temperature")],
                   [Sg.Spin(size=(3, 3), pad=(30, 0), values=[i for i in range(0, 101, 2)], initial_value=4, key="-Ti-"),
                    Sg.Spin(size=(3, 3), pad=(30, 0), values=[i for i in range(0, 101, 2)], initial_value=4, key="-Tn-"),
                    Sg.Spin(size=(3, 3), pad=(30, 0), values=[i for i in range(0, 101, 2)], initial_value=4, key="-Te-")]]

    # layout for wind - velocity frame
    win_vel_layout = [[Sg.Text("U\u2099", pad=((80, 0), (10, 0)), tooltip="Neutral wind"), Sg.Text("V\u1d62", pad=((110, 0), (10, 0)),
                       tooltip="Ion velocity")], [Sg.Spin(pad=(70, 0), size=(3, 3), values=[i for i in range(0, 101, 2)],
                                                  initial_value=4, key="-Un-"),
                      Sg.Spin(pad=(20, 0), size=(3, 3), values=[i for i in range(0, 101, 2)], initial_value=4, key="-Vi-")]]

    # layout for electric and magnetic field frame
    fields_layout = [[Sg.Text("E", pad=((80, 0), (10, 0)), tooltip="Electric field"), Sg.Text("B", pad=((120, 0), (10, 0)),
                      tooltip="Magnetic field")], [Sg.Spin(pad=(70, 0), size=(3, 3), values=[i for i in range(0, 101, 2)],
                                                   initial_value=4, key="-E-"),
                     Sg.Spin(pad=(20, 0), size=(3, 3), values=[i for i in range(0, 101, 2)], initial_value=4, key="-B-")]]

    # used in plotting choices frame
    col1 = [[Sg.Checkbox("Plot Heating Rates Absolute Error", default=False, tooltip="Plots heating rates absolute error", key="-HR_abs-")],
            [Sg.Checkbox("Plot Heating Rates Relative Error", default=False, tooltip="Plots heating rates relative error", key="-HR_rel-")],
            [Sg.Checkbox("Plot Heating Rates Relative Error Contributions", default=False,
             tooltip="Plots heating rates relative error with contributions", key="-HR_con-")],
            [Sg.Checkbox("Plot Collision Frequencies Absolute Error", default=False, tooltip="Plots collision frequencies absolute error",
             key="-COL_abs-")],
            [Sg.Checkbox("Plot Collision Frequencies Relative Error", default=False, tooltip="Plots collision frequencies relative error",
             key="-COL_rel-")],
            [Sg.Checkbox("Plot Collision Frequencies Relative Error Contributions", default=False,
             tooltip="Plots collision frequencies relative error with contributions", key="-COL_con-")],
            [Sg.Checkbox("Plot Conductivities Absolute Error", default=False, tooltip="Plots conductivities absolute error", key="-CON_abs-")],
            [Sg.Checkbox("Plot Conductivities Relative Error", default=False, tooltip="Plots conductivities relative error", key="-CON_rel-")],
            [Sg.Checkbox("Plot Conductivities Relative Error Contributions", default=False,
             tooltip="Plots conductivities relative error with contributions", key="-CON_con-")],
            [Sg.Checkbox("Plot Currents Absolute Error", default=False, tooltip="Plots currents absolute error", key="-CUR_abs-")],
            [Sg.Checkbox("Plot Currents Relative Error", default=False, tooltip="Plots currents relative error", key="-CUR_rel-")],
            [Sg.Checkbox("Plot Currents Relative Contributions", default=False, tooltip="Plots currents relative error with contributions",
             key="-CUR_con-")],
            [Sg.Checkbox("Plot Cross Sections Absolute Error", default=False, tooltip="Plots cross sections absolute error", key="-CR_abs-")],
            [Sg.Checkbox("Plot Cross Sections Relative Error", default=False, tooltip="Plots cross sections relative error", key="-CR_rel-")],
            [Sg.Checkbox("Plot Cross Sections Relative Error Contributions", default=False,
             tooltip="Plots cross sections relative error with contributions", key="-CR_con-")]]

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
                [Sg.Text("_"*15, pad=((20, 0), (0, 0))), Sg.Text("_"*30, pad=((105, 0), (0, 0)))],
                [Sg.Column(col2, pad=((0, 0), (0, 100)), scrollable=False), Sg.Column(col1, size=(350, 250), pad=((10, 0), (0, 0)), scrollable=True)]]

    templay2 = [[Sg.Text("Products", pad=((30, 0), (30, 0))), Sg.Text("Errors", pad=((200, 0), (30, 0)))],
                [Sg.Text("_"*15, pad=((20, 0), (0, 0))), Sg.Text("_"*30, pad=((105, 0), (0, 0)))],
                [Sg.Column(col3, scrollable=False), Sg.Column(col4, scrollable=False)]]

    templay3 = [[Sg.Text("min altitude(km)", pad=((30, 0), (15, 0))), Sg.Text("max altitude(km)", pad=((30, 10), (15, 0)))],
                [Sg.InputCombo(values=[i for i in range(100, 601, 10)], pad=((40, 0), (10, 20)), size=(10, 1), default_value="110",
                 key="-min_alt_la-"),
                 Sg.InputCombo(values=[i for i in range(100, 601, 10)], pad=((45, 0), (10, 20)), size=(11, 1), default_value="200",
                 key="-max_alt_la-")],
                [Sg.Text("Products", pad=((30, 0), (30, 0))), Sg.Text("Errors", pad=((200, 0), (30, 0)))],
                [Sg.Text("_"*15, pad=((20, 0), (0, 0))), Sg.Text("_"*30, pad=((105, 0), (0, 0)))],
                [Sg.Column(col5, scrollable=False), Sg.Column(col6, scrollable=False)]]

    templay4 = [[Sg.TabGroup([[Sg.Tab("Vertical Profile Plots", templay1), Sg.Tab("Map Profile (Lat-Lon) Plots", templay2),
                                  Sg.Tab("Map Profile (Lat-Alt) Plots", templay3)]], key="-TABGROUP1-")]]

    # frame used in densities
    den_frame = [Sg.Frame("Densities error(%)", densities_layout, pad=((0, 0), (15, 10)))]
    temp_frame = [Sg.Frame("Temperatures error(%)", temp_layout, pad=((0, 0), (15, 10)))]
    wv_frame = [Sg.Frame("Neutral wind - Ion velocity error(%)", win_vel_layout, size=(20, 1), pad=((0, 0), (15, 10)))]
    fields_frame = [Sg.Frame("Electric - Magnetic field error(%)", fields_layout, size=(40, 1), pad=((0, 0), (15, 10)))]

    # layout used in percentage errors frame
    perc_errors_layout = [den_frame, temp_frame, wv_frame, fields_frame]

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
                   [Sg.Frame("Percentage Errors", perc_errors_layout, pad=((0, 0), (10, 4))),
                    Sg.Frame("Choose plots", templay4, pad=((0, 0), (40, 30)))],
                   [Sg.Checkbox("Science study errors", default=True, tooltip="True errors", pad=((0, 800), (0, 10)), key="-ERROR-")],
                   [Sg.Text("Choose Profile")],
                   [Sg.TabGroup([[Sg.Tab("Vertical Profile", vert_layout), Sg.Tab("Map Profile (Lat-Lon)", map_layout),
                                  Sg.Tab("Map Profile (Lat-Alt)", map2_latout)]], key="-TABGROUP-")],
                   [button("Calculate Products", "Click to start calculation"), button("Calculate Error", "Click to start calculation"),
                    button("Help", "Click for info"),
                    button("Exit", "Click to exit program")]]
    # ################################################# MAIN LAYOUT END ################################################
    # ##################################################################################################################
    # ################################################# Create Window ##################################################
    window = Sg.Window("Joule Heating Error Propagation GUI", main_layout, element_justification="c", keep_on_top=False)
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
        b_error = float(values["-B-"])
        e_error = float(values["-E-"])
        no_error = float(values["-O-"])
        no2_error = float(values["-O2-"])
        nn2_error = float(values["-N2-"])
        nop_error = float(values["-O+-"])
        no2p_error = float(values["-O2+-"])
        nnop_error = float(values["-NO+-"])
        ne_error = float(values["-e-"])
        te_error = float(values["-Te-"])
        ti_error = float(values["-Ti-"])
        tn_error = float(values["-Tn-"])
        un_error = float(values["-Un-"])
        vi_error = float(values["-Vi-"])

        ERROR_FLAG = values["-ERROR-"]
        if event == "Calculate Products" and values["-TABGROUP-"] == "Vertical Profile":
            user_lat = values["-LAT-"]
            user_lon = values["-LON-"]
            min_alt = values["-min_alt-"]
            max_alt = values["-max_alt-"]
            user_time = int(values["-TIME_vert-"])
            models_input(file_name=user_file_name, timer=user_time, lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon])
            products(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon])
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
            products(pressure_level=user_lev)
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
            products(lon_value=lon_dictionary[user_lon])
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
                error(error_flag=ERROR_FLAG, B_error=b_error, E_error=e_error, NO_error=no_error, NO2_error=no2_error, NN2_error=nn2_error,
                      NOp_error=nop_error, NO2p_error=no2p_error, NNOp_error=nnop_error, Ne_error=ne_error, Te_error=te_error, Ti_error=ti_error,
                      Tn_error=tn_error, Un_error=un_error, Vi_error=vi_error, lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon])
                if values["-TABGROUP1-"] == "Vertical Profile Plots":
                    min_alt = values["-min_alt-"]
                    max_alt = values["-max_alt-"]
                    if values["-COL_abs-"]:
                        plot_collisions_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                              max_alt=max_alt)
                    if values["-COL_rel-"]:
                        plot_collisions_rel_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                                  max_alt=max_alt)
                    if values["-COL_con-"]:
                        plot_collisions_contr(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                              max_alt=max_alt)
                    if values["-HR_abs-"]:
                        plot_heating_rates_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                                 max_alt=max_alt)
                    if values["-HR_rel-"]:
                        plot_heating_rates_rel_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                                     max_alt=max_alt)
                    if values["-HR_con-"]:
                        plot_heating_rates_contr(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                                 max_alt=max_alt)
                    if values["-CON_abs-"]:
                        plot_conductivities_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                                  max_alt=max_alt)
                    if values["-CON_rel-"]:
                        plot_conductivities_rel_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                                      max_alt=max_alt)
                    if values["-CON_con-"]:
                        plot_conductivities_contr(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                                  max_alt=max_alt)
                    if values["-CUR_abs-"]:
                        plot_currents_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt, max_alt=max_alt)
                    if values["-CUR_rel-"]:
                        plot_currents_rel_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                                max_alt=max_alt)
                    if values["-CUR_con-"]:
                        plot_currents_contr(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt, max_alt=max_alt)
                    if values["-CR_abs-"]:
                        plot_csections_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt, max_alt=max_alt)
                    if values["-CR_rel-"]:
                        plot_csections_rel_error(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt,
                                                 max_alt=max_alt)
                    if values["-CR_con-"]:
                        plot_csections_contr(lat_value=lat_dictionary[user_lat], lon_value=lon_dictionary[user_lon], min_alt=min_alt, max_alt=max_alt)
        if event == "Calculate Error" and values["-TABGROUP-"] == "Map Profile (Lat-Lon)":
            user_lev = values["-Pr_level-"]
            error(error_flag=ERROR_FLAG, B_error=b_error, E_error=e_error, NO_error=no_error, NO2_error=no2_error, NN2_error=nn2_error,
                  NOp_error=nop_error, NO2p_error=no2p_error, NNOp_error=nnop_error, Ne_error=ne_error, Te_error=te_error, Ti_error=ti_error,
                  Tn_error=tn_error, Un_error=un_error, Vi_error=vi_error, pressure_level=user_lev)
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
        if event == "Calculate Error" and values["-TABGROUP-"] == "Map Profile (Lat-Alt)":
            user_lon = values["-Lon_map2-"]
            min_alt_la = values["-min_alt_la-"]
            max_alt_la = values["-max_alt_la-"]
            error(error_flag=ERROR_FLAG, B_error=b_error, E_error=e_error, NO_error=no_error, NO2_error=no2_error, NN2_error=nn2_error,
                  NOp_error=nop_error, NO2p_error=no2p_error, NNOp_error=nnop_error, Ne_error=ne_error, Te_error=te_error, Ti_error=ti_error,
                  Tn_error=tn_error, Un_error=un_error, Vi_error=vi_error, lon_value=lon_dictionary[user_lon])
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
