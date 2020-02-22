"""
Created on Sun Jan  12 18:00:33 2020

@author: Panagiotis Pyrnaris, Christos Psychalas
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

m_Oplus = 2.66 * (10 ** (-26))            # kg
m_O2plus = 2 * m_Oplus                    # kg
m_NOplus = m_Oplus + 2.32 * (10 ** (-26))  # kg

q_e = 1.60217657 * (10 ** (-19))  # coulomb

m_e = 9.11 * (10 ** (-31))  # kg

cubic = 10 ** 6

# inputMSISE = pd.read_csv("lat80_lon180_varHeight_MSISE00.csv")
# inputIRI16 = pd.read_csv("lat80_lon180_varHeight_IRI16.csv")
# inputIGRF = pd.read_csv("lat80_lon180_varHeight_IGRF.csv")

input = pd.read_csv("lat70_lon0_varHeight_ap_7_f107_90_201503170000Output.csv")

#input2 = pd.read_csv("lat70_lon0_12_March17_varHeight_Pedersen_OurCalc1.csv")
#Pedersen_k2 = input2["Pedersen_S/m"]
#Hall_k2 = input2["Hall S/m"]
#sigma_02 = input2["Parallel S/m"]

time = input["Epoch(UTCG)"]
glat = input["Lat_GEOD(deg)"]
glon = input["Lon_GEOD(deg)"]
height = input["Height_WGS84 (km)"]
Te = input["Te_K"]
Ti = input["Ti_K"]
Ne = input["Ne_cm^-3"]  #
Op = input["Op_cm^-3"]  #
O2p = input["O2p_cm^-3"]  #
NOp = input["NOp_cm^-3"]  #
Tn = input["Tn_K"]
O = input["O_cm^-3"]  #
O2 = input["O2_cm^-3"]  #
N2 = input["N2_cm^-3"]  #
B = input["B_Tesla"]

VO2p = []  # sec^-1
VOp = []   # sec^-1
VNOp = []  # sec^-1

Vei = []   # sec^-1
Ven = []   # sec^-1
Ve_total = []  # sec^-1

Omega_e = []    # sec^-1
Omega_Op = []   # sec^-1
Omega_O2p = []  # sec^-1
Omega_NOp = []  # sec^-1

Pedersen_k = []  # S/m
sigma_0 = []    # S/m
Hall_k = []     # S/m

for i in range(0, input.shape[0]):

    Tr = (Ti[i] + Tn[i]) / 2  # in kelvin

    #  ion-neutral collisions
    if Tr > 800:

        VO2p_O2 = 2.59 * (10 ** (-11)) * O2[i] * np.sqrt(Tr) * ((1 - 0.073 * np.log10(Tr)) ** 2)
        VOp_O = 3.67 * O[i] * (10 ** (-11)) * np.sqrt(Tr) * ((1 - 0.064 * np.log10(Tr)) ** 2)

    elif (Tr <= 800) & (Tr > 235):

        VO2p_O2 = 0
        VOp_O = 3.67 * O[i] * (10 ** (-11)) * np.sqrt(Tr) * ((1 - 0.064 * np.log10(Tr)) ** 2)

    else:

        VO2p_O2 = 0
        VOp_O = 0

    VO2p_O = 2.31 * (10 ** (-10)) * O[i]
    VO2p_N2 = 4.13 * (10 ** (-10)) * N2[i]

    VO2p.append(VO2p_O2 + VO2p_O + VO2p_N2)

    VOp_O2 = 6.64 * (10 ** (-10)) * O2[i]
    VOp_N2 = 6.82 * (10 ** (-10)) * N2[i]

    VOp.append(VOp_O + VOp_O2 + VOp_N2)

    VNOp_O2 = 4.27 * (10 ** (-10)) * O2[i]
    VNOp_O = 2.44 * (10 ** (-10)) * O[i]
    VNOp_N2 = 4.39 * (10 ** (-10)) * N2[i]

    VNOp.append(VNOp_O2 + VNOp_O + VNOp_N2)

    # electron-neutral collisions
    Ve_N2 = 2.33 * (10 ** (-11)) * N2[i] * (1 - 1.21 * (10 ** (-4)) * Te[i]) * Te[i]
    Ve_O2 = 1.82 * (10 ** (-10)) * O2[i] * (1 + 3.6 * (10 ** (-2)) * np.sqrt(Te[i])) * np.sqrt(Te[i])
    Ve_O = 8.9 * (10 ** (-11)) * O[i] * (1 + 5.7 * (10 ** (-4)) * Te[i]) * np.sqrt(Te[i])

    Ven.append(Ve_N2 + Ve_O2 + Ve_O)

    # electron-ion collisions
    Ve_Op = (54.5 * Op[i]) / (Te[i] ** (3/2))
    Ve_O2p = (54.5 * O2p[i]) / (Te[i] ** (3/2))
    Ve_NOp = (54.5 * NOp[i]) / (Te[i] ** (3/2))

    Vei.append(Ve_Op + Ve_O2p + Ve_NOp)

    # total electron collisions
    Ve_total.append(Ven[i] + Vei[i])

    # GyroFrequencies
    Omega_e.append(q_e * B[i] / m_e)
    Omega_Op.append(q_e * B[i] / m_Oplus)
    Omega_O2p.append(q_e * B[i] / m_O2plus)
    Omega_NOp.append(q_e * B[i] / m_NOplus)

    # electron conductivity(sigma_e) / (NOT in S/m)
    sigma_e = (Ne[i] * q_e ** 2) / (m_e * Ven[i])

    # ion conductivities / (NOT in S/m)
    sigma_Op = (Op[i] * q_e ** 2) / (m_Oplus * VOp[i])
    sigma_O2p = (O2p[i] * q_e ** 2) / (m_O2plus * VO2p[i])
    sigma_NOp = (NOp[i] * q_e ** 2) / (m_NOplus * VNOp[i])

    # Parallel conductivity(sigma_0) / (conversion to S/m using cubic constant)
    sigma_0.append(((Ne[i] * q_e ** 2) / (m_e * Ve_total[i])) * cubic)

    # Hall conductivity(sigma_2) /  (conversion to S/m using cubic constant)
    Hall_k.append(((sigma_e * Ven[i] * Omega_e[i]) / (Ven[i] ** 2 + Omega_e[i] ** 2)
                - (sigma_Op * VOp[i] * Omega_Op[i]) / (VOp[i] ** 2 + Omega_Op[i] ** 2)
                - (sigma_O2p * VO2p[i] * Omega_O2p[i]) / (VO2p[i] ** 2 + Omega_O2p[i] ** 2)
                - (sigma_NOp * VNOp[i] * Omega_NOp[i]) / (VNOp[i] ** 2 + Omega_NOp[i] ** 2)) * cubic)

    # Pedersen conductivity(sigma_1) / (conversion to S/m using cubic constant)
    Pedersen_k.append(((sigma_Op * (VOp[i] ** 2)) / (VOp[i] ** 2 + Omega_Op[i] ** 2)
                 + (sigma_O2p * (VO2p[i] ** 2)) / (VO2p[i] ** 2 + Omega_O2p[i] ** 2)
                 + (sigma_NOp * (VNOp[i] ** 2)) / (VNOp[i] ** 2 + Omega_NOp[i] ** 2)
                 + (sigma_e * (Ven[i] ** 2)) / (Ven[i] ** 2 + Omega_e[i] ** 2)) * cubic)


# prepare and export
data = {"Epoch(UTCG)": time,
        "Lat_GEOD(deg)": glat,
        "Lon_GEOD(deg)": glon,
        "Height_WGS84 (km)": height,
        "Pedersen_S/m": Pedersen_k,
        "Hall S/m": Hall_k,
        "Parallel S/m": sigma_0,
        # "Pedersen_r_S/m": Pedersen_r,
        }

Pedersen = pd.DataFrame(data)
Pedersen.to_csv("lat70_lon0_varHeight_ap_7_f107_90_201503170000Output_Pedersen_OurCalc1.csv", index=None)


# plotting
fig = plt.figure()
plt.grid()
plt.yticks(np.arange(0, 400, 25))
plt.xlabel('frequency (Hz)')
plt.ylabel('Altitude (km)')
plt.title("Collisions(March 17 (00:00) glat=70 glon=0 ap_7 )")
plt.xscale('log')

#plt.plot(Pedersen_k, height, color='blue', label='Pedersen')
#plt.plot(Hall_k, height, color='red', label='Hall')
#plt.plot(sigma_0, height, color='green', label='Parallel')

plt.plot(Ven, height, linestyle=':', color='blue', label='e-neutral collisions')
plt.plot(Ve_total, height, linestyle=':', color='purple', label='total e collisions')
plt.plot(VNOp, height, linestyle=':', color='red', label='NOplus collisions')
plt.plot(VO2p, height, linestyle=':', color='green', label='O2plus collisions')
plt.plot(VOp, height, linestyle=':', color='black', label='Oplus collisions')

#plt.plot(Pedersen_k2, height, linestyle=':',  color='blue', label='Pedersen(12:00')
#plt.plot(Hall_k2, height, linestyle=':',  color='red', label='Hall(12:00)')
#plt.plot(sigma_02, height, linestyle=':',  color='green', label='Parallel(12:00)')

plt.legend(loc='best')

plt.show()

# end
