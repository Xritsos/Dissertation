"""
Created on Thu Jan  2 19:35:17 2020

@author: Panagiotis Pyrnaris, Christos Psichalas
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# constants for Collision Frequency Calculation
c1 = 2.33 * 10 ** (-11)
c2 = 1.21 * 10 ** (-4)
c3 = 1.82 * 10 ** (-10)
c4 = 3.6 * 10 ** (-2)
c5 = 8.9 * 10 ** (-11)
c6 = 5.7 * 10 ** (-4)
c7 = 4.27 * 10 ** 10
c8 = 2.44 * 10 ** 10
c9 = 2.44 * 10 ** 10
c10 = 6.82 * 10 ** 10
c11 = 3.67 * 10 ** (-11)
c12 = 6.82 * 10 ** 10
c13 = 2.59 * 10 ** (-11)
c14 = 2.31 * 10 ** 10
c15 = 4.13 * 10 ** 10
c16 = 0.073
c17 = 0.064
fb = 1.7  # Burnside Factor, the factor that connect theoretical and practical estimations of Oplus [from 1.7-0.3 to 1.7+0.7]

ArO = 15.9994  # atomic mass O=16
ArN = 14.0067  # atomic mass N=14
NAvog = 6.02214076e23
q_e = 1.60217657e-19  # Coulomb
m_NOplus = (ArO + ArN) / (NAvog * 1000)
m_O2plus = 2 * ArO / (NAvog * 1000)
m_Oplus = 1 * ArO / (NAvog * 1000)  # kg
m_e = 9.11e-37  # kg
CubicCm2CubicM = 10 ** (6)

# read inputs and convert to m^-3
inputMSISE = pd.read_csv("lat80_lon0_varHeight_MSISE00.csv")
inputIRI16 = pd.read_csv("lat80_lon0_varHeight_IRI16.csv")
inputIGRF = pd.read_csv("lat80_lon0_varHeight_IGRF.csv")

time = inputIRI16["Epoch(UTCG)"]
glat = inputIRI16["Lat_GEOD(deg)"]
glon = inputIRI16["Lon_GEOD(deg)"]
height = inputIRI16["Height_WGS84 (km)"]
Te = inputIRI16["Te_iri16_K"]
Ti = inputIRI16["Ti_iri16_K"]
Ne = inputIRI16["ne_iri16_cm-3"]  #
Op = inputIRI16["O+_iri16_cm-3"]  #
O2p = inputIRI16["O2+_iri16_cm-3"]  #
NOp = inputIRI16["NO+_iri16_cm-3"]  #
Tn = inputMSISE["Tn_msise00_K"]
O = inputMSISE["O_msise00_cm-3"]  #
O2 = inputMSISE["O2_msise00_cm-3"]  #
N2 = inputMSISE["N2_msise00_cm-3"]  #
#H=inputMSISE[]
B = inputIGRF["B_Tesla"]

VO2p = []
VOp = []
VNOp = []
Ve = []
Omega_e = []
Omega_Op = []
Omega_O2p = []
Omega_NOp = []
r_e = []
r_Op = []
r_O2p = []
r_NOp = []
k_e = []
k_Op = []
k_O2p = []
k_NOp = []
Pedersen_k = []
for i in range(0, inputMSISE.shape[0]):
    if (height[i] <= 200):
        # collision frequencies
        # VO2p=VO2p_O2+VO2p_O+VO2p_N2
        G = Ti[i] + Tn[i]
        VO2p_O2 = 0
        if (height[i] <= 147) & (height[i] >= 80):
            VO2p_O2 = 8.2 * 10 ** (-10) * O2[i]
        if (height[i] <= 200) & (height[i] > 147):
            VO2p_O2 = 3.4 * 10 ** (-13) * O2[i] * np.sqrt(G) * (10.6 - 0.76 * np.log10(G)) ** 2
        VO2p_O = 7.51 * 10 ** (-22) * O[i]
        VO2p_N2 = 8.925 * 10 ** (-22) * N2[i]
        VO2p.append(VO2p_O2 + VO2p_O + VO2p_N2)

        # VOp=VOp_O+VOp_O2+VOp_N2
        VOp_O = 0
        if (height[i] <= 106) & (height[i] >= 80):
            VOp_O = 8.6 * 10 ** (-10) * O[i]
        if (height[i] <= 200) & (height[i] > 106):
            VOp_O = 4.7 * 10 ** (-13) * O[i] * np.sqrt(G) * (10.5 - 0.67 * np.log10(G)) ** 2
        VOp_O2 = 1.007 * 10 ** (-21)*O2[i]
        VOp_N2 = 1.081 * 10 ** (-21) * N2[i]
        VOp.append(VOp_O + VOp_O2 + VOp_N2)

        # VNOp=VNOp_O2+VNOp_O+VNOp_N2
        VNOp_O2 = 8.357 * 10 ** (-22) * O2[i]
        VNOp_O = 7.593 * 10 ** (-22) * O[i]
        VNOp_N2 = 9.0623 * 10 ** (-22) * N2[i]
        VNOp.append(VNOp_O2 + VNOp_O + VNOp_N2)

        # Ve=Ve_N2+Ve_O2+Ve_O
        Ve_N2 = 2.33 * 10 ** (-11) * N2[i] * (1 - 1.21 * 10 ** (-4) * Te[i]) * Te[i]
        Ve_O2 = 1.82 * 10 ** (-10) * O2[i] * (1 + 3.6 * 10 ** (-2) * np.sqrt(Te[i])) * np.sqrt(Te[i])
        Ve_O = 8.2 * 10 ** (-10) * O[i] * np.sqrt(Te[i])
       # Ve_H = 4.5 * 10 **(-9) * H[i] * np.sqrt(Te[i]) * (1 - 1.35 * 10 ** (-4) * Te)
        Ve.append(Ve_N2 + Ve_O2 + Ve_O)

        # GyroFrequencies
        Omega_e.append(q_e * B[i] / m_e)
        Omega_Op.append(q_e * B[i] / m_Oplus)
        Omega_O2p.append(q_e * B[i] / m_O2plus)
        Omega_NOp.append(q_e * B[i] / m_NOplus)

        # ratios
        # apo tis aplopoiiseis twn klasmatwn eite doulepsoume me to r eite me to K prokuptei to idio apotelesma
        r_e.append(Ve[i] / Omega_e[i])
        r_Op.append(VOp[i] / Omega_Op[i])
        r_O2p.append(VO2p[i] / Omega_O2p[i])
        r_NOp.append(VNOp[i] / Omega_NOp[i])

        k_e.append(1 / r_e[i])
        k_Op.append(1 / r_Op[i])
        k_O2p.append(1 / r_O2p[i])
        k_NOp.append(1 / r_NOp[i])

        # pedersen Conductivity
        Pedersen_k.append(q_e / B[i] * \
                        (Op[i] * (k_Op[i] / (1 + k_Op[i] ** 2)) + \
                         O2p[i] * (k_O2p[i] / (1 + k_O2p[i] ** 2)) + \
                         NOp[i] * (k_NOp[i] / (1 + k_NOp[i] ** 2)) + \
                         Ne[i] * (k_e[i] / (1 + k_e[i] ** 2))))
    else:
        VO2p.append(None)
        VOp.append(None)
        VNOp.append(None)
        Ve.append(None)
        Omega_e.append(None)
        Omega_Op.append(None)
        Omega_O2p.append(None)
        Omega_NOp.append(None)
        r_e.append(None)
        r_Op.append(None)
        r_O2p.append(None)
        r_NOp.append(None)
        k_e.append(None)
        k_Op.append(None)
        k_O2p.append(None)
        k_NOp.append(None)
        Pedersen_k.append(None)



# prepare and export
data = {"Epoch(UTCG)": time,
        "Lat_GEOD(deg)": glat,
        "Lon_GEOD(deg)": glon,
        "Height_WGS84 (km)": height,
        "Pedersen_S/m": Pedersen_k,
        # "Pedersen_r_S/m": Pedersen_r,
        }

Pedersen = pd.DataFrame(data)
Pedersen.to_csv("lat80_lon0_varHeight_Pedersen_OurCalc.csv", index=None)

# plotting
fig = plt.figure()
plt.grid()
plt.yticks(np.arange(0, 400, 25))
plt.xlabel('Pedersen Conductivity (S/m)')
plt.ylabel('Altitude (km)')
plt.title("Pedersen Conductivity")
plt.xscale('log')

plt.plot(Pedersen_k, height, color='blue')

plt.show()

# end
