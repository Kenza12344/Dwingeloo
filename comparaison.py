import matplotlib.pyplot as plt
import numpy as np
import os

def load_antenna_data(file):
    npzfile = np.load(file, allow_pickle=False)
    data = {
        'theta': npzfile['theta'],
        'e_plane': npzfile['e_plane'],
        'h_plane': npzfile['h_plane'],
        'theta_step': npzfile['theta_step'],
        'gain': npzfile['gain']
    }
    return data

def calculate_parameters(data):
    max_e = np.max(np.abs(data['e_plane']))
    max_h = np.max(np.abs(data['h_plane']))
    max_c = np.max(np.abs(data['e_plane']) + np.abs(data['h_plane'])) / np.sqrt(2)
    theta_rad = np.pi * data['theta'] / 180
    return max_e, max_h, max_c, theta_rad

def cotan(x):
    return -np.tan(x + np.pi / 2)

def ill_eff(field, theta, theta_step_rad, dt_angle):
    theta_edge_rad = np.pi * dt_angle / 180
    full_sphere = np.pi

    u1 = np.abs(np.sum(field[theta <= theta_edge_rad] * np.tan(theta[theta <= theta_edge_rad] / 2) * theta_step_rad)) ** 2
    u2 = np.sum(np.abs(field[theta <= full_sphere]) ** 2 * np.sin(theta[theta <= full_sphere]) * theta_step_rad)
    ill = 2 * (cotan(theta_edge_rad / 2) ** 2) * (u1 / u2)
    return ill

def plot_s11_comparison(freqs1, s11_mag1, freqs2, s11_mag2):
    fig, ax = plt.subplots()
    ax.plot(freqs1 / 1e6, 20 * np.log10(s11_mag1), "-", color="darkblue", label="Antenna 1")
    ax.plot(freqs2 / 1e6, 20 * np.log10(s11_mag2), "-", color="darkgreen", label="Antenna 2")
    ax.set_xlabel("freq. (MHz)")
    ax.set_ylabel("S11 (dB)")
    ax.axvline(1300, color='grey', linestyle='--')
    ax.axvline(1420, color='grey', linestyle='--')
    ax.axvline(1612, color='grey', linestyle='--')
    ax.legend()
    ax.set_title("S11 Comparison")
    plt.show()

def plot_vswr_comparison(freqs1, vswr1, freqs2, vswr2):
    fig, ax = plt.subplots()
    ax.plot(freqs1 / 1e6, vswr1, "-", color="darkblue", label="Antenna 1")
    ax.plot(freqs2 / 1e6, vswr2, "-", color="darkgreen", label="Antenna 2")
    ax.set_xlabel("freq. (MHz)")
    ax.set_ylabel("VSWR")
    ax.axvline(1300, color='grey', linestyle='--')
    ax.axvline(1420, color='grey', linestyle='--')
    ax.axvline(1612, color='grey', linestyle='--')
    ax.legend()
    ax.set_title("VSWR Comparison")
    plt.show()

def plot_gain_comparison(freqs1, gains1, freqs2, gains2):
    fig, ax = plt.subplots()
    ax.plot(freqs1 / 1e6, gains1, ".-", color="darkblue", label="Antenna 1")
    ax.plot(freqs2 / 1e6, gains2, ".-", color="darkgreen", label="Antenna 2")
    ax.set_xlabel("freq. (MHz)")
    ax.set_ylabel("gain (dB)")
    ax.axvline(1300, color='grey', linestyle='--')
    ax.axvline(1420, color='grey', linestyle='--')
    ax.axvline(1612, color='grey', linestyle='--')
    ax.legend()
    ax.set_title("Gain Comparison")
    plt.show()

def plot_beamwidth_comparison(freqs1, beamwidth1, freqs2, beamwidth2, dt_angle1, dt_angle2):
    fig, ax = plt.subplots()
    ax.plot(freqs1 / 1e6, beamwidth1[:, 0], ".-", color="darkblue", label="Antenna 1 -3 dB")
    ax.plot(freqs1 / 1e6, beamwidth1[:, 1], ".-", color="lightblue", label="Antenna 1 -10 dB")
    ax.plot(freqs1 / 1e6, beamwidth1[:, 2], ".-", color="purple", label="Antenna 1 -20 dB")
    ax.plot(freqs2 / 1e6, beamwidth2[:, 0], ".--", color="darkgreen", label="Antenna 2 -3 dB")
    ax.plot(freqs2 / 1e6, beamwidth2[:, 1], ".--", color="lightgreen", label="Antenna 2 -10 dB")
    ax.plot(freqs2 / 1e6, beamwidth2[:, 2], ".--", color="pink", label="Antenna 2 -20 dB")
    ax.axhline(dt_angle1, color="black", linestyle='--')
    ax.axhline(dt_angle2, color="grey", linestyle='--')
    ax.axvline(1300, color='grey', linestyle='--')
    ax.axvline(1420, color='grey', linestyle='--')
    ax.axvline(1612, color='grey', linestyle='--')
    ax.set_ylim(0, 90)
    ax.set_xlabel("freq. (MHz)")
    ax.set_ylabel("Beamwidth (deg.)")
    ax.set_title("Beamwidth Comparison")
    ax.legend()
    plt.show()

def plot_illumination_comparison(freqs1, nt_e1, nt_h1, nt_r1, nt_l1, freqs2, nt_e2, nt_h2, nt_r2, nt_l2):
    fig, ax = plt.subplots()
    ax.plot(freqs1 / 1e6, nt_e1, ".-", color="darkblue", label="Antenna 1 E-plane")
    ax.plot(freqs1 / 1e6, nt_h1, ".-", color="lightblue", label="Antenna 1 H-plane")
    ax.plot(freqs1 / 1e6, nt_r1, ".-", color="green", label="Antenna 1 RHCP-plane")
    ax.plot(freqs1 / 1e6, nt_l1, ".-", color="lightgreen", label="Antenna 1 LHCP-plane")
    ax.plot(freqs2 / 1e6, nt_e2, ".--", color="darkgreen", label="Antenna 2 E-plane")
    ax.plot(freqs2 / 1e6, nt_h2, ".--", color="lightgreen", label="Antenna 2 H-plane")
    ax.plot(freqs2 / 1e6, nt_r2, ".--", color="purple", label="Antenna 2 RHCP-plane")
    ax.plot(freqs2 / 1e6, nt_l2, ".--", color="pink", label="Antenna 2 LHCP-plane")
    ax.axvline(1300, color='grey', linestyle='--')
    ax.axvline(1420, color='grey', linestyle='--')
    ax.axvline(1612, color='grey', linestyle='--')
    ax.set_ylim(0, 1)
    ax.set_xlabel("freq. (MHz)")
    ax.set_ylabel("illumination efficiency")
    ax.set_title("Illumination Efficiency Comparison")
    ax.legend()
    plt.show()

def calculate_vswr(s11_mag):
    return (1 + s11_mag) / (1 - s11_mag)

# Load data for the first antenna
sim_freq1 = 1300  # in MHz
sim_horn1 = "v2"  # v1 = corrugated horn, v2 = DT feed
file1 = "C:\\Users\\BUNICE\\Downloads\\sonde\\sonde-BandL-" + sim_horn1 + "_" + str(int(sim_freq1*1000)) + "k.npz"
data1 = load_antenna_data(file1)
max_e1, max_h1, max_c1, theta_rad1 = calculate_parameters(data1)

# Load data for the second antenna
sim_freq2 = 1420  # in MHz
file2 = f"C:\\Users\\BUNICE\\Desktop\\stage\\HORN-STOCKERT2\\hornStockert2_{int(sim_freq2*1000)}k.npz"
data2 = load_antenna_data(file2)
max_e2, max_h2, max_c2, theta_rad2 = calculate_parameters(data2)

# Calculate dt_angle for both antennas
f_D1 = 0.48
dt_angle1 = 180 * np.arctan2(f_D1 / 2, f_D1 ** 2 - 1 / 16) / np.pi
f_D2 = 0.3
dt_angle2 = 180 * np.arctan2(f_D2 / 2, f_D2 ** 2 - 1 / 16) / np.pi

# Frequency range for both antennas
f_range1 = np.arange(1100, 1810, 10)
f_range2 = np.arange(1200, 1660, 10)

# Initialize arrays to store the results for both antennas
freqs1 = np.zeros(len(f_range1))
s11_mag1 = np.zeros(len(f_range1))
s11_complex1 = np.zeros(len(f_range1)).astype(complex)
gains1 = np.zeros(len(f_range1))
beamwidth1 = np.zeros((len(f_range1), 3))
nt_e1 = np.zeros(len(f_range1))
nt_h1 = np.zeros(len(f_range1))
nt_r1 = np.zeros(len(f_range1))
nt_l1 = np.zeros(len(f_range1))

freqs2 = np.zeros(len(f_range2))
s11_mag2 = np.zeros(len(f_range2))
s11_complex2 = np.zeros(len(f_range2)).astype(complex)
gains2 = np.zeros(len(f_range2))
beamwidth2 = np.zeros((len(f_range2), 3))
nt_e2 = np.zeros(len(f_range2))
nt_h2 = np.zeros(len(f_range2))
nt_r2 = np.zeros(len(f_range2))
nt_l2 = np.zeros(len(f_range2))

# Process the frequency range for the first antenna
for i in range(len(f_range1)):
    f = f_range1[i]
    file = "C:\\Users\\BUNICE\\Downloads\\sonde\\sonde-BandL-" + sim_horn1 + "_" + str(int(f * 1000)) + "k.npz"
    npzfile = np.load(file, allow_pickle=False)
    freqs1[i] = npzfile['freq']
    s11_mag1[i] = npzfile['s11_mag']
    s11_complex1[i] = npzfile['s11']
    gains1[i] = npzfile['gain']
    beamwidth1[i, :] = npzfile['beamwidth']
    max_c = (np.max(np.abs(npzfile['e_plane']) + np.abs(npzfile['h_plane'])))
    r_plane = np.abs(npzfile['e_plane'] + 1j * npzfile['h_plane']) / (np.sqrt(2) * max_c)
    l_plane = np.abs(npzfile['e_plane'] - 1j * npzfile['h_plane']) / (np.sqrt(2) * max_c)
    nt_e1[i] = ill_eff(npzfile['e_plane'], np.pi * data1['theta'] / 180, np.pi * data1['theta_step'] / 180, dt_angle1)
    nt_h1[i] = ill_eff(npzfile['h_plane'], np.pi * data1['theta'] / 180, np.pi * data1['theta_step'] / 180, dt_angle1)
    nt_r1[i] = ill_eff(r_plane, np.pi * data1['theta'] / 180, np.pi * data1['theta_step'] / 180, dt_angle1)
    nt_l1[i] = ill_eff(l_plane, np.pi * data1['theta'] / 180, np.pi * data1['theta_step'] / 180, dt_angle1)

# Process the frequency range for the second antenna
for i in range(len(f_range2)):
    f = f_range2[i]
    file = f"C:\\Users\\BUNICE\\Desktop\\stage\\HORN-STOCKERT2\\hornStockert2_{int(f*1000)}k.npz"
    npzfile = np.load(file, allow_pickle=False)
    freqs2[i] = npzfile['freq']
    s11_mag2[i] = npzfile['s11_mag']
    s11_complex2[i] = npzfile['s11']
    gains2[i] = npzfile['gain']
    beamwidth2[i, :] = npzfile['beamwidth']
    max_c = (np.max(np.abs(npzfile['e_plane']) + np.abs(npzfile['h_plane'])))
    r_plane = np.abs(npzfile['e_plane'] + 1j * npzfile['h_plane']) / (np.sqrt(2) * max_c)
    l_plane = np.abs(npzfile['e_plane'] - 1j * npzfile['h_plane']) / (np.sqrt(2) * max_c)
    nt_e2[i] = ill_eff(npzfile['e_plane'], np.pi * data2['theta'] / 180, np.pi * data2['theta_step'] / 180, dt_angle2)
    nt_h2[i] = ill_eff(npzfile['h_plane'], np.pi * data2['theta'] / 180, np.pi * data2['theta_step'] / 180, dt_angle2)
    nt_r2[i] = ill_eff(r_plane, np.pi * data2['theta'] / 180, np.pi * data2['theta_step'] / 180, dt_angle2)
    nt_l2[i] = ill_eff(l_plane, np.pi * data2['theta'] / 180, np.pi * data2['theta_step'] / 180, dt_angle2)

# Calculate VSWR for both antennas
vswr1 = calculate_vswr(s11_mag1)
vswr2 = calculate_vswr(s11_mag2)

# Plot comparisons
plot_s11_comparison(freqs1, s11_mag1, freqs2, s11_mag2)
plot_vswr_comparison(freqs1, vswr1, freqs2, vswr2)
plot_gain_comparison(freqs1, gains1, freqs2, gains2)
plot_beamwidth_comparison(freqs1, beamwidth1, freqs2, beamwidth2, dt_angle1, dt_angle2)
plot_illumination_comparison(freqs1, nt_e1, nt_h1, nt_r1, nt_l1, freqs2, nt_e2, nt_h2, nt_r2, nt_l2)
