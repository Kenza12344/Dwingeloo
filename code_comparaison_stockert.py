import matplotlib.pyplot as plt
import numpy as np

# Liste des chemins vers les fichiers de géométrie et leurs étiquettes
geometry_files = [
    ("C:\\Users\\BUNICE\\Desktop\\stage\\HORN-STOCKERT2\\hornStockert2.ge", "Geometry 1"),
    ("C:\\Users\\BUNICE\\Desktop\\stage\\horne-stockert32\\HORNESTOCKERT32.ge", "Geometry 2"),
    ("C:\\Users\\BUNICE\\Desktop\\stage\\horne-stockert-34\\HORNE-STOCKERT34.ge", "Geometry 3"),
    ("C:\\Users\\BUNICE\\Desktop\\stage\\horne-stockert-39\\HORNE-STOCKERT39.ge", "Geometry 4"),
    ("C:\\Users\\BUNICE\\Desktop\\stage\\horne-stockert-36\\HORNE-STOCKERT36.ge", "Geometry 4"),
    ("C:\\Users\\BUNICE\\Desktop\\stage\\hornstockert311\\horne-stockert-311.ge", "Geometry 5"),
    
]

# Fonction pour lire et traiter les coordonnées de géométrie
def read_geometry(file_path):
    guide_coordinates = []
    geometry_coordinates = []

    with open(file_path, 'r') as file:
        lines = file.readlines()

    start_idx = 0
    for i, line in enumerate(lines):
        if "GEOMETRY" in line:
            start_idx = i + 1
            break

    guide_lines_indices = [0, 2]
    for i in guide_lines_indices:
        parts = lines[start_idx + i].split()
        if len(parts) >= 3:
            x = float(parts[1].replace('D', 'E'))
            z = float(parts[2].replace('D', 'E'))
            guide_coordinates.append((z, x))

    for line in lines[start_idx + 3:]:
        parts = line.split()
        if len(parts) >= 3 and parts[0].isdigit() and int(parts[0]) > 1:
            x = float(parts[1].replace('D', 'E'))
            z = float(parts[2].replace('D', 'E'))
            geometry_coordinates.append((z, x))

    return np.array(guide_coordinates), np.array(geometry_coordinates)

# Fonction pour charger les résultats de simulation
def load_simulation(file):
    npzfile = np.load(file, allow_pickle=False)
    return npzfile

def cotan(x):
    return -np.tan(x + np.pi / 2)

def ill_eff(field, theta, theta_step_rad, dt_angle):
    theta_edge_rad = np.pi * dt_angle / 180
    full_sphere = np.pi

    u1 = np.abs(np.sum(field[theta <= theta_edge_rad] * np.tan(theta[theta <= theta_edge_rad] / 2) * theta_step_rad)) ** 2
    u2 = np.sum(np.abs(field[theta <= full_sphere]) ** 2 * np.sin(theta[theta <= full_sphere]) * theta_step_rad)
    ill = 2 * (cotan(theta_edge_rad / 2) ** 2) * (u1 / u2)
    return ill

# Variables pour stocker les données de simulation de différentes géométries
sim_data = []
dt_angle = 180 * np.arctan2(0.3 / 2, 0.3 ** 2 - 1 / 16) / np.pi  # Ajustez f_D selon vos besoins

# Boucle sur chaque fichier de géométrie
for geom_file, label in geometry_files:
    guide_coords, geometry_coords = read_geometry(geom_file)

    # Charger les résultats de simulation pour chaque fréquence
    f_range = np.arange(1200, 1660, 10)
    freqs = np.zeros(len(f_range))
    gains = np.zeros(len(f_range))
    directivity = np.zeros(len(f_range))
    beamwidth = np.zeros((len(f_range), 3))
    s11_mag = np.zeros(len(f_range))
    s11_complex = np.zeros(len(f_range)).astype(complex)
    nt_e = np.zeros(len(f_range))
    nt_h = np.zeros(len(f_range))
    nt_r = np.zeros(len(f_range))
    nt_l = np.zeros(len(f_range))

    for i in range(len(f_range)):
        f = f_range[i]
        simulation_file = f"{geom_file[:-3]}_{int(f*1000)}k.npz"  # Générer le chemin correct pour les fichiers de simulation
        npzfile = load_simulation(simulation_file)
        freqs[i] = npzfile['freq']
        gains[i] = npzfile['gain']
        directivity[i] = npzfile['directivity']
        s11_mag[i] = npzfile['s11_mag']
        s11_complex[i] = npzfile['s11']
        beamwidth[i, :] = npzfile['beamwidth']
        theta = npzfile['theta']
        theta_step = npzfile['theta_step']
        e_plane = npzfile['e_plane']
        h_plane = npzfile['h_plane']
        max_e = np.max(np.abs(e_plane))
        max_h = np.max(np.abs(h_plane))
        max_c = (np.max(np.abs(e_plane) + np.abs(h_plane)))
        r_plane = np.abs(e_plane + 1j * h_plane) / (np.sqrt(2) * max_c)
        l_plane = np.abs(e_plane - 1j * h_plane) / (np.sqrt(2) * max_c)
        nt_e[i] = ill_eff(e_plane, np.pi * theta / 180, np.pi * theta_step / 180, dt_angle)
        nt_h[i] = ill_eff(h_plane, np.pi * theta / 180, np.pi * theta_step / 180, dt_angle)
        nt_r[i] = ill_eff(r_plane, np.pi * theta / 180, np.pi * theta_step / 180, dt_angle)
        nt_l[i] = ill_eff(l_plane, np.pi * theta / 180, np.pi * theta_step / 180, dt_angle)

    sim_data.append({
        "label": label,
        "freqs": freqs,
        "gains": gains,
        "directivity": directivity,
        "beamwidth": beamwidth,
        "s11_mag": s11_mag,
        "nt_e": nt_e,
        "nt_h": nt_h,
        "theta": theta,
        "e_plane": e_plane,
        "h_plane": h_plane,
        "max_e": max_e,
        "max_h": max_h,
    })

# Tracer les résultats pour chaque géométrie

# S11
plt.figure()
for data in sim_data:
    plt.plot(data["freqs"] / 1e6, 20 * np.log10(data["s11_mag"]), label=f"{data['label']}")
plt.xlabel("Frequency (MHz)")
plt.ylabel("S11 (dB)")
plt.axvline(1300)
plt.axvline(1420)
plt.axvline(1612)
plt.title("S11 for Different Geometries")
plt.legend()
plt.show()

# Gains and Directivity
plt.figure()
for data in sim_data:
    plt.plot(data["freqs"] / 1e6, data["gains"], label=f"{data['label']} Gain")
    plt.plot(data["freqs"] / 1e6, data["directivity"], label=f"{data['label']} Directivity")
plt.axvline(1300)
plt.axvline(1420)
plt.axvline(1612)
plt.xlabel("Frequency (MHz)")
plt.ylabel("Gain (dB)")
plt.title("Gain and Directivity for Different Geometries")
plt.legend()
plt.show()

# Beamwidth
plt.figure()
for data in sim_data:
    plt.plot(data["freqs"] / 1e6, data["beamwidth"][:, 0], label=f"{data['label']} -3 dB")
    plt.plot(data["freqs"] / 1e6, data["beamwidth"][:, 1], label=f"{data['label']} -10 dB")
    plt.plot(data["freqs"] / 1e6, data["beamwidth"][:, 2], label=f"{data['label']} -20 dB")
plt.axhline(dt_angle, color="black")
plt.axvline(1300)
plt.axvline(1420)
plt.axvline(1612)
plt.ylim(0, 100)
plt.xlabel("Frequency (MHz)")
plt.ylabel("Beamwidth (deg.)")
plt.title("Beamwidth for Different Geometries")
plt.legend()
plt.show()

# Illumination Efficiency
plt.figure()
for data in sim_data:
    plt.plot(data["freqs"] / 1e6, data["nt_e"], label=f"{data['label']} E-plane")
    plt.plot(data["freqs"] / 1e6, data["nt_h"], label=f"{data['label']} H-plane")
    plt.plot(data["freqs"] / 1e6, (data["nt_h"] + data["nt_e"]) / 2, label=f"{data['label']} E/H-plane")
plt.axvline(1300)
plt.axvline(1420)
plt.axvline(1612)
plt.ylim(0, 1)
plt.xlabel("Frequency (MHz)")
plt.ylabel("Illumination Efficiency")
plt.title("Illumination Efficiency for Different Geometries")
plt.legend()
plt.show()

# Radiation Patterns (Normalized Power)
plt.figure()
for data in sim_data:
    theta = data["theta"]
    e_plane = data["e_plane"]
    h_plane = data["h_plane"]
    max_e = data["max_e"]
    max_h = data["max_h"]
    plt.plot(theta, 20 * np.log10(np.abs(e_plane) / max_e), label=f"{data['label']} E-plane")
    plt.plot(theta, 20 * np.log10(np.abs(h_plane) / max_h), label=f"{data['label']} H-plane")
plt.axvline(dt_angle, color='grey')
plt.axvline(-dt_angle, color='grey')
plt.ylim(-45, 1)
plt.xlabel("Angle (deg.)")
plt.ylabel("Normalized Power (dB)")
plt.title("Radiation Patterns for Different Geometries")
plt.legend()
plt.show()
