import matplotlib.pyplot as plt
import numpy as np

# Chemin vers le fichier de géométrie
file_path = "C:\\Users\\BUNICE\\Desktop\\stage\\HORN-STOCKERT2\\hornStockert2.ge"

# Initialiser une liste pour stocker les coordonnées de la géométrie et du guide d'onde
guide_coordinates = []
geometry_coordinates = []

# Lire les lignes à partir du fichier
with open(file_path, 'r') as file:
    lines = file.readlines()

# Identifier le début de la section géométrie
start_idx = 0
for i, line in enumerate(lines):
    if "GEOMETRY" in line:
        start_idx = i + 1
        break

# Extraire spécifiquement les coordonnées du guide d'onde (1ère et 3ème lignes après "GEOMETRY")
guide_lines_indices = [0, 2]  # Indices pour la 1ère et la 3ème ligne
for i in guide_lines_indices:
    parts = lines[start_idx + i].split()
    if len(parts) >= 3:
        x = float(parts[1].replace('D', 'E'))
        z = float(parts[2].replace('D', 'E'))
        guide_coordinates.append((z, x))

# Traiter le reste de la géométrie
for line in lines[start_idx + 3:]:
    parts = line.split()
    if len(parts) >= 3 and parts[0].isdigit() and int(parts[0]) > 1:
        x = float(parts[1].replace('D', 'E'))
        z = float(parts[2].replace('D', 'E'))
        geometry_coordinates.append((z, x))

# Créer la symétrie par rapport à l'axe des Z
guide_coords = np.array(guide_coordinates)
geometry_coords = np.array(geometry_coordinates)
guide_coords_sym = np.copy(guide_coords)
guide_coords_sym[:, 1] = -guide_coords_sym[:, 1]
geometry_coords_sym = np.copy(geometry_coords)
geometry_coords_sym[:, 1] = -geometry_coords_sym[:, 1]

# Tracer la géométrie avec symétrie
plt.figure(figsize=(10, 5))
plt.plot(guide_coords[:, 0], guide_coords[:, 1], 'ro-', label='Guide d\'onde original')
plt.plot(guide_coords_sym[:, 0], guide_coords_sym[:, 1], 'ro-', label='Guide d\'onde symétrique')
plt.plot(geometry_coords[:, 0], geometry_coords[:, 1], 'bo-', label='Géométrie originale')
plt.plot(geometry_coords_sym[:, 0], geometry_coords_sym[:, 1], 'bo-', label='Géométrie symétrique')
plt.xlabel('Z (m)')
plt.ylabel('X (m)')
plt.legend()
plt.title('Géométrie complète de l\'antenne avec symétrie')
plt.grid(True)
plt.axis('equal')
plt.show()

# Charger les résultats de simulation
sim_freq = 1420
file = f"C:\\Users\\BUNICE\\Desktop\\stage\\HORN-STOCKERT2\\hornStockert2_{int(sim_freq*1000)}k.npz"
npzfile = np.load(file, allow_pickle=False)

# Extraction des variables nécessaires
theta = npzfile['theta']
e_plane = npzfile['e_plane']
h_plane = npzfile['h_plane']
max_e = np.max(np.abs(e_plane))
max_h = np.max(np.abs(h_plane))
theta_step = npzfile['theta_step']

# Calculs et affichage des résultats
f_D = 0.3
dt_angle = 180 * np.arctan2(f_D/2, f_D**2 - 1/16) / np.pi

# Afficher les patrons de rayonnement
fig, ax = plt.subplots()
ax.plot(theta, 20 * np.log10(np.abs(e_plane) / max_e), color="darkblue", label="E-plane")
ax.plot(-theta, 20 * np.log10(np.abs(e_plane) / max_e), color="darkblue")
ax.plot(theta, 20 * np.log10(np.abs(h_plane) / max_h), color="darkgreen", label="H-plane")
ax.plot(-theta, 20 * np.log10(np.abs(h_plane) / max_h), color="darkgreen")
ax.set_ylim(-45, 1)
ax.axvline(-dt_angle, color='grey')
ax.axvline(dt_angle, color='grey')
ax.set_title(f"Stockert horn (feed), {sim_freq} MHz")
ax.set_xlabel("angle (deg.)")
ax.set_ylabel("normalized power (dB)")
ax.legend()
plt.show()

# Calcul de la puissance gaussienne
taper = -17
alpha = -taper / 20 * np.log(10)
Nt = 2 * (1 - np.exp(-alpha))**2 / (alpha * (1 - np.exp(-2 * alpha)))

# Deuxième graphique : distribution de puissance normalisée avec profil gaussien
fig, ax = plt.subplots()
ax.plot(theta, np.abs(e_plane)**2 / np.max(np.abs(e_plane)**2), "-", color="darkblue", label="E-plane")
ax.plot(theta, np.abs(h_plane)**2 / np.max(np.abs(h_plane)**2), "-", color="brown", label="H-plane")
ax.plot(theta, np.exp(-alpha * (theta / dt_angle)**2 )**2, color="grey", label="Gaussian")
ax.axvline(dt_angle, color='grey', label='Gaussian cutoff angle')
ax.set_xlim(0, 90)
ax.set_title(f"Stockert horn (feed), {sim_freq} MHz")
ax.set_xlabel("angle (deg.)")
ax.set_ylabel("normalized power")
ax.legend()
plt.show()

rh_plane = np.abs(e_plane-1j*h_plane)/(np.sqrt(2)*max_e)
lh_plane = np.abs(e_plane+1j*h_plane)/(np.sqrt(2)*max_e)

fig, ax = plt.subplots()
# ax.plot(angles1, power1-power1[0], ".-", label="DT doc");
# ax.plot(angles2, power2, ".-", label="N2UO paper");
ax.plot(theta, 20*np.log10(np.abs(e_plane)/max_e), color="red", label="E-plane");
ax.plot(theta, 20*np.log10(np.abs(h_plane)/max_h), color="purple", label="H-plane");
#ax.plot(theta, 20*np.log10(rh_plane), color="orange", label="RHCP")
#ax.plot(theta, 20*np.log10(lh_plane), color="purple", label="LHCP")
ax.plot(theta, 20*np.log10( np.exp(-alpha*((theta/84.76))**2) ), color="grey", label="Gaussian")
ax.set_xlim(0,90)
ax.set_ylim(-21,0)
plt.axvline(x=dt_angle, color="darkblue");
ax.set_title("Stockert horn (feed), "+str(sim_freq)+" Mhz");
ax.set_xlabel("angle (deg.)")
ax.set_ylabel("normlized power (dB)")
ax.legend();
fig, ax = plt.subplots()
ax.plot(theta, 180*np.angle(e_plane)/np.pi, "-", color="darkblue", label="E-plane");
ax.plot(theta, 180*np.angle(h_plane)/np.pi, "-", color="brown", label="H-plane");
ax.axvline(dt_angle);
ax.set_xlim(0,90)
ax.set_ylim(-180,180)
ax.set_title("Stockert horn (feed), "+str(sim_freq)+" Mhz");
ax.set_xlabel("angle (deg.)")
ax.set_ylabel("phase (deg.)")
ax.legend();
def cotan(x):
    return -np.tan(x + np.pi/2)
def ill_eff(field, theta, theta_step_rad, open_angle):
    theta_edge_rad = np.pi*dt_angle/180

    full_sphere = np.pi;
    
    u1 = np.abs( np.sum(field[theta<=open_angle] * np.tan(theta[theta<=open_angle]/2) * theta_step_rad) )**2 
    u2 = np.sum( np.abs(field[theta<=full_sphere])**2 * np.sin(theta[theta<=full_sphere]) * theta_step_rad )

    ill = 2 * (cotan(theta_edge_rad/2)**2) * (u1/u2)
    return(ill)

ill = ill_eff(e_plane, np.pi*theta/180, np.pi*theta_step/180, np.pi*dt_angle/180)
print("E:",ill)

ill = ill_eff(h_plane, np.pi*theta/180, np.pi*theta_step/180, np.pi*dt_angle/180)
print("H:", ill)
field = e_plane[theta<=dt_angle]
theta_rad = np.pi*theta[theta<=dt_angle]/180
theta_step_rad = np.pi*theta_step/180

theta_full_rad = np.pi*theta[theta<180]/180
field_full = e_plane[theta<180]

u_dish = np.sum( np.abs(field)**2 * np.sin(theta_rad) * theta_step_rad )
u_full = np.sum( np.abs(field_full)**2 * np.sin(theta_full_rad) * theta_step_rad )

ns = u_dish/u_full

print("E:",ns)
print("Tspill:",(1-ns)*290)
f_range = np.arange(1200,1660,10)
print(f_range)
freqs = np.zeros(len(f_range))
gains = np.zeros(len(f_range))
directivity = np.zeros(len(f_range))
beamwidth = np.zeros( (len(f_range),3))
s11_mag = np.zeros(len(f_range))
s11_complex = np.zeros(len(f_range)).astype(complex)
nt_e = np.zeros(len(f_range))
nt_h = np.zeros(len(f_range))
nt_r = np.zeros(len(f_range))
nt_l = np.zeros(len(f_range))

for i in range(len(f_range)):
    f = f_range[i]
    file = f"C:\\Users\\BUNICE\\Desktop\\stage\\HORN-STOCKERT2\\hornStockert2_{int(f*1000)}k.npz"
    npzfile = np.load(file, allow_pickle=False)
    gains[i] = npzfile['gain']
    directivity[i] = npzfile['directivity']
    freqs[i] = npzfile['freq']
    s11_mag[i] = npzfile['s11_mag']
    s11_complex[i] = npzfile['s11']
    beamwidth[i,:] = npzfile['beamwidth']
    max_c = (np.max( np.abs(npzfile['e_plane']) + np.abs(npzfile['h_plane'])))
    r_plane = np.abs(npzfile['e_plane']+1j*npzfile['h_plane'])/(np.sqrt(2)*max_c)
    l_plane = np.abs(npzfile['e_plane']-1j*npzfile['h_plane'])/(np.sqrt(2)*max_c)
    nt_e[i] = ill_eff(npzfile['e_plane'],np.pi*theta/180, np.pi*theta_step/180, np.pi*dt_angle/180)
    nt_h[i] = ill_eff(npzfile['h_plane'],np.pi*theta/180, np.pi*theta_step/180, np.pi*dt_angle/180)
    nt_r[i] = ill_eff(r_plane,np.pi*theta/180, np.pi*theta_step/180, np.pi*dt_angle/180)
    nt_l[i] = ill_eff(l_plane,np.pi*theta/180, np.pi*theta_step/180, np.pi*dt_angle/180)
fig, ax = plt.subplots()
ax.plot(freqs/1e6, 20*np.log10(s11_mag), "-", color="darkblue");
ax.set_xlabel("freq. (Mhz)")
ax.set_ylabel("S11 (dB)");
ax.axvline(1300)
ax.axvline(1420);
ax.axvline(1612);
vswr = (1+s11_mag)/(1-s11_mag)

fig, ax = plt.subplots()
ax.plot(freqs/1e6, vswr, "-", color="darkblue");
ax.set_xlabel("freq. (Mhz)")
ax.set_ylabel("VSWR");
fig, ax = plt.subplots()
ax.plot(freqs/1e6, gains, ".-", color="darkblue");
ax.plot(freqs/1e6, directivity, ".-", color="grey");
ax.axvline(1300)
ax.axvline(1420)
ax.axvline(1612);
ax.set_xlabel("freq. (Mhz)")
ax.set_ylabel("gain (dB)");
ax.set_title("Stockert horn (feed)");
fig, ax = plt.subplots()
ax.plot(freqs/1e6, beamwidth[:,0], ".-", color="darkblue", label="-3 dB");
ax.plot(freqs/1e6, beamwidth[:,1], ".-", color="lightblue", label="-10 dB");
ax.plot(freqs/1e6, beamwidth[:,2], ".-", color="purple", label="-20 dB");
ax.axhline(dt_angle, color="black");
ax.axvline(1300)
ax.axvline(1420)
ax.axvline(1612);
ax.set_ylim(0,100)
ax.set_xlabel("freq. (Mhz)")
ax.set_ylabel("Beamwidth (deg.)");
ax.set_title("Stockert horn (feed)");
ax.legend();
fig, ax = plt.subplots()
ax.plot(freqs/1e6, nt_e, ".-", color="darkblue", label="E-plane");
ax.plot(freqs/1e6, nt_h, ".-", color="lightblue", label="H-plane");
# ax.plot(freqs/1e6, nt_r, ".-", color="green", label="RHCP-plane");
# ax.plot(freqs/1e6, nt_l, ".-", color="lightgreen", label="LHCP-plane");
ax.plot(freqs/1e6, (nt_h+nt_e)/2, ".-", color="orange", label="E/H-plane");
ax.legend()
ax.axvline(1300)
ax.axvline(1420)
ax.axvline(1612);
ax.set_ylim(0,1)
ax.set_xlabel("freq. (Mhz)")
ax.set_ylabel("illumination efficiency")
ax.set_title("Stockert horn (feed)");


