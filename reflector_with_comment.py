import numpy as np  # Importe la bibliothèque numpy pour les opérations numériques
import matplotlib.pyplot as plt  # Importe la bibliothèque matplotlib pour la création de graphiques

# Définit le chemin d'accès au fichier .ge contenant les données de géométrie
file_path = "C:\\Users\\BUNICE\\Desktop\\stage\\para-horne-stockert\\para-horn-stockert.ge"

# Définit une fonction pour lire les coordonnées de géométrie à partir du fichier
def read_geometry_coordinates(file_path):
    with open(file_path, 'r') as file:  # Ouvre le fichier en mode lecture
        lines = file.readlines()  # Lit toutes les lignes du fichier et les stocke dans une liste

    # Initialisation de l'index pour trouver le début de la section géométrie
    start_idx = 0
    for i, line in enumerate(lines):  # Itère sur chaque ligne avec son indice
        if "GEOMETRY" in line:  # Vérifie si le mot "GEOMETRY" est présent dans la ligne
            start_idx = i + 1  # Définit l'index de début de la section géométrie
            break  # Arrête la boucle après avoir trouvé le début de la section géométrie

    # Initialisation de la liste pour stocker les coordonnées de la géométrie
    geometry_coordinates = []
    for line in lines[start_idx + 3:]:  # Itère à partir de trois lignes après le début de la section géométrie
        parts = line.split()  # Divise chaque ligne en plusieurs parties basées sur les espaces
        if len(parts) >= 3 and parts[0].isdigit() and int(parts[0]) > 1:  # Vérifie si la ligne contient des coordonnées valides
            x = float(parts[1].replace('D', 'E'))  # Remplace 'D' par 'E' pour convertir en format float (not. scientifique)
            z = float(parts[2].replace('D', 'E'))  # Idem pour la coordonnée z
            geometry_coordinates.append((z, x))  # Ajoute un tuple de coordonnées à la liste
    return geometry_coordinates  # Retourne la liste des coordonnées

# Appelle la fonction pour lire les coordonnées à partir du chemin de fichier spécifié
geometry_coordinates = read_geometry_coordinates(file_path)

# Sélectionne les deux dernières coordonnées qui représentent le réflecteur
reflector_coordinates = geometry_coordinates[-2:]

# Affiche les coordonnées extraites pour le réflecteur
print("Coordonnées du réflecteur extraites :")
print(reflector_coordinates)

# Décompose les coordonnées du réflecteur en sommet, point haut et point bas
sommet = reflector_coordinates[0]
point_haut = reflector_coordinates[1]
point_bas = (point_haut[0], -point_haut[1])  # Crée un point bas en inversant le signe de la coordonnée y du point haut

# Calcule le coefficient 'a' de l'équation parabolique (x = ay^2 + by + c)
c = sommet[0]  # Le terme constant c est la coordonnée x du sommet
a = (point_haut[0] - c) / point_haut[1]**2  # Calcule le coefficient a basé sur le sommet et le point haut

# Définit une fonction pour calculer les valeurs x de la parabole basée sur y
def parabole(y, a, c):
    return a * y**2 + c  # Retourne la valeur de x en utilisant l'équation parabolique

# Génère un ensemble de valeurs y pour tracer la parabole, entre le point haut et bas
y_values = np.linspace(point_bas[1], point_haut[1], 1000)  # Crée 1000 points linéairement espacés entre y bas et y haut
x_values = parabole(y_values, a, c)  # Calcule les correspondants x pour les valeurs y générées

# Configure le graphique pour la visualisation de la parabole
plt.figure(figsize=(10, 8))  # Définit la taille du graphique
plt.plot(x_values, y_values, label='Parabole')  # Trace la parabole
plt.scatter(*sommet, color='red', zorder=5, label='Sommet')  # Ajoute le sommet sur le graphique
plt.scatter(*point_haut, color='green', zorder=5, label='Point haut')  # Ajoute le point haut sur le graphique
plt.scatter(*point_bas, color='blue', zorder=5, label='Point bas')  # Ajoute le point bas sur le graphique

# Ajoute des légendes et des détails au graphique
plt.title('Réflecteur parabolique')  # Définit le titre du graphique
plt.xlabel('X')  # Définit le label de l'axe des x
plt.ylabel('Y')  # Définit le label de l'axe des y
plt.grid(True)  # Active la grille pour une meilleure visibilité
plt.legend()  # Affiche les légendes
plt.axis('equal')  # Assure que les axes sont à échelle égale

# Affiche le graphique
plt.show()  # Affiche le graphique à l'écran




