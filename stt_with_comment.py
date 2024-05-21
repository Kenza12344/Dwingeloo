import matplotlib.pyplot as plt  # Importing the library for graph visualization
import numpy as np  # Importing the library for numerical array manipulation

# Path to the file containing the antenna geometry data
file_path = "C:\\Users\\BUNICE\\Desktop\\stage\\HORN-STOCKERT2\\hornStockert2.ge"

# Initializing lists to store the waveguide and geometry coordinates
guide_coordinates = []
geometry_coordinates = []

# Reading data from the file
with open(file_path, 'r') as file:
    lines = file.readlines()  # Reading all lines from the file and storing them in a list

# Finding the start index of the 'GEOMETRY' section in the file
start_idx = 0
for i, line in enumerate(lines):
    if "GEOMETRY" in line:  # Checking for the presence of the keyword "GEOMETRY"
        start_idx = i + 1  # Saving the index of the next line
        break

# Extracting waveguide coordinates from specified lines after "GEOMETRY"
guide_lines_indices = [0, 2]  # Line indices to use for the waveguide
for i in guide_lines_indices:
    parts = lines[start_idx + i].split()  # Splitting the line into parts
    if len(parts) >= 3:  # Ensuring the line has at least three elements
        x = float(parts[1].replace('D', 'E'))  # Converting Fortran D format to standard Python E format for floats
        z = float(parts[2].replace('D', 'E'))
        guide_coordinates.append((z, x))  # Adding coordinates (inverting X and Z) to the list

# Extracting geometry coordinates for the rest of the lines
for line in lines[start_idx + 3:]:
    parts = line.split()  # Splitting the line into parts
    if len(parts) >= 3 and parts[0].isdigit() and int(parts[0]) > 1:  # Checking that the line is valid for extraction
        x = float(parts[1].replace('D', 'E'))  # Same conversion as above
        z = float(parts[2].replace('D', 'E'))
        geometry_coordinates.append((z, x))  # Adding inverted coordinates to the list

# Creating symmetrical coordinates for plotting
guide_coords = np.array(guide_coordinates)  # Converting the list into a numpy array
geometry_coords = np.array(geometry_coordinates)
guide_coords_sym = np.copy(guide_coords)  # Copying the array for symmetry
guide_coords_sym[:, 1] = -guide_coords_sym[:, 1]  # Inverting X coordinates for symmetry
geometry_coords_sym = np.copy(geometry_coords)
geometry_coords_sym[:, 1] = -geometry_coords_sym[:, 1]

# Plotting the complete antenna geometry with symmetry
plt.figure(figsize=(10, 5))  # Setting the figure size
plt.plot(guide_coords[:, 0], guide_coords[:, 1], 'ro-', label='Waveguide')  # Plotting waveguide coordinates
plt.plot(guide_coords_sym[:, 0], guide_coords_sym[:, 1], 'ro-')  # Plotting the symmetry of the waveguide
plt.plot(geometry_coords[:, 0], geometry_coords[:, 1], 'bo-', label='Geometry')  # Plotting geometry coordinates
plt.plot(geometry_coords_sym[:, 0], geometry_coords_sym[:, 1], 'bo-')  # Plotting the symmetry of the geometry
plt.xlabel('Z (m)')  # X-axis label
plt.ylabel('X (m)')  # Y-axis label
plt.legend()  # Displaying the legend
plt.title('Complete Antenna Geometry with Symmetry')  # Graph title
plt.grid(True)  # Enabling the grid
plt.axis('equal')  # Scaling the axes for a uniform ratio
plt.show()  # Displaying the graph


