import numpy as np
import matplotlib.pyplot as plt

def read_coordinates(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    start_idx = 0
    geometry_coordinates = []
    waveguide_coordinates = []
    
    # Find the start of the geometry section
    for i, line in enumerate(lines):
        if "GEOMETRY" in line:
            start_idx = i + 1
            break
    
    # Extract waveguide coordinates
    waveguide_line_indices = [0, 2]
    for i in waveguide_line_indices:
        parts = lines[start_idx + i].split()
        if len(parts) >= 3:
            x = float(parts[1].replace('D', 'E'))
            z = float(parts[2].replace('D', 'E'))
            waveguide_coordinates.append((z, x))
    
    # Extract the rest of the geometry
    for line in lines[start_idx + 3:]:
        parts = line.split()
        if len(parts) >= 3 and parts[0].isdigit() and int(parts[0]) > 1:
            x = float(parts[1].replace('D', 'E'))
            z = float(parts[2].replace('D', 'E'))
            geometry_coordinates.append((z, x))
    
    return waveguide_coordinates, geometry_coordinates

def parabola(y, focal_length, c):
    return (-1 / (4 * focal_length)) * y**2 + c

# Read coordinates from the file
file_path = "C:\\Users\\BUNICE\\Desktop\\stage\\PARA\\para-horne-stockert-dwing.ge"
waveguide_coords, geometry_coords = read_coordinates(file_path)

reflector_coords = geometry_coords[-2:]
vertex = reflector_coords[0]
high_point = reflector_coords[1]
low_point = (high_point[0], -high_point[1])

c = vertex[0]
focal_length_min = 11000  # Minimum focal length
focal_length_max = 13000  # Maximum focal length

# Generate focal length values within the range with a step of 100 mm
focal_lengths = np.arange(focal_length_min, focal_length_max + 100, 100)

# Create a plot for each focal length
for focal_length in focal_lengths:
    plt.figure(figsize=(15, 8))

    # Plot the waveguide and the geometry
    for coords, color, label in [(waveguide_coords, 'r', 'Waveguide'), (geometry_coords[:-2], 'b', 'Geometry')]:
        coords = np.array(coords)
        plt.plot(coords[:, 0], coords[:, 1], marker='o', linestyle='-', color=color, label=label)
        plt.plot(coords[:, 0], -coords[:, 1], marker='o', linestyle='-', color=color)  

    # Generate y values for plotting the parabola
    y_values = np.linspace(high_point[1], low_point[1], 1000)
    x_values = parabola(y_values, focal_length, c)

    # Plot the parabola of the reflector
    plt.plot(x_values, y_values, linestyle='-', label=f'Focal Length: {int(focal_length)} mm')

    # Titles and labels
    plt.title(f'Antenna Geometry with Parabolic Reflector Profile (Focal Length: {int(focal_length)} mm)')
    plt.xlabel('Z (m)')
    plt.ylabel('X (m)')
    plt.legend()

    # Display the graph
    plt.grid(True)
    plt.axis('equal')
    plt.show()
