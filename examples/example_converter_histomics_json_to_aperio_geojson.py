import json
import math

def convert_to_qupath_geojson(input_filename, output_filename):
    """
    Converts a JSON annotation format to QuPath-compatible GeoJSON format.
    Args:
        input_filename (str): Path to the input JSON annotation file.
        output_filename (str): Path to save the output GeoJSON file.
    """
    try:
        with open(input_filename, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: File '{input_filename}' not found.")
        return
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from file '{input_filename}'.")
        return
      
    # Initialize the structure for a GeoJSON FeatureCollection
    geojson_output = {
        "type": "FeatureCollection",
        "features": []
    }
    
    # The input data is a list of annotations
    for annotation_item in data:
        # Check if the required keys exist
        if "annotation" not in annotation_item or "elements" not in annotation_item["annotation"]:
            print(f"Skipping item due to missing keys: {annotation_item}")
            continue

       # Iterate through each geometric element in the annotatio
        for element in annotation_item['annotation']['elements']:
            geometry_type = element.get('type', '')
            coordinates = []
            if geometry_type == 'polyline' and 'points' in element:
                # Extract coordinates for the polyline
                coordinates = [point[:2] for point in element['points']]
                # Ensure the polygon is closed by making sure the first and last points are the same
                if element.get('closed') and coordinates and coordinates[0] != coordinates[-1]:
                    coordinates.append(coordinates[0])  # Close the polygon

            elif geometry_type == 'ellipse' and 'center' in element and 'radiusX' in element and 'radiusY' in element:
                
                # Convert ellipse to polygon with vertices
                center = element['center']
                radiusX = element['radiusX']
                radiusY = element['radiusY']
                num_points = 50  # Number of vertices to approximate the ellipse
                # Generate points to outline shape of an ellipse using trigonometric functions(parametric equation of an ellipse)
                coordinates = [
                    [
                        center[0] + radiusX * math.cos(2 * math.pi * i / num_points),
                        center[1] + radiusY * math.sin(2 * math.pi * i / num_points)
                    ]
                    for i in range(num_points)
                ]
                coordinates.append(coordinates[0])  # Close the polygon
            else:
                print(f"Unsupported geometry type: {geometry_type}")
                continue

            # Validate the number of points in the polygon
            if len(coordinates) < 4:
                print(f"Skipping feature with insufficient points: {coordinates}")
            continue
               geometry = {
                "type": "Polygon",
                "coordinates": [coordinates]  # Linear ring for a polygon
            }

            properties = {
                "classification": annotation_item['annotation'].get('name', 'DefaultClassification'),
                "name": annotation_item['annotation'].get('name', 'Unnamed'),
                "id": element.get('id', ''),
                "fillColor": element.get('fillColor', '#FFFFFF'),
                "lineColor": element.get('lineColor', '#000000'),
                "lineWidth": element.get('lineWidth', 1)
            }

            feature = {
                "type": "Feature",
                "geometry": geometry,
                "properties": properties
            }

            geojson_output["features"].append(feature)

    try:
        with open(output_filename, 'w') as f:
            json.dump(geojson_output, f, indent=2)
        print(f"Successfully converted '{input_filename}' to '{output_filename}'.")
    except IOError:
        print(f"Error: Could not write to file '{output_filename}'.")


# Example Usage
convert_to_qupath_geojson(r"file_path.json",
                          r"your_result.geojson")
