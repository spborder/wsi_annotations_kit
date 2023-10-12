import json
import numpy as np
import girder_client

from wsi_annotations_kit import wsi_annotations_kit as wak
from shapely.geometry import Polygon, Point

#NAMES = ['non-gs-glomerulus', 'gs-glomerulus']

# Open and load the JSON data from the file
with open('VAN0016-LK-208-2-PAS-glomerulus-annotations.json', 'r') as json_file:
    data = json.load(json_file)

# Initialize spot_annotationprints
spot_annotations = wak.Annotation()
spot_annotations.add_names(["Glomeruli"])


# Iterate through each "geometry" object and extract the "coordinates"
for feature in data:
    coordinates = feature["geometry"]["coordinates"][0]
    scaled_coordinates = []
    #coordinates = [tuple(coord) for coord in coordinates]

    """ for coord_set in coordinates:
        for coord in coord_set:
            # Assuming coord[0] is the y-coordinate and coord[1] is the x-coordinate
            scaled_y = int(coord[0])
            scaled_x = int(coord[1])
            scaled_coordinates.append((scaled_x, scaled_y)) """

    # Create a Polygon from the scaled coordinates
    spot_poly = Polygon(coordinates)

    # Add the Polygon as a shape to spot_annotations
    spot_annotations.add_shape(
        poly=spot_poly,
        box_crs=[0, 0],
        structure="Glomeruli",  
        name=None,
        properties=None
    )

# Create a Histomics object from spot_annotations
annot = wak.Histomics(spot_annotations)
annotations_json = json.dumps(annot.json, indent=4)

output_file_path = 'annot.json'

# Write the JSON-formatted string to the file
with open(output_file_path, 'w') as output_file:
    output_file.write(annotations_json)

# You can now use 'annot' for further processing or export it as needed
        
#girder_folder_id = '647f32c9435c92704a565d1a'
itemID = '6522f8472a1551572931e58e'
girderApiUrl = "http://athena.rc.ufl.edu/api/v1/"
gc = girder_client.GirderClient(apiUrl=girderApiUrl)
girderToken = gc.get('token/session')['token']
gc.setToken(girderToken)
response = gc.post(f'annotation/item/{itemID}', data=json.dumps(annot.json), headers={'X-HTTP-Method': 'POST','Content-Type':'application/json'})   
print(response)         




