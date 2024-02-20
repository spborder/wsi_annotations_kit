"""

Text file conversion from XML to JSON

"""

import json
import lxml.etree as ET
import uuid

# Input file will be input_file.xml
# Output file will be output_file.json

input_file_name = 'C:\\Users\\Sam\\Downloads\\output_file.xml'
output_file_name = 'C:\\Users\\Sam\\Downloads\\output_file_2.json'

def add_element(coords):

    element_dict = {
        'type':'polyline',
        'points':[i+[0] for i in coords],
        'id':uuid.uuid4().hex[:24],
        'closed':True
    }

    return element_dict

tree = ET.parse(input_file_name)
structures_in_xml = tree.getroot().findall('Annotation')
output_file = []    
for ann_idx in range(0,len(structures_in_xml)):

    this_structure = tree.getroot().findall(f'Annotation[@Id="{str(ann_idx+1)}"]/Regions/Region')

    output_file.append(
        {
            'annotation':{
                'name':f'Layer_{ann_idx}',
                'elements':[]
            }
        }
    )

    for obj in this_structure:
        vertices = obj.findall('./Vertices/Vertex')
        coords = []
        for vert in vertices:
            coords.append([
                int(float(vert.attrib['X'])),
                int(float(vert.attrib['Y']))
            ])
        
        output_file[ann_idx]['annotation']['elements'].append(add_element(coords))

with open(output_file_name,'w') as f:
    json.dump(output_file,f)





