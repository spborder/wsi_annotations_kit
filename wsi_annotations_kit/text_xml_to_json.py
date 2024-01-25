"""

Text file conversion from XML to JSON

"""

import json
import lxml.etree as ET
import uuid

# Input file will be input_file.xml
# Output file will be output_file.json

input_file_name = 'input_file.xml'
output_file_name = 'output_file.json'

def add_element(output_file,ann_idx,coords):

    element_dict = {
        'type':'polyline',
        'points':[i+[0] for i in coords],
        'id':uuid.uuid4().hex[:24],
        'closed':True
    }

    output_file[ann_idx]['elements'].append(element_dict)



tree = ET.parse(input_file_name)
structures_in_xml = tree.getroot().findall('Annotation')
output_file = []    
for ann_idx in range(0,len(structures_in_xml)):

    this_structure = tree.getroot().findall(f'Annotation[@Id="{str(ann_idx+1)}]/Regions/Region')

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
                int(vert.attrib['X']),
                int(vert.attrib['Y'])
            ])
        
        add_element(output_file,ann_idx,coords)







