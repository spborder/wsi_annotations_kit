"""

Text file conversion from JSON to XML

"""

import json
import lxml.etree as ET


# Input file will be "input_file.json"
# output file will be "output_file.xml"

input_file_name = 'C:\\Users\\Sam\\Downloads\\Glomeruli.json'
output_file_name = 'C:\\Users\\Sam\\Downloads\\output_file.xml'

with open(input_file_name,'r') as f:
    input_json = json.load(f)

def xml_add_annotation(xml,id,layer_name=None):

    if layer_name is None:
        layer_name = f'Layer_{id}'

    # For adding more annotation colors to the current list (this needs some work)
    new_annotation = ET.SubElement(xml,'Annotation',attrib={'Type':'4','Visible':'1','ReadOnly':'0',
                                                            'Incremental':'0','LineColorReadOnly':'0',
                                                            'Id':str(id),'NameReadOnly':'0','LayerName':layer_name})
    
    regions = ET.SubElement(new_annotation,'Regions')

def xml_add_region(poly_coords,xml,layer_id,region_id=None):

    annotation = xml.find(f'Annotation[@Id="{layer_id}"]')
    regions = annotation.find('Regions')

    if region_id is None:
        region_id = len(regions.findall('Region'))+1

    region = ET.SubElement(regions,'Region',attrib={'NegativeROA':'0','ImageFocus':'-1',
                                                    'DisplayId':'1','InputRegionId':'0',
                                                    'Analyze':'0','Type':'0','Id':str(region_id)})
    
    vertices = ET.SubElement(region,'Vertices')

    # Scaling according to source crs for this object (crs should just be the y,x top left coordinates)
    # Note: The int conversion here might fail for very very large slides >32bit

    for pt in poly_coords:
        ET.SubElement(vertices,'Vertex',attrib={'X':str(pt[1]),
                                                'Y':str(pt[0]),
                                                'Z':'0'})
    
    ET.SubElement(vertices,'Vertex',attrib={'X':str(poly_coords[0][1]),
                                            'Y':str(poly_coords[0][0]),
                                            'Z':'0'})
    

# Initial xml colors, add more if needed
output_xml = ET.Element('Annotations')

if not type(input_json)==list:
    input_json = [input_json]

print(type(input_json))

for ann_idx,ann in enumerate(input_json):
    print(ann)
    if 'name' in ann['annotation']:
        structure_name = ann['annotation']['name']
    else:
        structure_name = f'Structure_{ann_idx}'

    xml_add_annotation(output_xml,ann_idx+1,structure_name)
    
    if 'elements' in ann['annotation']:
        for obj_idx,obj in enumerate(ann['annotation']['elements']):

            coords = obj['points']

            xml_add_region(coords,output_xml,ann_idx+1)

xml_string = ET.tostring(output_xml,encoding='unicode',pretty_print=True)
with open(output_file_name,'w') as f:
    f.write(xml_string)





