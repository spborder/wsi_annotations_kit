"""

Codes for generating, saving, and converting annotations for segmentations on whole slide images

"""

import os
import sys
import numpy as np

import lxml.etree as ET
import geojson
from geojson import Feature, dump
import json

from tqdm import tqdm

import shapely
from shapely.geometry import box, Polygon, shape
from skimage.measure import label, find_contours
import uuid


class Object:
    def __init__(self,
                 shape_structure,
                 location,
                 structure_name,
                 name):
    
        self.poly = shape_structure
        self.crs = location
        self.structure = structure_name
        self.name = name

        # If not specified just pick something general
        if self.structure is None:
            self.structure = 'Structure'
        if self.name is None:
            self.name = uuid.uuid4().hex[:24]


class AperioXML:
    def __init__(self,
                 annotations,
                 layer_names=None,
                 mpp=None,
                 verbose = True
                 ):
        
        self.annotations = annotations
        self.layer_names = layer_names
        self.mpp = mpp

        if self.layer_names is None:
            self.layer_names = {i:j for i,j in zip(list(self.annotations.objects.keys()),list(range(1,1+len(list(self.annotations.objects.keys())))))}
        
        self.xml_colors = [65280,65535,33023,255,16711680]

        if verbose:
            pbar = tqdm(list(self.annotations.objects.keys()))

        self.xml_create()
        for n_idx,n in enumerate(self.annotations.objects):

            if verbose:
                pbar.update(n_idx)
                pbar.set_description(f'Converting to Aperio XML format, on: {n}, {len(self.annotations.objects[n])} found')

            # Adding new layer to current annotations
            self.xml_add_annotation(id=self.layer_names[n],layer_name = n)
            for o in self.annotations.objects[n]:
                self.xml_add_region(o)

        if verbose:
            pbar.close()

    def xml_create(self):
        if not self.mpp is None:
            self.xml = ET.Element('Annotations',attrib={'MicronsPerPixel':str(self.mpp)})
        else:
            self.xml = ET.Element('Annotations')

    def xml_add_annotation(self,id,layer_name=None):

        if layer_name is None:
            layer_name = f'Layer_{id}'

        # For adding more annotation colors to the current list (this needs some work)
        if id>=len(self.xml_colors):
            self.xml_colors += [i+10 for i in self.xml_colors]
        new_annotation = ET.SubElement(self.xml,'Annotation',attrib={'Type':'4','Visible':'1','ReadOnly':'0',
                                                               'Incremental':'0','LineColorReadOnly':'0',
                                                               'LineColor':str(self.xml_colors[id-1]),
                                                               'Id':str(id),'NameReadOnly':'0','LayerName':layer_name})
        
        regions = ET.SubElement(new_annotation,'Regions')

    def xml_add_region(self,obj,region_id=None):

        poly = obj.poly
        crs = obj.crs
        name = obj.name
        id = self.layer_names[obj.structure]

        annotation = self.xml.find(f'Annotation[@Id="{id}"]')
        regions = annotation.find('Regions')

        if region_id is None:
            region_id = len(regions.findall('Region'))+1
        
        if name is None:
            name = ''

        region = ET.SubElement(regions,'Region',attrib={'NegativeROA':'0','ImageFocus':'-1',
                                                        'DisplayId':'1','InputRegionId':'0',
                                                        'Analyze':'0','Type':'0','Id':str(region_id),
                                                        'Text':name})
        
        vertices = ET.SubElement(region,'Vertices')

        # Outputting coordinates for this object
        poly_coords = list(poly.exterior.coords)

        # Scaling according to source crs for this object (crs should just be the y,x top left coordinates)
        # Note: The int conversion here might fail for very very large slides >32bit
        scaled_poly_coords = [[int(i[0]+crs[0]),int(i[1]+crs[1])] for i in poly_coords]

        for pt in scaled_poly_coords:
            ET.SubElement(vertices,'Vertex',attrib={'X':str(pt[1]),
                                                    'Y':str(pt[0]),
                                                    'Z':'0'})
        
        ET.SubElement(vertices,'Vertex',attrib={'X':str(scaled_poly_coords[0][1]),
                                                'Y':str(scaled_poly_coords[0][0]),
                                                'Z':'0'})
        
    
class GeoJSON:
    def __init__(self,
                 annotations,
                 verbose = True):
        
        self.annotations = annotations

        if verbose:
            pbar = tqdm(list(self.annotations.objects.keys()))

        self.geojson_create()
        for n_idx,n in enumerate(self.annotations.objects):
            
            if verbose:
                pbar.update(n_idx)
                pbar.set_description(f'Converting to GeoJSON on: {n}, {len(self.annotations.objects[n])} found')

            for o in self.annotations.objects[n]:
                self.geojson_add_region(o)

        if verbose:
            pbar.close()

    def geojson_create(self):
        self.geojson = {'type':'FeatureCollection','features':[]}
    
    def geojson_add_region(self, obj):
        
        # Shifting object into global crs
        poly_coords = list(obj.poly.exterior.coords)
        scaled_poly_coords = [(int(i[0]+obj.crs[0]),int(i[1]+obj.crs[1])) for i in poly_coords]

        new_poly = Polygon(scaled_poly_coords)


        self.geojson['features'].append(
            Feature(geometry=new_poly,properties={'label':obj.name,'structure':obj.structure})
        )
    

class Histomics:
    def __init__(self,
               annotations,
               verbose = True):
        
        self.annotations = annotations

        if verbose:
            pbar = tqdm(list(self.annotations.objects.keys()))

        self.json_create()
        for n_idx,n in enumerate(self.annotations.objects):
            structure_dict = {'name':n,'attributes':{},'elements':[]}

            if verbose:
                pbar.update(n_idx)
                pbar.set_description(f'Converting to Histomics format, on: {n}, {len(self.annotations.objects[n])} found')

            for o in self.annotations.objects[n]:
                structure_dict['elements'].append(self.json_add_region(o))
            self.json.append(structure_dict)

        if verbose:
            pbar.close()

    def json_create(self):
        self.json = []

    def json_add_region(self,obj):

        # Shifting object into global crs 
        poly_coords = list(obj.poly.exterior.coords)
        scaled_poly_coords = [(int(i[0]+obj.crs[0]),int(i[1]+obj.crs[1])) for i in poly_coords]

        new_poly = Polygon(scaled_poly_coords)

        new_struct_dict = {
            'type':'polyline',
            'points':list(new_poly.exterior.coords),
            'id':uuid.uuid4().hex[:24],
            'closed':True,
            'user': {
                'name': obj.name,
                'structure': obj.structure
            }
        }

        return new_struct_dict


class Annotation:
    def __init__(self,
                 mpp = None,
                 min_size = None):

        self.objects = {}
        self.structure_names = []

        self.mpp = mpp
        self.min_size = min_size
    
    def __str__(self):
        
        if len(list(self.objects.keys()))==0:
            print('Empty annotation object')
        else:
            print(f'Annotation object containing: {len(list(self.objects.keys()))}')
            for n in self.objects:
                print(f'{n}: {len(self.objects[n])}')

    def add_shape(self, poly, box_crs, structure = None, name = None):
        
        # poly = shapely polygon or other geometry
        # box_crs = upper left corner for location that this object is in 
        # structure = structure name (which structure is this in general)
        # name = shape name (individual object name)
        if structure not in self.objects:
            self.objects[structure] = [Object(poly,box_crs,structure,name)]
        else:
            self.objects[structure] = [Object(poly,box_crs,structure,name)]
    
    def add_names(self,names):
        
        # Initializing 
        for n in names:
            self.objects[n] = {}
            self.structure_names.append(n)
            
    def xml_save(self,filename, layer_ids=None):

        xml_annotations = AperioXML(self,layer_ids,self.mpp)
        xml_string = ET.tostring(xml_annotations.xml,encoding='unicode',pretty_print = True)
        with open(filename,'w') as f:
            f.write(xml_string)
            f.close()

    def geojson_save(self,filename):

        geojson_annotations = GeoJSON(self)
        with open(filename,'w') as f:
            dump(geojson_annotations.geojson,f)        
            f.close()

    def json_save(self,filename):

        json_annotations = Histomics(self)
        with open(filename,'w') as f:
            dump(json_annotations.json,f)
            f.close()

    def add_mask(self,mask,box_crs,mask_type,structure = None):

        # Adding a mask object to your set of annotations
        # type = either 'one-hot' or nothing
        # structure here should be a list or dictionary for aligning index/label with structure
        if mask_type == 'one-hot':
            # Expecting mask format [height, width, classes]
            for cls in range(np.shape(mask)[-1]):
                
                if not structure is None:
                    if type(structure) == list:
                        structure_name = structure[cls]
                    elif type(structure) == dict:
                        structure_name = structure[str(cls)]
                    elif type(structure) == str:
                        structure_name = structure
                    else:
                        # Add some raise statement here
                        raise ValueError
                    
                class_mask = mask[:,:,cls].copy()

                # Using label to get number of objects 
                labeled_mask, n_objects = label(class_mask,background=0,return_num=True)
                for i in range(n_objects):
                    
                    # Find contours where labeled mask is equal to i
                    obj_contours = find_contours(labeled_mask,i)

                    # This is in (rows,columns) format
                    poly_list = [(i[0],i[1]) for i in obj_contours]
                    # Making polygon from contours
                    obj_poly = Polygon(poly_list)
                    self.add_shape(obj_poly,box_crs,structure_name)
        else:
            # Expecting mask format [height, width]
            if type(structure)==list:
                n_struct = len(structure)
            elif type(structure)==dict:
                n_struct=len(list(structure.keys()))
            elif type(structure)==str:
                n_struct = 1
            else:
                raise ValueError
            
            for cls in range(n_struct):

                # Assuming background is set to zero
                if type(structure)==list:
                    struct_idx = cls+1
                    structure_name = structure[cls]
                elif type(structure)==dict:
                    struct_idx = structure[list(structure.keys())[cls]]
                    structure_name = list(structure.keys())[cls]
                elif type(structure)==str:
                    struct_idx = 1
                    structure_name = structure

                class_mask = (mask.copy()==struct_idx)

                # Using label to get number of objects
                labeled_mask, n_objects = label(class_mask,background=0,return_num=True)
                for i in range(n_objects):

                    # Find contours where labeled mask is equal to i
                    obj_contours = find_contours(labeled_mask==i+1)
                    if len(obj_contours)==1:
                        poly_list = obj_contours[0].tolist()
                        poly_list = [(int(i[0]),int(i[1])) for i in poly_list]
                    else:
                        # Find the largest one
                        contours_size = [np.shape(i)[0] for i in obj_contours]
                        poly_list = obj_contours[np.argmax(contours_size)].tolist()
                        poly_list = [(int(i[0]),int(i[1])) for i in poly_list]

                    if len(poly_list)>2:
                        # Making polygon from contours 
                        obj_poly = Polygon(poly_list)
                        self.add_shape(obj_poly,box_crs,structure_name)



class Converter:
    def __init__(self,
                 starting_file: str,
                 ann_dict: dict,
                 verbose = True):

        self.starting_file = starting_file
        self.file_ext = self.starting_file.split('.')[-1]
        self.verbose = verbose

        self.invalid_count = 0

        # ann_dict = dictionary like {'structure_name':id} 
        # where id can either be layer id for aperio, index for json, or the property key that has the structure name for geojson
        self.ann_dict = ann_dict

        self.annotation = self.ingest_annotations()

    def ingest_annotations(self):

        self.annotation = Annotation()
        self.annotation.add_names(list(self.ann_dict.keys()))

        if self.file_ext=='xml':
            
            # Verboseness just for number of layers for xmls
            if self.verbose:
                pbar = tqdm(list(self.ann_dict.keys()))

            # Add MPP here if there is one

            tree = ET.parse(self.starting_file)
            for idx, structure in enumerate(self.ann_dict):

                structures_in_xml = tree.getroot().findall(f'Annotations[@Id="{str(self.ann_dict[structure])}"]/Regions/Region')

                if self.verbose:
                    pbar.update(idx)
                    pbar.set_description(f'Working on: {structure}: {len(structures_in_xml)} found')
                
                for struct_idx,region in enumerate(structures_in_xml):
                    vertices = region.findall('./Vertices/Vertex')

                    coords = []
                    for vert in vertices:
                        pixel_coords = (np.float32(vert.attrib['X']),np.float32(vert.attrib['Y']))

                        coords.append(pixel_coords)

                    if len(coords)>2:
                        if 'Text' in region.attrib:
                            if not region.attrib['Text']=='':
                                name = region.attrib['Text']
                            else:
                                name = None
                        else:
                            name = None
                        
                        region_poly = Polygon(coords)

                        # Checking if this polygon is valid
                        checked_poly = self.check_validity(region_poly)

                        # Adding to current annotations if valid
                        if not checked_poly is None:
                            self.annotation.add_shape(
                                poly = checked_poly,
                                crs = [0,0],
                                structure = structure,
                                name = name)

        if self.file_ext == 'geojson':
            # Loading initial annotations 
            with open(self.starting_file) as f:
                geojson_polys = geojson.load(f)

            # Verboseness for each feature in geojson
            if self.verbose:
                pbar = tqdm(geojson_polys['features'])
            
            for f_idx,f in enumerate(geojson_polys['features']):

                if self.verbose:
                    pbar.update(f_idx)
                    pbar.set_description(f'On Feature: {f_idx} of {len(geojson_polys["features"])}')
                
                if 'structure' in f['properties']:
                    structure_name = f['properties']['structure']
                else:
                    structure_name = 'Structure'

                if 'name' in f['properties']:
                    name = f['properties']['name']
                else:
                    name = None
                
                self.annotation.add_shape(
                    poly = shape(f['geometry']),
                    box_crs = [0,0],
                    structure = structure_name,
                    name = name
                )

        if self.file_ext == 'json':

            with open(self.starting_file) as f:
                json_annotations = json.load(f)

            # Verboseness for number of structures
            if self.verbose:
                pbar = tqdm(json_annotations)

            for st_idx,st in enumerate(json_annotations):
                structure = st['name']

                if self.verbose:
                    pbar.update(st_idx)
                    pbar.set_description(f'Working on: {structure}, found: {len(st["elements"])}')

                for o in st['elements']:

                    coordinates = o['points']
                    coordinates = [(i[0],i[1]) for i in coordinates]

                    if 'user' in o:
                        if 'name' in o['user']:
                            name = o['user']['name']
                        else:
                            name = None
                    else:
                        name = None
                    
                    # Probably don't need to check these? 
                    region_poly = Polygon(coordinates)
                    checked_poly = self.check_validity(region_poly)

                    if not checked_poly is None:
                        self.annotation.add_shape(
                            poly = checked_poly,
                            crs = [0,0],
                            structure = structure,
                            name = name
                        )

        # Closing progress bar
        if self.verbose:
            pbar.close()
            print(f'Found: {self.invalid_count} Invalid Polygons in annotations')

    def check_validity(self,poly):

        if not poly.is_valid:
            mod_shape = shapely.geometry.base.geom_factory(shapely.geos.lgeos.GEOSMakeValid(poly._geom))

            if mod_shape.geom_type == 'GeometryCollection' or mod_shape.geom_type == 'MultiPolygon':
                mod_shape = [i for i in list(mod_shape.geoms) if i.geom_type=='Polygon']
                if len(mod_shape)>1:
                    struct_areas = [i.area for i in mod_shape]
                    struct_poly = mod_shape[np.argmax(struct_areas)]
                elif len(mod_shape)==1:
                    struct_poly = mod_shape[0]
                else:
                    struct_poly = None

            elif mod_shape.geom_type == 'LineString':
                mod_shape = poly.buffer(0)
                mod_shape = shapely.geometry.base.geom_factory(shapely.geos.lgeos.GEOSMakeValid(mod_shape._geom))
                
                if mod_shape.geom_type == 'GeometryCollection' or mod_shape.geom_type == 'MultiPolygon':
                    mod_shape = [i for i in list(mod_shape.geoms) if i.geom_type=='Polygon']
                    if len(mod_shape)>1:
                        struct_areas = [i.area for i in mod_shape]
                        struct_poly = mod_shape[np.argmax(struct_areas)]
                    elif len(mod_shape)==1:
                        struct_poly = mod_shape[0]
                    else:
                        struct_poly = None
                else:
                    struct_poly = None

            else:
                struct_poly = None
        else:
            struct_poly = None

        if not struct_poly is None:
            self.invalid_count +=1

        return struct_poly
        
    def xml_save(self,filename):

        self.annotation.xml_save(filename,self.ann_dict)

    def geojson_save(self,filename):
        
        self.annotation.geojson_save(filename)

    def json_save(self,filename):

        self.annotation.json_save(filename)







