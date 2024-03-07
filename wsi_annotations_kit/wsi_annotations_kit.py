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
from math import ceil, floor

import shapely
from shapely.geometry import box, Polygon, shape
from shapely.validation import make_valid
from shapely.ops import unary_union
from skimage.measure import label, find_contours
from skimage.segmentation import clear_border
import uuid

#import matplotlib.pyplot as plt

class Object:
    def __init__(self,
                 shape_structure,
                 location,
                 structure_name,
                 name,
                 properties = None):
    
        self.poly = shape_structure
        self.crs = location
        self.structure = structure_name
        self.name = name
        self.properties = properties

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

        if obj.properties is None:
            self.geojson['features'].append(
                Feature(geometry=new_poly,properties={'label':obj.name,'structure':obj.structure})
            )
        else:
            prop_dict = obj.properties
            prop_dict['label'] = obj.name
            prop_dict['structure'] = obj.structure
            self.geojson['features'].append(
                Feature(geometry = new_poly, properties=prop_dict)
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
                new_region = self.json_add_region(o)
                if not new_region is None:
                    structure_dict['elements'].append(self.json_add_region(o))
            if len(structure_dict)>0:
                self.json.append({'annotation':structure_dict})

        if verbose:
            pbar.close()

    def json_create(self):
        self.json = []

    def json_add_region(self,obj):

        # Shifting object into global crs 
        poly_coords = list(obj.poly.exterior.coords)
        scaled_poly_coords = [(int(i[0]+obj.crs[0]),int(i[1]+obj.crs[1])) for i in poly_coords]

        new_poly = Polygon(scaled_poly_coords)

        if len(scaled_poly_coords)>0:
            if obj.properties is None:
                new_struct_dict = {
                    'type':'polyline',
                    'points':[list(i)+[0] for i in list(new_poly.exterior.coords)],
                    'id':uuid.uuid4().hex[:24],
                    'closed':True,
                    'user': {
                        'name': obj.name,
                        'structure': obj.structure
                    }
                }
            else:
                prop_dict = obj.properties
                prop_dict['name'] = obj.name
                prop_dict['structure'] = obj.structure
                new_struct_dict = {
                    'type':'polyline',
                    'points':[list(i)+[0] for i in list(new_poly.exterior.coords)],
                    'id':uuid.uuid4().hex[:24],
                    'closed':True,
                    'user': prop_dict
                }

            return new_struct_dict
        else:
            return None


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
            return 'Empty annotation object'
        else:
            return f'Annotation object containing: {len(list(self.objects.keys()))}'

    def add_shape(self, poly, box_crs, structure = None, name = None, properties = None):
        
        # poly = shapely polygon or other geometry
        # box_crs = upper left corner for location that this object is in 
        # structure = structure name (which structure is this in general)
        # name = shape name (individual object name)
        if structure not in self.objects:
            self.objects[structure] = [Object(poly,box_crs,structure,name,properties)]
        else:
            self.objects[structure].append(Object(poly,box_crs,structure,name,properties))
    
    def add_names(self,names):
        
        # Initializing 
        for n in names:
            self.objects[n] = []
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
                        poly_list = [(int(i[1]),int(i[0])) for i in poly_list]
                    else:
                        # Find the largest one
                        contours_size = [np.shape(i)[0] for i in obj_contours]
                        poly_list = obj_contours[np.argmax(contours_size)].tolist()
                        poly_list = [(int(i[1]),int(i[0])) for i in poly_list]

                    if len(poly_list)>2:
                        # Making polygon from contours 
                        obj_poly = Polygon(poly_list)
                        self.add_shape(obj_poly,box_crs,structure_name)


class Patch:
    def __init__(self,
                 left,
                 top,
                 right,
                 bottom,
                 patch_index,
                 edge_direction = None):
        
        self.left = left
        self.top = top
        self.right = right
        self.bottom = bottom
        self.patch_index = patch_index
        self.edge_direction = edge_direction

        # Determining if this patch is on an edge
        if self.patch_index[-1]==0:
            self.edge_patch = True

        else:
            self.edge_patch = False

        # Storing complete and incomplete objects for each structure
        self.complete_objects = {}
        self.incomplete_objects = {}

        self.n_incomplete = 0
    def __str__(self):
        return f'Patch at: {self.left}, {self.top}, {self.right}, {self.bottom}'



class AnnotationPatches(Annotation):
    def __init__(self,
                 mpp = None,
                 min_size = None,
                 clear_edges = False):
        super().__init__(mpp, min_size)

        self.patch_list = []
        self.stride = None
        self.n_patch = None

        # Hierarchy:
        #   - World
        #   - Super-Community (Maybe find another name)
        #   - Community
        #   - Neighborhood
        #   - Patch

        # Whether or not to exclude annotations that intersect with the edges (useful if only complete annotations within a region are desired)
        self.clear_edges = clear_edges

    def add_patch_shape(self, poly, patch_obj, structure = None, name = None, properties = None):
        
        # For this, just pass a list of all polys in the patch to make it easier
        patch_obj.complete_structures[structure] = []
        patch_obj.incomplete_structures[structure] = []

        outer_box = box(
            minx = 0,
            miny = 0,
            maxx = int(patch_obj.right - patch_obj.left),
            maxy = int(patch_obj.bottom - patch_obj.top)
        )

        for p in poly:
            if not p.intersects(outer_box):
                patch_obj.complete_objects[structure].append(
                    Object(
                        p, [patch_obj.left, patch_obj.top], structure, name, properties
                    )
                )
            else:
                patch_obj.incomplete_objects[structure].append(
                    Object(
                        p, [patch_obj.left, patch_obj.top], structure, name, properties
                    )
                )

        patch_obj.n_incomplete += sum([len(patch_obj.incomplete_objects[i]) for i in list(patch_obj.incomplete_objects.keys())])

        self.patch_list.append(patch_obj)

    def add_patch_mask(self, mask, patch_obj, mask_type, structure = None):

        # Adding a mask patch, different procedures for whether something is an edge patch or an interior patch
        # Can have 'one-hot-labeled' or 'labeled-one-hot' to specify different objects within a structure that intersect
        if 'one-hot' in mask_type:
            # Expecting [height,width,classes] shape
            for cls in range(np.shape(mask)[-1]):

                if not structure is None:
                    if type(structure)==list:
                        structure_name = structure[cls]

                    elif type(structure)==dict:
                        structure_name = structure[str(cls)]

                    elif type(structure)==str:
                        structure_name = structure
                    else:
                        print('Invalid "structure" kwarg input')
                        raise ValueError
                
                if structure_name not in self.structure_names:
                    self.objects[structure_name] = []
                    self.structure_names.append(structure_name)
                class_mask = np.uint8(mask[:,:,cls].copy())

                # If individual objects in this structure are already labeled (objects which are touching would otherwise be merged)
                if 'labeled' in mask_type:
                    label_copy = class_mask.copy()
                    class_mask = np.uint8(label_copy>0)

                # Checking for whether or not to clear edges
                if patch_obj.edge_patch:

                    # Determining edge direction to clear
                    edge_mask = np.zeros_like(class_mask)
                    # Iterating through different edge directions
                    for edge in patch_obj.edge_direction:
                        if edge==0:
                            edge_mask[0,:] = 1
                        elif edge==1:
                            edge_mask[:,0] = 1
                        elif edge==2:
                            edge_mask[-1,:] = 1
                        elif edge==3:
                            edge_mask[:,-1] = 1

                    if self.clear_edges:
                        # Using edge mask to clear certain borders of the current mask
                        class_mask = clear_border(class_mask,mask = edge_mask)

                # Finding the structures which intersect with other edges (which will be added to "incomplete_structures")
                cleared_borders = clear_border(class_mask)

                # Edge structures:
                incomplete_structures = class_mask - cleared_borders
                # Interior, contained structures
                complete_structures = cleared_borders.copy()

                if not 'labeled' in mask_type:
                    labeled_incomplete, n_incomplete = label(incomplete_structures,background = 0, return_num = True)
                    labeled_complete, n_complete = label(complete_structures,background=0,return_num=True)
                else:
                    labeled_incomplete = incomplete_structures * label_copy
                    n_incomplete = len(np.unique(labeled_incomplete))-1
                    labeled_complete = complete_structures * label_copy
                    n_incomplete = len(np.unique(labeled_complete))-1


                # Adding structure list to patch object
                patch_obj.incomplete_objects[structure_name] = []
                patch_obj.complete_objects[structure_name] = []

                for i in np.unique(labeled_incomplete).tolist()[1:]:
                    # Finding object contours:
                    obj_contours = find_contours(labeled_incomplete==i)

                    if len(obj_contours)==1:
                        obj_contours = obj_contours[0].tolist()
                    else:
                        # Finding the largest one
                        contours_size = [np.shape(i)[0] for i in obj_contours]
                        obj_contours = obj_contours[np.argmax(contours_size)].tolist()

                    poly_list = [(int(i[0]+patch_obj.left),int(i[1]+patch_obj.top)) for i in obj_contours]

                    if len(poly_list)>2:
                        obj_poly = Polygon(poly_list)

                        patch_obj.incomplete_objects[structure_name].append(
                            Object(
                                obj_poly, [0,0], structure_name
                            )
                        )
                
                # Repeating for complete structures
                for i in np.unique(labeled_complete).tolist()[1:]:
                    obj_contours = find_contours(labeled_complete==i)

                    if len(obj_contours)==1:
                        obj_contours = obj_contours[0].tolist()
                    else:
                        # Finding largest one
                        contours_size = [np.shape(i)[0] for i in obj_contours]
                        obj_contours = obj_contours[np.argmax(contours_size)].tolist()

                    poly_list = [(int(i[0]+patch_obj.left),int(i[1]+patch_obj.top)) for i in obj_contours]

                    if len(poly_list)>2:
                        obj_poly = Polygon(poly_list)
                        patch_obj.complete_objects[structure_name].append(
                            Object(
                                obj_poly, [0,0], structure_name
                            )
                        )

        else:
            # Expecting mask format [height, width]
            if type(structure)==list:
                n_struct = len(structure)
                for s in structure:
                    if s not in self.structure_names:
                        self.objects[s] = []
                        self.structure_names.append(s)

            elif type(structure)==dict:
                n_struct = len(list(structure.keys()))
                for s in list(structure.keys()):
                    if s not in self.structure_names:
                        self.objects[s] = []
                        self.structure_names.append(s)

            elif type(structure)==str:
                n_struct =  1
                if structure not in self.structure_names:
                    self.objects[structure] = []
                    self.structure_names.append(structure)
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
                
                class_mask = np.uint8(mask.copy()==struct_idx)
                if len(np.unique(class_mask).tolist())>1:
                    # Checking for whether or not to clear edges
                    if patch_obj.edge_patch:

                        # Determining edge direction to clear
                        edge_mask = np.zeros_like(class_mask)
                        # Iterating through different edge directions
                        for edge in patch_obj.edge_direction:
                            if edge==0:
                                edge_mask[0,:] = 1
                            elif edge==1:
                                edge_mask[:,0] = 1
                            elif edge==2:
                                edge_mask[-1,:] = 1
                            elif edge==3:
                                edge_mask[:,-1] = 1

                        if self.clear_edges:
                            class_mask = np.uint8(clear_border(class_mask, mask = edge_mask>0))

                    # Finding the structures which intersect with other edges
                    cleared_borders = clear_border(class_mask)

                    # Edge structures:
                    incomplete_structures = class_mask - cleared_borders

                    # Interior, contained structures
                    complete_structures = cleared_borders.copy()

                    labeled_incomplete, n_incomplete = label(incomplete_structures, background = 0, return_num=True)
                    labeled_complete, n_complete = label(complete_structures, background = 0, return_num = True)

                    # Adding structure list to patch object
                    patch_obj.incomplete_objects[structure_name] = []
                    patch_obj.complete_objects[structure_name] = []

                    if len(np.unique(labeled_incomplete).tolist())>1:
                        
                        for i in np.unique(labeled_incomplete).tolist()[1:]:
                            obj_contours = find_contours(np.float32(labeled_incomplete==i))

                            for cont in obj_contours:
                                poly_list = [(int(i[1]+patch_obj.left),int(i[0]+patch_obj.top)) for i in cont]

                                if len(poly_list)>2:
                                    obj_poly = Polygon(poly_list)

                                    if obj_poly.area>1 and obj_poly.is_valid:
                                        patch_obj.incomplete_objects[structure_name].append(
                                            Object(
                                                obj_poly, [0,0], structure_name, None, None
                                            )
                                        )
                    

                    if len(np.unique(labeled_complete).tolist())>1:
                        # Repeating for complete structures
                        for i in np.unique(labeled_complete).tolist()[1:]:
                            obj_contours = find_contours(labeled_complete==i)

                            for cont in obj_contours:

                                poly_list = [(int(i[1]+patch_obj.left),int(i[0]+patch_obj.top)) for i in cont]

                                if len(poly_list)>2:
                                    obj_poly = Polygon(poly_list)
                                    if obj_poly.area>1 and obj_poly.is_valid:
                                        patch_obj.complete_objects[structure_name].append(
                                            Object(
                                                obj_poly, [0,0], structure_name, None, None
                                            )
                                        )

        # Adding number of incomplete objects for each structure to the patch object
        patch_obj.n_incomplete += sum([len(patch_obj.incomplete_objects[i]) for i in patch_obj.incomplete_objects])

    def find_adjacent(self,patch_index):
        # Patch index, mixed with self.n_neighbors returns all possible intersecting patch indices
        #print(f'base_patch index: {patch_index}')
        neighbor_patches = []
        for n_x in range(-(self.n_neighbors),(self.n_neighbors+1)):
            for n_y in range(-(self.n_neighbors),(self.n_neighbors+1)):
                if patch_index[0]+n_y>=0 and patch_index[0]+n_y<self.n_patch[0]:
                    if patch_index[1]+n_x>=0 and patch_index[1]+n_x<self.n_patch[1]:
                        if patch_index[0]+n_y==0 or patch_index[1]+n_x==0:
                            distance_from_edge = 0
                        else:
                            distance_from_edge = np.minimum((self.n_patch[0]-1-(patch_index[0]+n_y)),(self.n_patch[1]-1-(patch_index[1]+n_x)))

                        neighbor_patches.append([patch_index[0]+n_y, patch_index[1]+n_x, distance_from_edge])

        #print(f'neighboring patches: {neighbor_patches}')

        return neighbor_patches

    def clean_patches(self, merge_method = 'union'):
        
        # First adding all the complete objects for each patch
        pre_objects = {}
        for patch in self.patch_list:
            complete_structures = list(patch.complete_objects.keys())
            for s in complete_structures:
                #print(f'structure: {s} has: {len(patch.complete_objects[s])} complete structures')
                if s not in pre_objects:
                    pre_objects[s] = []
                    if s not in self.structure_names:
                        self.objects[s] = []
                        self.structure_names.append(s)
                
                pre_objects[s].extend([
                    Object(
                        p.poly.simplify(0.5,preserve_topology=False), [0,0], s, None, None
                    )
                    for p in patch.complete_objects[s]
                ])

        # Post-processing annotations, merging intersecting, incomplete annotations from adjacent patches
        all_patches_with_incomplete = [i for i in self.patch_list if i.n_incomplete>0]
        # It's possible some patches will overlap so that an "incomplete" object will be "complete" in another patch.
        # Therefore it should be okay to leave some trailing "incomplete" patches.
        intersecting_groups = []
        for base_patch in all_patches_with_incomplete:
            
            # Finding the neighbors for this patch
            possible_neighbors = self.find_adjacent(base_patch.patch_index)
            # This is a group of that incomplete patch and all its possible neighbors with incomplete objects
            intersecting_groups.append([i for i in self.patch_list if i.patch_index in possible_neighbors])
        
        # Possibly not the most efficient method
        inc_pre_objects = {}
        for g_idx,neighborhood in enumerate(intersecting_groups):
            # Getting the structures which have incomplete objects within each patch in a group
            structures_with_incomplete = []
            all_incomplete_structures = []
            for p in neighborhood:
                incomplete_structures = [i for i in p.incomplete_objects if len(p.incomplete_objects[i])>0]
                structures_with_incomplete.append(incomplete_structures)
                all_incomplete_structures.extend(incomplete_structures)
            
            unique_structures = np.unique(all_incomplete_structures)
            for u_idx,u in enumerate(unique_structures):
                if u not in inc_pre_objects:
                    inc_pre_objects[u] = []

                # Getting all the patches in the group which also have incomplete objects for that structure
                p_with_u = [neighborhood[i] for i in range(len(neighborhood)) if u in structures_with_incomplete[i]]

                # Now getting incomplete objects for that structure in this sub-group and finding the ones that intersect
                incomplete_objects_in_structure = []
                for pu_idx,pu in enumerate(p_with_u):
                    #print(f'group: {g_idx}/{len(intersecting_groups)}, structure: {u_idx}/{len(unique_structures)}, patch: {pu_idx}/{len(p_with_u)}')
                    for i in pu.incomplete_objects[u]:
                        if i.poly.buffer(0).is_valid:
                            # Not adding incomplete objects that intersect with complete objects
                            if not any([i.poly.buffer(0).intersects(j.poly) for j in pre_objects[u]]):
                                if any([i.poly.buffer(0).overlaps(j) for j in incomplete_objects_in_structure]):
                                    combination = unary_union([i.poly.buffer(0)]+incomplete_objects_in_structure)
                                    if combination.geom_type=='Polygon':
                                        incomplete_objects_in_structure.append(combination)
                                    elif combination.geom_type=='MultiPolygon':
                                        for piece in combination.geoms:
                                            incomplete_objects_in_structure.append(piece)
                                    elif combination.geom_type=='GeometryCollection':
                                        for piece in combination.geoms:
                                            if piece.geom_type=='Polygon':
                                                incomplete_objects_in_structure.append(piece)
                                else:
                                    incomplete_objects_in_structure.append(i.poly.buffer(0))

                        else:
                            valid_poly = make_valid(i.poly.buffer(0))
                            if not any([valid_poly.intersects(j.poly) for j in pre_objects[u]]):
                                if any([valid_poly.overlaps(j) for j in incomplete_objects_in_structure]):
                                    combination = unary_union([valid_poly]+incomplete_objects_in_structure)
                                    if combination.geom_type=='Polygon':
                                        incomplete_objects_in_structure.append(combination)
                                    elif combination.geom_type=='MultiPolygon':
                                        for piece in combination.geoms:
                                            incomplete_objects_in_structure.append(piece)
                                    elif combination.geom_type=='GeometryCollection':
                                        for piece in combination.geoms:
                                            if piece.geom_type=='Polygon':
                                                incomplete_objects_in_structure.append(piece)
                                else:
                                    incomplete_objects_in_structure.append(i.poly.buffer(0))

                inc_pre_objects[u].extend([
                    Object(
                        inc, [0,0], u, None, None
                    )
                    for inc in incomplete_objects_in_structure
                ])

        # Adding merged, non-overlapping incomplete structures
        for st in self.structure_names:
            if st in inc_pre_objects:
                pre_objects[st].extend([
                    Object(
                        i.poly,[0,0],st,None,None
                    )
                    for i in inc_pre_objects[st]
                ])

        # Making sure all structures in pre_objects are added:
        for structure in pre_objects:
            # Merging structures which overlap (test)
            if merge_method=='union':
                merged_objects = unary_union([i.poly.buffer(0) for i in pre_objects[structure] if i.poly.buffer(0).is_valid])
            elif merge_method=='envelope':
                merged_objects = shapely.MultiPolygon([i.poly.buffer(0) for i in pre_objects[structure] if i.poly.buffer(0).is_valid]).envelope
            elif merge_method=='minimum_rotated_rectangle':
                merged_objects = shapely.MultiPolygon([i.poly.buffer(0) for i in pre_objects[structure] if i.poly.buffer(0).is_valid]).minimum_rotated_rectangle
            elif merge_method=='convex_hull':
                merged_objects = shapely.MultiPolygon([i.poly.buffer(0) for i in pre_objects[structure] if i.poly.buffer(0).is_valid]).convex_hull

            else:
                print('Invalid merge_method')
                print('Choose from: union, envelope, minimum_rotated_rectangle, or convex_hull')
                raise ValueError
            
            # Adding objects finally
            if merged_objects.geom_type=='Polygon':
                #print(structure)
                #print(merged_objects.area)
                self.objects[structure].append(
                    Object(
                        merged_objects.simplify(0.5,preserve_topology=False), [0,0], structure, None, {'merged': True}
                    )
                )
            elif merged_objects.geom_type=='MultiPolygon':
                # Iterating through and adding each polygon
                for obj in merged_objects.geoms:
                    self.objects[structure].append(
                        Object(
                            obj.simplify(0.5,preserve_topology=False), [0,0], structure, None, {'merged': True}
                        )
                    )
            elif merged_objects.geom_type=='GeometryCollection':
                # only add the polygons
                for obj in merged_objects.geoms:
                    if obj.geom_type=='Polygon':
                        self.objects[structure].append(
                            Object(
                                obj.simplify(0.5,preserve_topology=False), [0,0], structure, None, {'merged': True}
                            )
                        )

    def define_patches(self, region_crs, height, width, patch_height, patch_width, overlap_pct):
        # Used to increase efficiency, pre-allocating patch coordinates and adjacency
        # region_crs = [left, top], x,y coordinates for upper left hand portion of the region
        self.region_crs = region_crs
        self.height = height
        self.width = width
        self.patch_height = patch_height
        self.patch_width = patch_width
        self.overlap_pct = overlap_pct

        if self.overlap_pct<0 or self.overlap_pct>=1:
            print('Overlap percentage must be less than 1 and greater than or equal to 0.')
            raise ValueError

        # This will tell you how many patches intersect with each other from each direction
        if self.overlap_pct>0:
            self.n_neighbors = floor(1/(1-self.overlap_pct))
        else:
            self.n_neighbors = 1

        if height <= patch_height and width <= patch_width:
            # In this case there will be a single patch (0,0) with edges in all directions
            self.patch_list = [
                Patch(
                    left = region_crs[0],
                    top = region_crs[1],
                    right = region_crs[0]+width,
                    bottom = region_crs[1]+height,
                    patch_index = [0,0],
                    edge_direction = [0,1,2,3]
                )
            ]
            self.n_patch = [1,1]
            self.stride = [width,height]
            self.n_neighbors = 0

        else:

            # Defining coordinates for each patch and adding to patch_indices array
            if width >= patch_width:
                stride_x = int(patch_width*(1-overlap_pct))
                n_patch_x = 1+floor((width-patch_width)/stride_x)

            else:
                stride_x = int(patch_width)
                n_patch_x = 1
            
            if height >= patch_height:
                stride_y = int(patch_height*(1-overlap_pct))
                n_patch_y = 1+floor((height-patch_height)/stride_y)

            else:
                stride_y = int(patch_height)
                n_patch_y = 1

            self.stride = [stride_x, stride_y]

            # Defining start coordinates for each direction. This is for same size patches.
            col_starts = [int(region_crs[0]+(i*stride_x)) for i in range(0,n_patch_x)]
            row_starts = [int(region_crs[1]+(i*stride_y)) for i in range(0,n_patch_y)]

            # The last patch for each direction. Again, this is the same size as the other patches so there may be higher overlap in the last patches.
            if width >= patch_width:
                col_starts.append(int(width-patch_width))
            
            if height >= patch_height:
                row_starts.append(int(height-patch_height))

            self.n_patch = [len(row_starts),len(col_starts)]

            # Creating patch_indices array with patch coordinates.
            for r_idx,r in enumerate(row_starts):
                for c_idx,c in enumerate(col_starts):
                    
                    # This parameter determines distance from an edge, if it is equal to zero, this is an edge patch.
                    if r_idx==0 or c_idx==0:
                        distance_from_edge = 0
                    else:
                        distance_from_edge = np.minimum((len(row_starts)-1-r_idx),(len(col_starts)-1-c_idx))

                    if not distance_from_edge==0:
                        self.patch_list.append(
                            Patch(
                                left = c,
                                top = r,
                                right = c + self.patch_width,
                                bottom = r + self.patch_height,
                                patch_index = [r_idx,c_idx,distance_from_edge]
                            )
                        )
                    else:
                        if r_idx==0 and c_idx==0:
                            # top-left corner
                            edge_direction = [0,1]
                        elif r_idx==0 and c_idx<len(col_starts)-1:
                            # top edge
                            edge_direction = [0]
                        elif r_idx==0 and c_idx==len(col_starts)-1:
                            # top-right corner
                            edge_direction = [0,3]
                        elif r_idx<len(row_starts)-1 and c_idx==0:
                            # left edge
                            edge_direction = [1]
                        elif r_idx==len(row_starts)-1 and c_idx==0:
                            # bottom-left corner
                            edge_direction = [1,2]
                        elif r_idx==len(row_starts)-1 and c_idx<len(col_starts)-1:
                            # bottom edge
                            edge_direction = [2]
                        elif r_idx==len(row_starts)-1 and c_idx==len(col_starts)-1:
                            # bottom-right corner
                            edge_direction = [2,3]
                        elif r_idx<len(row_starts)-1 and c_idx==len(col_starts)-1:
                            # right edge
                            edge_direction = [3]

                        self.patch_list.append(
                            Patch(
                                left = c,
                                top = r,
                                right = c+self.patch_width,
                                bottom = r+self.patch_height,
                                patch_index = [r_idx, c_idx, distance_from_edge],
                                edge_direction = edge_direction
                            )
                        )

    def write_patches(self):
        
        # Adding patches to current objects
        if 'Patches' not in self.structure_names:
            self.objects['Patches'] = []
            self.structure_names.append('Patches')
        
        for p in self.patch_list:

            self.objects['Patches'].append(
                Object(
                    shape_structure = box(p.left, p.top, p.right, p.bottom),
                    location = [0,0],
                    structure_name = 'Patches',
                    name = None,
                    properties = None
                )
            )

    def __iter__(self):

        if not len(self.patch_list)==0:
            self.patch_idx = 0

            return self
        else:
            print('Need to call define_patches first!')
            raise AttributeError

    def __next__(self):

        if self.patch_idx ==len(self.patch_list):
            raise StopIteration
        else:
            new_patch = self.patch_list[self.patch_idx]
            self.patch_idx += 1

            return new_patch


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

        annotation = Annotation()
        annotation.add_names(list(self.ann_dict.keys()))

        if self.file_ext=='xml':
            
            # Verboseness just for number of layers for xmls
            if self.verbose:
                pbar = tqdm(list(self.ann_dict.keys()))

            # Add MPP here if there is one

            tree = ET.parse(self.starting_file)
            for idx, structure in enumerate(self.ann_dict):

                structures_in_xml = tree.getroot().findall(f'Annotation[@Id="{str(self.ann_dict[structure])}"]/Regions/Region')

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
                        checked_poly, add_to_annotation = self.check_validity(region_poly)

                        # Adding to current annotations if valid
                        if add_to_annotation:
                            annotation.add_shape(
                                poly = checked_poly,
                                box_crs = [0,0],
                                structure = structure,
                                name = name)
                        else:
                            print('Found invalid shape')

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
                elif 'classification' in f['properties']:
                    if 'name' in f['properties']['classification']:
                        structure_name = f['properties']['classification']['name']
                    else:
                        structure_name = f'Structure'
                else:
                    structure_name = 'Structure'
                
                if 'name' in f['properties']:
                    name = f['properties']['name']
                else:
                    name = None
                
                annotation.add_shape(
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
                if 'name' in st:
                    structure = st['name']
                elif 'properties' in st:
                    if 'classification' in st['properties']:
                        if 'name' in st['properties']['classification']:
                            structure = st['properties']['classification']['name']
                        else:
                            structure = f'Structure_{st_idx}'
                    else:
                        structure = f'Structure_{st_idx}'
                else:
                    structure = f'Structure_{st_idx}'

                if self.verbose:
                    pbar.update(st_idx)
                    if 'elements' in st:
                        pbar.set_description(f'Working on: {structure}, found: {len(st["elements"])}')
                    
                    
                if 'elements' in st:
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
                            annotation.add_shape(
                                poly = checked_poly,
                                box_crs = [0,0],
                                structure = structure,
                                name = name
                            )
                elif 'geometry' in st:
                    coordinates = st['geometry']['coordinates']
                    coordinates = np.squeeze(coordinates)
                    coordinates = [(i[0],i[1]) for i in coordinates]

                    if 'name' in st['properties']:
                        name = st['properties']['name']
                    else:
                        name = None

                    region_poly = Polygon(coordinates)
                    checked_poly, add_to_annotation = self.check_validity(region_poly)

                    if add_to_annotation:
                        annotation.add_shape(
                            poly = checked_poly,
                            box_crs = [0,0],
                            structure = structure,
                            name = name
                        )

        # Closing progress bar
        if self.verbose:
            pbar.close()
            print(f'Found: {self.invalid_count} Invalid Polygons in annotations')

        return annotation

    def check_validity(self,poly):

        if not poly.is_valid:
            mod_shape = make_valid(poly)

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
                mod_shape = make_valid(mod_shape)
                
                if mod_shape.geom_type == 'GeometryCollection' or mod_shape.geom_type == 'MultiPolygon':
                    mod_shape = [i for i in list(mod_shape.geoms) if i.geom_type=='Polygon']
                    if len(mod_shape)>1:
                        struct_areas = [i.area for i in mod_shape]
                        struct_poly = mod_shape[np.argmax(struct_areas)]
                    elif len(mod_shape)==1:
                        struct_poly = mod_shape[0]
                    else:
                        struct_poly = None

                elif mod_shape.geom_type=='Polygon':
                    struct_poly = mod_shape
                else:
                    struct_poly = None

            elif mod_shape.geom_type=='Polygon':
                struct_poly = mod_shape
            else:
                struct_poly = None
        else:
            struct_poly = poly

        if struct_poly is None:
            self.invalid_count +=1
            add_to_annotation = False
        else:
            add_to_annotation = True

        return struct_poly, add_to_annotation
        
    def xml_save(self,filename):

        self.annotation.xml_save(filename,self.ann_dict)

    def geojson_save(self,filename):
        
        self.annotation.geojson_save(filename)

    def json_save(self,filename):

        self.annotation.json_save(filename)







