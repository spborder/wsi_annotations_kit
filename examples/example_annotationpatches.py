"""

Testing AnnotationPatches functionality for generating patches and annotations from those patches


"""

import os
import sys
import numpy as np

from PIL import Image

sys.path.append('..')
import wsi_annotations_kit.wsi_annotations_kit.wsi_annotations_kit as wak


def preprocess_mask(mask_image,mask_type):

    gray_image = 255-np.uint8(np.mean(mask_image,axis=-1))
    classes = np.unique(gray_image).tolist()
    classes = [i for i in classes if i != 0]
    print(classes)
    n_class = len(classes)

    if mask_type=='one-hot':
        # Processing for one-hot
        processed_mask = np.zeros((np.shape(gray_image)[0],np.shape(gray_image)[1],n_class))
        for class_idx,cls in enumerate(classes):
            processed_mask[:,:,class_idx] = (gray_image==cls)

    else:
        # Processing for non-one-hot, categorical labels
        processed_mask = np.zeros_like(gray_image)
        for class_idx,cls in enumerate(classes):
            processed_mask[np.where(gray_image==cls)] = class_idx+1

    return processed_mask



def main():

    # Loading example object mask
    test_image_path = './examples/test_image.png'
    test_image_mask = np.array(Image.open(test_image_path))
    test_mask_type = 'non-one-hot'

    # Preprocessing mask to be in either one-hot or class label format
    processed_mask = preprocess_mask(test_image_mask,test_mask_type)
    print(f'Number of objects: {len(np.unique(processed_mask))-1}')
    n_struct = len(np.unique(processed_mask).tolist())

    # Initializing annotation object (single "patch")
    annotation = wak.Annotation()
    annotation.add_names([f'Structure{i}' for i in range(n_struct-1)])
    annotation.add_mask(
        mask = processed_mask,
        box_crs = [0,0],
        mask_type=test_mask_type,
        structure = [f'Structure{i}' for i in range(n_struct-1)]
    )

    # Testing saving annotations
    annotation.xml_save('./examples/test_xml.xml')
    annotation.geojson_save('./examples/test_geojson.geojson')
    annotation.json_save('./examples/test_json.json')

    # Now doing the same thing but using the AnnotationPatches method
    patch_annotation = wak.AnnotationPatches()
    patch_annotation.define_patches(
        region_crs = [0,0],
        height = np.shape(processed_mask)[0],
        width = np.shape(processed_mask)[1],
        patch_height = 256,
        patch_width = 256,
        overlap_pct = 0.25
    )

    patch_annotation = iter(patch_annotation)
    keep_iterating = True
    while keep_iterating:
        try:
            new_patch = next(patch_annotation)
            mask_region = processed_mask[new_patch.top:new_patch.bottom, new_patch.left:new_patch.right]

            patch_annotation.add_patch_mask(
                mask = mask_region,
                patch_obj = new_patch,
                mask_type = test_mask_type,
                structure = [f'Structure{i}' for i in range(n_struct-1)]
            )

        except StopIteration:
            keep_iterating = False

    patch_annotation.clean_patches()
    # Testing saving annotations
    patch_annotation.xml_save('./examples/patch_test_xml.xml')
    patch_annotation.geojson_save('./examples/patch_test_geojson.geojson')
    patch_annotation.json_save('./examples/patch_test_json.json')






if __name__=='__main__':
    main()
































