"""

Testing file conversion


"""

from wsi_annotations_kit import wsi_annotations_kit as wak

test_file = 'examples/test_annotation.xml'
ann_dict = {'Glomeruli':1}

test_converter = wak.Converter(test_file,ann_dict)
test_converter.json_save(test_file.replace('.xml','.json'))




