"""

Testing file conversion


"""

from wsi_annotations_kit import wsi_annotations_kit as wak

test_file = 'C:\\Users\\Sam\\Downloads\\XY04_IU-21-020F.xml'
ann_dict = {'Cortical interstitium':1,'Medullary interstitium':2,'Glomeruli':3,'Sclerotic Glomeruli':4,'Tubules':5,'Arteries and Arterioles':6}

test_converter = wak.Converter(test_file,ann_dict)
test_converter.json_save(test_file.replace('.xml','.json'))




