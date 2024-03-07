import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="wsi-annotations-kit",
    version="1.4.1",
    author="Sam Border",
    author_email="sam.border2256@gmail.com",
    description="Utility functions for generating, saving, and converting annotation files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/spborder/wsi_annotations_kit",
    install_requires=[
        "lxml",
        "geojson",
        "shapely==2.0.1",
        "tqdm",
        "numpy",
        "uuid",
        "scikit-image"
    ],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)