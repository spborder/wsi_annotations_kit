import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="wsi-annotations-kit",
    version="1.4.8",
    author="Sam Border",
    author_email="sam.border2256@gmail.com",
    description="Utility functions for generating, saving, and converting annotation files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/spborder/wsi_annotations_kit",
    install_requires=[
        "lxml==4.9.2",
        "geojson==3.0.1",
        "shapely==2.0.1",
        "tqdm==4.65.0",
        "numpy>=1.20.0",
        "uuid",
        "scikit-image==0.20.0"
    ],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)