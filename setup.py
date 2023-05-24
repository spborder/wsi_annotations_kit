import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="wsi_annotations_kit",
    version="1.0.2",
    author="Sam Border",
    author_email="sam.border2256@gmail.com",
    description="Utility functions for generating, saving, and converting annotation files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/spborder/wsi_annotations_kit",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)