# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open("README.md") as f:
    readme = f.read()

with open("LICENSE") as f:
    license = f.read()

setup(
    name="liputils",
    packages = ["liputils"],
    version="0.1",
    description="Picks individual fatty acids from individual complex lipids",
    long_description=readme,
    author="Stefano Manzini",
    author_email="stefano.manzini@gmail.com",
    url="https://github.com/Stemanz/liputils",
    download_url="https://github.com/Stemanz/liputils/archive/master.zip",
    license=license,
    packages=find_packages(exclude=("tests", "docs", "images")),
    keywords = ['lipids', 'lipids residues', 'lipidomics'],
    install_requires=[
          'numpy',
          'pandas',
      ],
)
