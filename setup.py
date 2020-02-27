# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open("README.md", "r") as f:
    readme = f.read()

setup(
    name="liputils",
    version="0.12",
    description="A small Python package to manipulate complex lipids.",
    long_description=readme,
    long_description_content_type="text/markdown",
    author="Stefano Manzini",
    author_email="stefano.manzini@gmail.com",
    url="https://github.com/Stemanz/liputils",
    download_url="https://github.com/Stemanz/liputils/archive/master.zip",
    license="GPL-3.0",
    packages=find_packages(exclude=("tests", "docs", "images")),
    keywords = ['lipids', 'lipids residues', 'lipidomics'],
    install_requires=[
          'numpy',
          'pandas',
      ],
)
