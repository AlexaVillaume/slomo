from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="slomo",
    author="Asher Wasserman",
    author_email="asher.d.wasserman@gmail.com",
    description="Jeans models for slow-rotator galaxies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/adwasser/Slomo.jl",
    install_requires=["julia"],
    setup_requires=["setuptools_scm"],
    use_scm_version=True,
    packages=find_packages(),
    include_package_data=True,
)
