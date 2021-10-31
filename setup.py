from setuptools import setup

setup(
    name='lsdo_mesh',
    packages=[
        'lsdo_mesh',
    ],
    install_requires=[
        'numpy',
        'gmsh',
    ],
    version_config=True,
)