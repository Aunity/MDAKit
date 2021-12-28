from setuptools import setup, find_packages

packages = [
    'MDAKit',
    'MDAKit.lib',
    'MDAKit.tools',
    'MDAKit.utils',
]
install_requires = [
        'tqdm>=4.62.2',
        'numpy>=1.21.2',
        'pandas>=1.3.2',
        'MDAnalysis>=2.0.0',
        'mdtraj>=1.9.6',
]
setup(
    name = 'MDAKit',
    version='1.0.0',
    author='Maohua Yang',
    author_email='maohuay@hotmail.com',
    description=('A kit with lots of tools for Molecular dynamic simulation analysis.'),
    license=None,
    keywords='molecular dynamic, analysis',
    url='https://github.com/Aunity/MDAKit',
    packages=find_packages(),
    #packages=packages,
    entry_points={'console_scripts': [
         'mdakit = MDAKit.main:main',
     ]},

)
