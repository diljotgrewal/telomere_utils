import versioneer
from setuptools import setup, find_packages

setup(
    name='telomere_utils',
    packages=find_packages(),
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='python utilities',
    author='Diljot Grewal',
    author_email='diljot.grewal@gmail.com',
    entry_points={
        'console_scripts': [
            'telomere_utils = telomere_utils.main:main',
        ]
    },
    package_data={'': ['*.py', '*.R', '*.npz', "*.yaml", "data/*"]}
)
