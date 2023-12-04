from setuptools import setup, find_packages


with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='fluvp',
    version='1.0.0',
    author='viralInformatics',
    author_email='ssdx202203@163.com',
    description='A tool for predicting influenza virus virulence level from sequence data.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/viralInformatics/FluVirulencePredictor',
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'fluvp=fluvp.fluvp:main',
        ],
    },
    package_data={
        '': ['data/*.*', 'model/*.*','tests/*.*'],
    },
    install_requires=required,
    zip_safe=False,
 
    classifiers=[
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    keywords='influenza virulence prediction bioinformatics',
)
