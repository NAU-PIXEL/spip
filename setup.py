import setuptools

with open('README.md', 'r') as f:
    long_description = f.read()
with open('requirements.txt', 'r') as f:
    requirements = f.read().strip('\n').split('\n')

package_data = {
    '': ['data/*'],
    }

setuptools.setup(
    name='spip',
    version='1.1.3',
    author='Aurélien Stcherbinine',
    author_email='aurelien.stcherbinine@nau.edu',
    description='Spacecraft Pixel footprint Projection',
    long_description=long_description,
    long_description_content_type='text/markdown',
    project_urls={
        'Source' : 'https://github.com/NAU-PIXEL/spip',
    },
    packages=setuptools.find_packages(),
    package_data=package_data,
    python_requires='>=3.6',
    setup_requires=['wheel'],
    install_requires=requirements,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy'
        ]
    )
