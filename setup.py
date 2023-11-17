from setuptools import setup, find_packages
setup(
    name='hicool',
    version='0.1.0',
    author='wzl',
    author_email='401244520@qq.com',  # Replace with your email
    description='Analysis and Visualization package for Massive Hi-C data.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/401244520/hicool',  # Replace with the URL of your package/repository
    packages=find_packages(),  # Automatically find all packages and subpackages
    # install_requires=[
        # cooltools
    #     # Add your package dependencies here
    #     # 'numpy',
    #     # 'pandas',
    # ],
    # classifiers=[
    #     # Trove classifiers
    #     # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    #     'Development Status :: 3 - Alpha',
    #     'Intended Audience :: Developers',
    #     'Topic :: Software Development :: Build Tools',
    #     'License :: OSI Approved :: MIT License',  # Again, pick a license
    #     'Programming Language :: Python :: 3',  # Specify which pyhton versions that you want to support
    #     'Programming Language :: Python :: 3.10',
    # ],
    # python_requires='>=3.10',  # Minimum version requirement of the package
    # extras_require={
    #     'dev': [
    #         'check-manifest',
    #         'twine',
    #         # 'pytest>=3.7',  # Optional dependencies for development
    #     ],
    # },
    # entry_points={
    #     'console_scripts': [
    #         'hicool-cli=hicool.cli:main',  # Replace with your package's cli entry point if any
    #     ],
    # },
)