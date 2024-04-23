from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = ('A computational framework for associating cancer cell-intrinsic driver genes with'
               ' -extrinsic cell-cell communication')
LONG_DESCRIPTION = ('For systematic investigation of the relationship between cancer cell-intrinsic and -extrinsic factors'
                    ', we develop a general computational framework, Driver2Comm, for associating driver genes of '
                    'cancer cells with CCC in the TME using single-cell transcriptomics data. ')
setup(
    name="Driver2Comm",
    version=VERSION,
    keywords=["scRNA-seq", "cell-cell communication", "cancer driver gene", "tumor microenvironment"],
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    url="https://github.com/xierunzhi",
    author="Runzhi Xie, Yuxuan Hu*, Gao Lin*",
    author_email="rzxie1998@gmail.com",
    maintainer="Runzhi Xie",
    maintainer_email="rzxie1998@gmail.com",
    packages=find_packages(),
    include_package_data=True,
    platforms="any",
    license="MIT Licence",
    install_requires=[
        "numpy",
        "pandas",
        "statsmodels",
        "networkx",
        "matplotlib",
        "scipy"
    ]
)
