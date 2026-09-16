from setuptools import setup, find_packages
version = {}
with open("bactipipe/__version__.py") as f:
    exec(f.read(), version)

setup(
    name="bactipipe",
    version=version["__version__"],
    packages=find_packages(),  # Finds all packages under bactipipe/
    description="Bacterial WGS data analysis pipeline",
    author="Maurice Byukusenge",
    author_email="bmaurice@psu.edu",
    install_requires=[
        "Pillow>=10.0",
        "boto3>=1.40.6",
        "biopython>=1.85",
        "colorama",
        "openpyxl>=3.1",
        "pandas>=2.0",
        "psutil>=7.0",
        "pyfastx",
        "reportlab>=4.4",
        "SpeciesFinder==4.0.4",
        "tabulate>=0.9",
        "textwrap3>=0.9",
        "tqdm>=4.67",
    ],
    entry_points={
        "console_scripts": [
            "bactipipe=bactipipe.cli:main"
        ]
    },
    include_package_data=True,
    package_data={
        "bactipipe": ["data/lambda.fasta", "data/pathogenic_bacteria.txt", "data/*.tsv"],
    },
    python_requires=">=3.10",
    url="https://github.com/mauricebyuka/bactiPipe",  # optional
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
