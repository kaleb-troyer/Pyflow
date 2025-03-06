from setuptools import find_packages, setup

with open("pyflow/readme.md", "r") as f:
    long_description = f.read()

setup(
    name="htf-pyflow",
    version="0.0.10",
    description="",
    package_dir={"": "pyflow"},
    packages=find_packages(where="pyflow"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    author="KalebTroyer",
    author_email="kaleb.troyer@gmail.com",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    install_requires=[""],
    extras_require={
        "dev": ["twine>=5.1.1"],
    },
    python_requires=">=3.12",
)
