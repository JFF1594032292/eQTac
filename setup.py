import setuptools

with open("README.md", "r") as f1:
    long_description = f1.read()

setuptools.setup(
    name="eQTac",
    version="1.0.17",
    author="Jiang Feng",
    author_email="1594032292@qq.com",
    description="The eQTac method.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JFF1594032292/eQTac",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
    ],
    python_requires='==3.8',
    install_requires=[
        'numpy == 1.22.4',
        'pandas == 1.4.3',
        'pybedtools == 0.8.2',
        'pysam == 0.16.0.1',
        'rpy2 == 3.5.11',
        'scipy == 1.8.1'],
)
