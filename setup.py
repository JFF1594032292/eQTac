import setuptools

with open("README.md", "r") as f1:
    long_description = f1.read()

setuptools.setup(
    name="eQTac",
    version="1.0.9",
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
    python_requires='>=3.7',
    install_requires=[
        'numpy >= 1.21.6',
        'pandas >= 1.2.3',
        'pybedtools >= 0.8.1',
        'pysam >= 0.15.3',
        'rpy2 >= 3.5.11',
        'scipy >= 1.7.3'],
)
