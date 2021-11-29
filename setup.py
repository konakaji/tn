import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tn",
    version="0.0.1",
    author="kouhei nakaji",
    author_email="kohei.nakaji@keio.jp",
    description="Tool for haar tensor network computation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/konakaji/tn",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "matplotlib==3.5.0",
        "networkx==2.6.3"
    ],
    python_requires='>=3.7',
)
