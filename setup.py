import setuptools
import os

setuptools.setup(

    name="pathogenseq_scripts",
    version="0.1.3",
    packages=["pathogenseq_scripts"],
    license="MIT",
    long_description="Utilities to get from fastq files to a variant matrix",
    scripts= ["scripts/%s" % x for x in os.listdir("scripts")],

)
