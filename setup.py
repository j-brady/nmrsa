import sys
from distutils.core import setup

# Update file headers with running version of python
files = ["scripts/makeYamlFiles.py", "scripts/proc_with_yaml.py", "scripts/spec.py"]
for i in files:

    with open(i, "r") as f:
        lines = f.readlines()
        if lines[0].startswith("#!"):
            lines.pop(0)
            f_str = "".join(lines)
        else:
            f_str = "".join(lines)

        s = "#!" + sys.executable + "\n" + f_str
    with open(i, "w") as f:
        f.write(s)


setup(
    name="nmrsa",
    version="",
    description="NMR processing scripts",
    author="Jacob Brady and Rui Huang",
    packages=["nmrsa"],
    scripts=files,
    package_data={"nmrsa": ["templates/*"]},
)
