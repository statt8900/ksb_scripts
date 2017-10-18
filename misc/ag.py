from ase.visualize import view
from ase.io import read
import sys

view(read(sys.argv[1]))