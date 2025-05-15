from setuptools import Extension, setup
from os import sep, getenv
from glob import glob

from Cython.Build import cythonize

# Define Cython extensions to build. The pyx files are assumed to be in src/ctfr/implementations.
IMPLEMENTATIONS_SOURCE_DIR = f"src{sep}ctfr_bm{sep}implementations"
IMPLEMENTATIONS_MODULE = "ctfr_bm.implementations"

def get_cy_extensions():
    method_cy_source_paths = glob(f"{IMPLEMENTATIONS_SOURCE_DIR}{sep}*.pyx")
    return [Extension(f"{IMPLEMENTATIONS_MODULE}.{path.split(sep)[-1].split('.')[0]}", [path]) for path in method_cy_source_paths]
    
extensions = get_cy_extensions()

annotate = bool(int(getenv("ANNOTATE", 0)))

print("Building from .pyx files.")
from Cython.Build import cythonize
compiler_directives = {"language_level": 3}
ext_modules = cythonize(
    extensions,
    annotate=annotate,
    compiler_directives=compiler_directives
)

setup(
    name="ctfr_bm",
    ext_modules=cythonize(ext_modules)
)
