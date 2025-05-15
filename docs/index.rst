Implementations code
====================

The code for the implementations of the CTFR algorithms (both baseline and ``ctfr``) is available in the following links:

.. include:: code_tree.rst

Note that the implementations above don't include the code for ``ctfr.ctfr_from_specs``, which wraps every ``ctfr`` algorithm. It prepares the input data, calls the corresponding method and normalizes the output representation (if applicable). Its code is available in the following link:

.. toctree:: 
   :caption: Source code (core)
   
   code/core

.. toctree::
   :maxdepth: 1
   :caption: See also

   Github (benchmarks) <https://github.com/b-boechat/ctfr_benchmarks>
   Github (ctfr) <https://github.com/b-boechat/ctfr>
   Documentation (ctfr) <https://ctfr.readthedocs.io/en/latest/>