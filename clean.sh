#!/bin/bash
python -m pip uninstall HybPSF
python3 -m pip uninstall HybPSF
python3.10 -m pip uninstall HybPSF
rm -r build
rm -r __pychahe__
rm -r HybPSF.egg-info build
make clean -C HybPSF/staridf
rm HybPSF/staridf/*.o
make clean -C HybPSF/source
rm HybPSF/source/*.o
make clean -C HybPSF/cfitsio
rm HybPSF/cfitsio/*.o
rm -r HybPSF/include
rm -r HybPSF/lib
rm -r HybPSF/__pychahe__
rm HybPSF/*.so