MOD_NAME=vwake

python setup.py build_ext --inplace
swig -python ${MOD_NAME}.i

echo "Install Complete"

python pycallc.py

