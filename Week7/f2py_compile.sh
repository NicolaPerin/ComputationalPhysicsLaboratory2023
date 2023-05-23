python3 -m numpy.f2py isingmetro.f90 -m isingf -h ising.pyf
python3 -m numpy.f2py -c ising.pyf isingmetro.f90
