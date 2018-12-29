nmrPipe -in ./$1/test.fid \
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5 \
| nmrPipe  -fn ZF -zf 3 \
| nmrPipe  -fn FT -verb \
| nmrPipe  -fn PS -p0 0. -p1 0.0 -di \
  -verb -ov -out ./$1/test.ft
