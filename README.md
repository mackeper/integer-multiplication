# Integer Multiplication in n log n
Marcus Östling

## Project plan
The purpose of this project is to a implement the algorithm presented in the
paper “Integer multiplication in time O(n log n)
” by David Harvey and Joris
van der Hoeven, and then compare it to a simple integer multiplication using
a FFT in one dimension.
###  Decide on parameters (week 51)
The new algorithm have a couple on parameters that I have to decide on.
These parameters will be chosen so that it is possible to run the algorithms on
8 GB of ram. Input sizes, output size, intermediate results and some space for
the OS as well.

###  Implement simple FFT multiplication (week 52)
Implement a simple FFT integer multiplication, imiliar to Schönhage-Strassen,.
Given two integers, split them up into coefficients for a polynomial, use FFT,
apply pointwise multiplication, inverse FFT back into one integer.

#### Progress
The planned input size of 2^30 was harder to achieve than predicted. This algorithm
can handle 2^25 right now, coefficients can only be around 16 bits (in 8 byte double)
and the coefficients is complex values so they need 2x8 bytes (two doubles).

###  Implement new algorithm (week 1-3)
This algorithm is inspired by Schönhage-Strassen: Given two integers, split
them up into coefficients for a polynomial, use the algorithm given by the
paper, convert back into one integer.
###  Run algorithms and compare (week 4)
Run both algorithms on the same inputs, up to large integer (size will be
decided during week 51), and measure their run-times.
### Write report (week 5-6)
All results should be done and I will write a complete report for the project.
