# magnetostatics
Magnetostatics using brute-force application of the Biot-Savart law to wire elements. 

### Compatibility notes
This notebook was originally written in Julia 0.6.x. There have been some changes in both syntax and definition of base functions. The following adaptions have been made to this notebook to make it compatible with newer versions (Julia 1.x).

* `transpose()` $\rightarrow$ `permutedims()` for data manipulation, same goes for `adjoint()`, which can also be denoted by a single quotatio mark (').
* `sum(x, 1)` $\rightarrow$ `sum(x, dims = 1)` for summing over a specific dimension of an array `x`
* `1.^2 + z` $\rightarrow$ `1^2 .+ z` if `z` is an array of sorts
* `linspace(start, stop, n)` is deprecated, use `range(start, stop, length = n)` instead
* `zeros()` used to take any array as its argument, and return a zeros array with the same dimensions. No longer! You have to specify the dimensions. So instead of `zeros(A)` do `zeros(size(A))` if you want a zero array with the same size as `A`.
* It's not allowed anymore to set a range in an array equal to a scalar as in `A[start : end] = b`. Instead use `A[start : end] = b * ones(length( A[start : end] ))` (there's probably a neater way to do this, but it works).
* Creating an unintialized array no longer works with `Array{T, N}(dims)`, but you need to use `Array{T, N}(undef, dims)` instead.

<font color='red'>*Caution, with these changes this calculation is no longer compatible with Julia 0.6.x.*</font>