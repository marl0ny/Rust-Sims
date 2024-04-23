use crate::complex::*;


// TODO: somehow make the array generic
// instead of limiting this to just Complex<T>
pub fn square_transpose_in_place<T: Copy>(
    array: &mut [Complex<T>], n: usize) {
    for i in 0..n {
        for j in i+1..n {
            let tmp = array[i*n + j];
            array[i*n + j] = array[j*n + i];
            array[j*n + i] = tmp;
        }
    }
}
