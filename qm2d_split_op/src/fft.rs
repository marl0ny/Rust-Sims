use crate::complex::*;

/* Reverse bit sort an array, where the size of the array
must be a power of two.
*/
fn reverse_bit_sort<T: Copy>(array: &mut [Complex<T>], n: usize) {
    let mut u: usize;
    let mut d: usize;
    let mut rev: usize;
    for i in 0..n {
        u = 1;
        d = n >> 1;
        rev = 0;
        while u < n {
            rev += d*((i&u)/u);
            u <<= 1;
            d >>= 1;
        }
        if rev >= i {
            let tmp = array[i];
            array[i] = array[rev];
            array[rev] = tmp;
        } 
    }
}

/* This function implements the iterative in place radix-2 
Cooley-Tukey Fast Fourier Transform Algorithm. The size of the input
array must be a power of two, or else bad things will happen. There
are currently no checks done to ensure this.

References:

Wikipedia - Cooley–Tukey FFT algorithm
https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm

MathWorld Wolfram - Fast Fourier Transform:
http://mathworld.wolfram.com/FastFourierTransform.html 

William Press et al.
12.2 Fast Fourier Transform (FFT) - Numerical Recipes
https://websites.pmc.ucsc.edu/~fnimmo/eart290c_17/NumericalRecipesinF77.pdf

*/
pub fn base_f32_fft_in_place(array: &mut [Complex<f32>], 
                             size: usize, is_inverse: bool) {
    reverse_bit_sort(array, size);
    let mut block_size: usize = 2;
    let sgn: f64 = if is_inverse {-1.0} else {1.0};
    while block_size <= size {
        let e1: Complex<f64> = Complex {
            real: f64::cos(2.0*std::f64::consts::PI
                           *1.0/(block_size as f64)),
            imag: sgn*f64::sin(2.0*std::f64::consts::PI
                               *1.0/(block_size as f64)),
        };
        let mut j: usize = 0;
        while j < size {
            let mut e: Complex<f64> = Complex {real: 1.0, imag: 0.0};
            for i in 0..block_size/2 {
                let even: Complex<f64> = array[j + i].into();
                let odd: Complex<f64> = array[j + i + block_size/2].into();
                let s: f64 = if is_inverse && block_size == size 
                    {1.0/(size as f64)} else {1.0};
                array[j + i] = (s*(even + odd*e)).into();
                array[j + i + block_size/2] = (s*(even - odd*e)).into();
                e = e*e1;
            }
            j += block_size;
        }
        block_size *= 2;
    }
}

pub fn fft_in_place(array: &mut [Complex<f32>], size: usize) {
    base_f32_fft_in_place(array, size, false);
}

pub fn ifft_in_place(array: &mut [Complex<f32>], size: usize) {
    base_f32_fft_in_place(array, size, true);
}

/* Perform the fft algorithm on each row of an array.
Rows are placed into separate groups, where each group is
handled by its own thread.

Multithreading reference:
https://doc.rust-lang.org/book/ch16-01-threads.html
https://doc.rust-lang.org/book/ch16-02-message-passing.html
*/
pub fn horizontal_square_fft(is_inverse: bool,
                             array: &mut [Complex<f32>],
                             width: usize, // Width of the square array
                             th_total: usize // Number of threads to use
                            ) {
    /*
    // https://doc.rust-lang.org/std/thread/fn.scope.html 
    std::thread::scope(|s| {
        for th_index in 0..th_total {
            let elem_count: usize = (width*width)/th_total;
            let slice: &mut [Complex<f32>] 
                = &mut array[th_index*elem_count..(th_index+1)*elem_count];
            s.spawn(|| {
                for i in 0..width/th_total {
                    if is_inverse {
                        ifft_in_place(slice, width);
                    } else {
                        fft_in_place(slice, width);
                    }
                }
            });
        }
    });*/
    let mut receivers = std::vec::Vec::<
        std::sync::mpsc::Receiver<std::vec::Vec<Complex<f32>>>
        >::with_capacity(th_total);
    for th_index in 0..th_total {
        let (tx, rx) = std::sync::mpsc::channel();
        let mut v
            = std::vec::Vec::<Complex<f32>>::
            with_capacity(width*width/th_total);
        for i in th_index*width*width/th_total
            ..(th_index + 1)*width*width/th_total {
            v.push(array[i]);
        }
        std::thread::spawn(move || {
            for i in 0..width/th_total {
                if is_inverse {
                    ifft_in_place(&mut v.as_mut_slice()[i*width..(i+1)*width],
                                  width);
                } else {
                    fft_in_place(&mut v.as_mut_slice()[i*width..(i+1)*width],
                                 width);
                }
            }
            tx.send(v).unwrap();
        });
        receivers.push(rx);
    }
    let mut th_index: usize = th_total - 1;
    while let Some(r) = receivers.pop() {
        let v = r.recv().unwrap();
        for i in th_index*width/th_total..(th_index + 1)*width/th_total {
            let i_get = i - th_index*width/th_total;
            for j in 0..width {
                array[width*i + j] = v[i_get*width + j];
            }
        }
        th_index = if th_index == 0 {th_index} else {th_index-1};
    }
}
