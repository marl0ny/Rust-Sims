
use qm2d_split_op::constants::*;
use qm2d_split_op::transpose::*;
use qm2d_split_op::fft::*;
use qm2d_split_op::complex::*;
use qm2d_split_op::wavepacket::*;
use qm2d_split_op::bitmap::*;
use qm2d_split_op::domain_coloring::*;
use qm2d_split_op::serialize::*;

use std::env;


/* Initialize the V(x, y) term of the Schrodinger equation. */
fn init_potential(potential: &mut [Complex<f32>]) {
    let mut potential_low_res = std::vec::Vec::<u8>::with_capacity(32*16);
    for i in 0..(W_LOW_RES*H_LOW_RES + H_LOW_RES) {
        let c: u8 = POTENTIAL_ASCII[i];
        if c != b'\n' {
            potential_low_res.push(c);
        }
    }
    for i in 0..N { // Height
        for j in 0..N { // Width
            let d_i: usize = i/(N/H_LOW_RES);
            let d_j: usize = j/(N/W_LOW_RES);
            let c: u8 = potential_low_res[d_i*W_LOW_RES + d_j];
            let re_phi: f32 = if c == b'#' {
                (b'.' - c) as f32} else {0.0};
            let im_phi: f32 = if c == b'I' {
                (c - b'.') as f32} else {0.0};
            potential[(N - 1 - i)*N + j] = Complex {
                // real: 0.0*re_phi,
                real: 0.1*re_phi,
                imag: -im_phi,
            };
        }
    }
}

/* fn init_vector_potential(v_x: &mut [Complex<f32>],
                            v_y: &mut [Complex<f32>]) {
    for i in 0..N {
        for j in 0..N {
            v_x[i*N + j] = Complex {real: 0.0, imag: 0.0};
            v_y[i*N + j] = Complex {real: 0.0, imag: 0.0};
        }
    }

}*/

/* Initialize the square of the momentum values that correspond to the
real-space simulation domain. These are shifted to match the fft output. */
fn init_momentum_squared(p_squared: &mut [f32]) {
    for i in 0..N {
        for j in 0..N {
            let i_shift: i32 = if i < N/2 {i as i32} 
                else {-(N as i32) + (i as i32)};
            let j_shift: i32 = if j < N/2 {j as i32}
                else {-(N as i32) + (j as i32)};
            let px: f32 
                = 2.0*std::f32::consts::PI*(i_shift as f32)/(N as f32);
            let py: f32 
                = 2.0*std::f32::consts::PI*(j_shift as f32)/(N as f32);
            p_squared[i*N + j] = px*px + py*py;
        }
    }
}


/* Propagate the wave function psi with free-space periodic boundary
conditions for time step dt:
    |psi(dt)> = exp(-i*p_squared*dt/2)|psi(0)>.*/
fn propagate_kinetic(psi: &mut [Complex<f32>], 
                     p_squared: &[f32], dt: Complex<f32>,
                     use_mt: bool) {
    if use_mt {
        horizontal_square_fft(false, psi, N, TH_COUNT);
        square_transpose_in_place(psi, N);
        horizontal_square_fft(false, psi, N, TH_COUNT);
        square_transpose_in_place(psi, N);
    } else {
        for i in 0..N {
            fft_in_place(&mut psi[i*N..(i+1)*N], N);
        }
        square_transpose_in_place(psi, N);
        for i in 0..N {
            fft_in_place(&mut psi[i*N..(i+1)*N], N);
        }
        square_transpose_in_place(psi, N);
    }
    for i in 0..N {
        for j in 0..N {
            psi[i*N + j] = psi[i*N + j]*c64exp(
                Complex {real: 0.0, imag: -0.5*p_squared[i*N + j]} * dt);
        }
    }
    if use_mt {
        horizontal_square_fft(true, psi, N, TH_COUNT);
        square_transpose_in_place(psi, N);
        horizontal_square_fft(true, psi, N, TH_COUNT);
        square_transpose_in_place(psi, N);
    } else {
        for i in 0..N {
            ifft_in_place(&mut psi[i*N..(i+1)*N], N);
        }
        square_transpose_in_place(psi, N);
        for i in 0..N {
            ifft_in_place(&mut psi[i*N..(i+1)*N], N);
        }
        square_transpose_in_place(psi, N);
    }
}

struct Nonlinear {
    square: f32,
}

/* Apply the transformation 
    exp(-i*(potential + nonlinear(psi))*dt) on psi */
fn propagate_spatial_terms(
    psi: &mut [Complex<f32>], 
    potential: &[Complex<f32>],
    // vec_potential_x: &[Complex<f32>],
    // vec_potential_y: &[Complex<f32>],
    nonlinear: Nonlinear,
    dt: Complex<f32>) {
    for i in 0..N {
        for j in 0..N {
            let ij: usize = i*N + j;
            let psi_ij = psi[ij];
            let potential_ij = potential[ij];
            let nonlinear_term 
                = (psi_ij*psi_ij.conj()).scale(nonlinear.square);
            psi[i*N + j] = psi_ij*c64exp(
                Complex {real: 0.0, imag: -1.0} 
                * (potential_ij + nonlinear_term)* dt);
        }
    }
}

/* Dampen the wavefunction inside a region, where the probability
 * current inside this region is used to compute the decay. 
 *
 * References:
 *
 * Wikipedia - Probability current
 * https://en.wikipedia.org/wiki/Probability_current
 *
 * Widipedia - Perfectly matched layer
 * https://en.wikipedia.org/wiki/Perfectly_matched_layer
 */
fn dampen(psi: &mut [Complex<f32>], dt: f32) {
    let modi = |val: usize| {
        return val % N;
    };
    let mut jx = std::vec::Vec::<f32>::with_capacity(N*N);
    let mut jy = std::vec::Vec::<f32>::with_capacity(N*N);
    for i in 0..N { // height
        for j in 0..N { // width
            let y = (i as f32)/(N as f32); 
            let abs_psi2 = (psi[i*N + j]*psi[i*N + j].conj()).real;
            if y > 0.9 && abs_psi2 > 1e-30 {
                let ddx_psi = 
                    psi[N*i + modi(j+1)] - psi[N*i + modi(j)];
                let ddy_psi = 
                    psi[N*modi(i+1) + j] - psi[N*modi(i) + j];
                let val = 0.05 - f32::abs(y - 0.95);
                jx.push(val*(psi[i*N + j]*ddx_psi).imag);
                jy.push(val*(psi[i*N + j]*ddy_psi).imag);
            } else {
                jx.push(0.0);
                jy.push(0.0);
            }
        }
    }
    for i in 0..N {
        for j in 0..N {
            let damp_factor
                = f32::exp(-0.35*dt*f32::sqrt(jx[N*i + j]*jx[N*i + j]
                                + jy[N*i + j]*jy[N*i + j]));
            psi[N*i + j].real *= damp_factor;
            psi[N*i + j].imag *= damp_factor;
        }
    }
}

fn fill_pixel_data(pixels: &mut [u8], pixel_offset: usize,
    psi: & [Complex<f32>], psi_brightness: f64,
    phi: & [Complex<f32>], phi_brightness: f64,
    w: usize, h: usize) {
    for i in 0..h {
        for j in 0..w {
            let index: usize = i*w + j;
            let abs_val2: f64 = psi[index].length_squared() as f64;
            let c_psi: Color = argument_to_color(psi[index].arg() as f64);
            let phi_val: f64 = (phi[index].real as f64)*phi_brightness;
            let c = Color {
                r: phi_val + c_psi.r*abs_val2*psi_brightness,
                g: phi_val + c_psi.g*abs_val2*psi_brightness,
                b: phi_val + c_psi.b*abs_val2*psi_brightness};
            pixels[pixel_offset + 3*index + 2] = c.r as u8;
            pixels[pixel_offset + 3*index + 1] = c.g as u8;
            pixels[pixel_offset + 3*index] = c.b as u8;
        }
    }
}

fn main() {

    let mut boxed_pixels: Box<[u8; 54 + 3*N*N]> 
        = Box::new([0; 54 + 3*N*N]);
    let info: BitmapInfo = BitmapInfo {
        total_file_size: ((54 + 3*N*N) as i32),
        data_offset: 54,
        header_size: 40,
        width: N as i32,
        height: N as i32,
        plane_count: 1,
        bits_per_pixel: 24,
        compression_method: 0,
        image_size: ((3*N*N) as u32),
        horizontal_resolution: 100,
        vertical_resolution: 100,
        color_palette_count: 16777216,
        important_colors_count: 0,
    };
    fill_bitmap_header(&mut *boxed_pixels, info);
    
    let mut psi_vec 
        = std::vec::Vec::<Complex<f32>>::with_capacity(N*N);
    let mut potential_vec 
        = std::vec::Vec::<Complex<f32>>::with_capacity(N*N);
    // let mut vector_potential_x_vec
    //     = std::vec::Vec::<Complex<f32>>::with_capacity(N*N);
    // let mut vector_potential_y_vec
    //     = std::vec::Vec::<Complex<f32>>::with_capacity(N*N);
    let mut p_squared_vec
        = std::vec::Vec::<f32>::with_capacity(N*N);
    for _ in 0..N*N {
        psi_vec.push(Complex {real: 0.0, imag: 0.0});
        potential_vec.push(Complex {real: 0.0, imag: 0.0});
        // vector_potential_x_vec.push(Complex {real: 0.0, imag: 0.0});
        // vector_potential_y_vec.push(Complex {real: 0.0, imag: 0.0});
        p_squared_vec.push(0.0);
    }
    // On a device from around 2016:
    // The glsl/js version ran at 10 fps for 1024x1024.
    // This Rust implementation took 17m22.741s of real time
    // for 2400 frames and a grid size of 1024x1024,
    // corresponding to 2.3 fps.

    // On a device released late 2020:
    // The glsl/js version ran at 27-28 fps for 1024x1024.
    // Rust version took 5:48.50 minutes of cpu time for 3000 steps
    // and a grid size of 1024x1024, corresponding to
    // 9.13 steps/s.
    let dt = Complex {real: RE_DT, imag: IM_DT};

    // https://doc.rust-lang.org/book/
    //   ch12-01-accepting-command-line-arguments.html
    let mut input_args: Vec<String> = env::args().collect();
    if input_args.len() > 1 {
        let fname = input_args.pop();
        match load_f32_simulation_data(psi_vec.as_mut_slice(),
                                       potential_vec.as_mut_slice(),
                                       N as u32, N as u32,
                                       fname.unwrap()) {
            Ok(a) => a,
            Err(e) => println!("{}", e),
        };
    } else {
        init_wave_packet(psi_vec.as_mut_slice(), N, N,
                         WavePacket {a: 25.0, x0: 0.5, y0: 0.2,
                         sx: 0.07, sy: 0.07, 
                         nx: 0.0, 
                         // ny: 50.0*(N as f32)/512.0,
                         ny: 60.0,
                        });
        init_potential(potential_vec.as_mut_slice());
    }
    init_momentum_squared(p_squared_vec.as_mut_slice());

    for i in 0..NUMBER_OF_STEPS {
        propagate_spatial_terms(psi_vec.as_mut_slice(), 
                                potential_vec.as_slice(),
                                Nonlinear {square: 0.0},
                                dt.scale(0.5));
        propagate_kinetic(psi_vec.as_mut_slice(),
                          p_squared_vec.as_slice(), dt, true);
        dampen(psi_vec.as_mut_slice(), dt.real);
        propagate_spatial_terms(psi_vec.as_mut_slice(),
                                potential_vec.as_slice(),
                                Nonlinear {square: 0.0}, 
                                dt.scale(0.5));
        let at_every_step: usize = 3;
        if i % at_every_step == 0 {
            fill_pixel_data(&mut *boxed_pixels, 54,
                            psi_vec.as_slice(), 12.0, 
                            potential_vec.as_slice(), 100.0, N, N);
            let frame_number: usize = i/at_every_step;
            let prefix: String = String::from(SAVE_DIRECTORY);
            let number_str: String = if frame_number < 10 {
                "000".to_string() + &frame_number.to_string()
            } else if frame_number < 100 {
                "00".to_string() + &frame_number.to_string()
            } else if frame_number < 1000 {
                "0".to_string() + &frame_number.to_string()
            } else {
                frame_number.to_string()
            };
            let filename: String = prefix + &number_str + &".bmp".to_string();
            println!("Saving {}", filename);
            let _ = make_bitmap_file(filename, &mut *boxed_pixels);
        }
    }
    let _ = save_f32_simulation_data("last_state.bin".to_string(),
                                     psi_vec.as_slice(),
                                     potential_vec.as_slice(),
                                     N as u32, N as u32);
}
