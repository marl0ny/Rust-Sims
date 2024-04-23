/* Serialize simulation data.
*/
use crate::complex::*;
use crate::copy_to_bytes::*;

pub fn save_f32_simulation_data(filename: std::string::String,
             psi: &[Complex<f32>],
             potential: &[Complex<f32>],
             simulation_width: u32,
             simulation_height: u32
             ) -> std::io::Result<()> {
    let sizeof_complex: usize = 2*4;
    let header_size: usize = 8;
    let width: u32 = simulation_width;
    let height: u32 = simulation_height;
    let total_size: usize 
        = header_size + 2*sizeof_complex*((width*height) as usize);
    let mut bytes = std::vec::Vec::<u8>::with_capacity(total_size);
    for _ in 0..total_size {
        bytes.push(0);
    };
    let mut offset: usize = 0;
    copy_u32(bytes.as_mut_slice(), width, &mut offset);
    copy_u32(bytes.as_mut_slice(), height, &mut offset);
    for i in 0..width*height {
        copy_f32(bytes.as_mut_slice(), psi[i as usize].real, &mut offset);
        copy_f32(bytes.as_mut_slice(), psi[i as usize].imag, &mut offset);
        copy_f32(bytes.as_mut_slice(), 
                 potential[i as usize].real, &mut offset);
        copy_f32(bytes.as_mut_slice(), 
                 potential[i as usize].imag, &mut offset);
    }
    use std::io::Write;
    let mut file = std::fs::File::create(filename)?;
    file.write_all(bytes.as_slice())?;
    Ok(())
}


pub fn load_f32_simulation_data(
    psi: &mut [Complex<f32>],
    potential: &mut [Complex<f32>],
    simulation_width: u32,
    simulation_height: u32,
    filename: std::string::String,
    ) -> std::io::Result<()> {
    let sizeof_complex: usize = 2*4;
    let header_size: usize = 8;
    let width: u32 = simulation_width;
    let height: u32 = simulation_height;
    let total_size: usize 
        = header_size + 2*sizeof_complex*((width*height) as usize);
    use std::io::prelude::*;
    let mut bytes = std::vec::Vec::<u8>::with_capacity(total_size);
    let file = std::io::BufReader::new(std::fs::File::open(filename)?);
    let mut size = 0;
    // https://doc.rust-lang.org/std/io/trait.Read.html#method.bytes
    for b in file.bytes() {
        // https://doc.rust-lang.org/
        // rust-by-example/error/result/early_returns.html
        match b {
            Ok(val) => bytes.push(val),
            Err(e) => return Err(e),
        };
        size += 1;
        if size > total_size {
            break; // TODO!
            // return Err::<&str, ()>(());
        }
    }
    let width2: u32 = write_u32(&bytes.as_slice()[0..4]);
    let height2: u32 = write_u32(&bytes.as_slice()[4..8]);
    // println!("{}, {}", width2, height2);
    if width2 == width && height2 == height {
        let arr = &bytes.as_slice()[8..total_size];
        let s = sizeof_complex*2;
        for i in 0..width*height {
            let k = s*(i as usize);
            psi[i as usize] = Complex {
                real: write_f32(&arr[k..k+4]),
                imag: write_f32(&arr[k+4..k+8])};
            potential[i as usize] = Complex {
                real: write_f32(&arr[k+8..k+12]), 
                imag: write_f32(&arr[k+12..k+16])};
        }
    } else {
        // TODO
    }
    Ok(())

}

fn write_f32(bytes: &[u8]) -> f32 {
    let mut val_arr: [u8; 4] = [0; 4];
    for i in 0..4 {
        val_arr[i] = bytes[i];
    }
    return f32::from_ne_bytes(val_arr);
}

fn write_u32(bytes: &[u8]) -> u32 {
    let mut val_arr: [u8; 4] = [0; 4];
    for i in 0..4 {
        val_arr[i] = bytes[i];
    }
    return u32::from_ne_bytes(val_arr);
}

