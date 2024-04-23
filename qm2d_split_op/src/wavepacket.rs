use crate::complex::*;

pub struct WavePacket {
    pub a: f32, // amplitude
    // initial x and y positions (x: [0, 1], y: [0, 1])
    pub x0: f32, pub y0: f32,
    pub sx: f32, pub sy: f32, // x and y standard deviations
    pub nx: f32, pub ny: f32, // Wavenumber in the x and y direction
}

pub fn init_wave_packet(
    array: &mut [Complex<f32>],
    width: usize, height: usize,
    w: WavePacket) {
    for i in 0..height {
        for j in 0..width  {
            let x: f32 = (j as f32)/(width as f32);
            let y: f32 = (i as f32)/(height as f32);
            let xt: f32 = x - w.x0;
            let yt: f32 = y - w.y0;
            let abs_val: f32 = w.a
                *f32::exp(-0.5*xt*xt/(w.sx*w.sx))
                *f32::exp(-0.5*yt*yt/(w.sy*w.sy));
            let nr = w.nx*x + w.ny*y;
            array[i*width + j] = Complex {
                real: abs_val*f32::cos(2.0*std::f32::consts::PI*nr),
                imag: abs_val*f32::sin(2.0*std::f32::consts::PI*nr),
            };

        }
    }
}
