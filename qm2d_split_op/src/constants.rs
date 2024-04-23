/* Constants used for the simulation. Change these before rebuilding
the project.
*/

// Simulation side length in number of pixels (currently only squares used).
// THIS NEEDS TO BE A POWER OF TWO, or else bad things may happen.
pub const N: usize = 1024;

pub const NUMBER_OF_STEPS: usize = 3000;

// The time step used. This is a complex value.
pub const RE_DT: f32 = 0.5;
pub const IM_DT: f32 = 0.0;

// Number of threads to use in threaded computations
pub const TH_COUNT: usize = 8;

// Where to save output images
pub const SAVE_DIRECTORY: &str = "./";

pub const W_LOW_RES: usize = 32;
pub const H_LOW_RES: usize = 32;

// ASCII representation of the potential.
pub const POTENTIAL_ASCII: [u8; W_LOW_RES*H_LOW_RES + H_LOW_RES] = *b"
................................
................................
................................
................................
................................
................................
................................
................................
................................
................................
................................
................................
................................
................................
................................
................................
................................
................................
##############.##.##############
................................
................................
................................
................................
................................
................................
................................
................................
................................
................................
................................
................................
................................";
