pub mod constants {

    // Simulation side length in number of pixels (currently only squares used)
    pub const N: usize = 1024;

    pub const NUMBER_OF_STEPS: usize = 3000;
    // The timestep used. This is a complex value.
    pub const RE_DT: f32 = 0.5;
    pub const IM_DT: f32 = 0.0;

    // Where to save output images
    pub const SAVE_DIRECTORY: &str = "./";

    // Number of threads to use for threaded computations
    pub const TH_COUNT: usize = 8;

    pub const W_LOW_RES: usize = 32;
    pub const H_LOW_RES: usize = 32;

}
