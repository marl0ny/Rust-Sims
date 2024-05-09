# Quantum Mechanics in 2D with the Split Operator method

Numerically solve for the time-dependent Schrodinger equation in 2D,
using the split operator method. No dependencies required, just type
`cargo build` and `cargo run`.
This will output a series of bmp images where each image is a frame of the 
simulation. Note that the file `constants.rs` contains compile time constants
that control the parameters of the simulation. While these remain constant when the simulation is running, change these values beforehand and rebuild to modify the simulation output.

This project is a refactor of [this gist](https://gist.github.com/marl0ny/81a2e5498a05f50040f4d928ad805ef6), where everything was originally contained in a single massive file. Special thanks to the commentators of [this discussion](https://www.reddit.com/r/rust/comments/1c8tsig/comment/l0guur7), especially VorpalWay's and LiesAreFunny's helpful suggestions, as well as gaolaowai's [pull request](https://github.com/marl0ny/Rust-Sims/pull/1) for properly fixing file organization.

## References:

### Split-Operator Method:

 - James Schloss. [The Split Operator Method - Arcane Algorithm Archive.](https://www.algorithm-archive.org/contents/split-operator_method/split-operator_method.html)

### Fast Fourier Transform (Used in the Split-Operator method):

 - [Wikipedia - Cooley–Tukey FFT algorithm](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm)

 - [MathWorld Wolfram - Fast Fourier Transform](http://mathworld.wolfram.com/FastFourierTransform.html)

 - William Press et al. [12.2 Fast Fourier Transform (FFT) - in Numerical Recipes](https://websites.pmc.ucsc.edu/~fnimmo/eart290c_17/NumericalRecipesinF77.pdf)

### Domain coloring method for visualizing complex-valued functions:

 - [Wikipedia - Domain coloring](https://en.wikipedia.org/wiki/Domain_coloring)

 - [Wikipedia - Hue](https://en.wikipedia.org/wiki/Hue)

 - [https://en.wikipedia.org/wiki/Hue#/media/File:HSV-RGB-comparison.svg](https://en.wikipedia.org/wiki/Hue#/media/File:HSV-RGB-comparison.svg)

### Bitmap file format:

 - [Wikipedia - BMP file format](https://en.wikipedia.org/wiki/BMP_file_format)