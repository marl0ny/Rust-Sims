pub struct Color {
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

/* Function that converts a hue angle to its corresponding color.

References:

Wikipedia - Domain coloring
https://en.wikipedia.org/wiki/Domain_coloring

Wikipedia - Hue
https://en.wikipedia.org/wiki/Hue

https://en.wikipedia.org/wiki/Hue#/media/File:HSV-RGB-comparison.svg

 */
pub fn argument_to_color(arg_val: f64) -> Color {
    let pi: f64 = std::f64::consts::PI;
    let max_col: f64 = 1.0;
    let min_col: f64 = 50.0/255.0;
    let col_range: f64 = max_col - min_col;
    if arg_val <= pi/3.0 && arg_val >= 0.0 {
        return Color {
            r: max_col,
            g: min_col + col_range*arg_val/(pi/3.0), 
            b: min_col};
    } else if arg_val > pi/3.0 && arg_val <= 2.0*pi/3.0 {
        return Color {
            r: max_col - col_range*(arg_val - pi/3.0)/(pi/3.0),
            g: max_col, 
            b: min_col};
    } else if arg_val > 2.0*pi/3.0 && arg_val <= pi {
        return Color {
            r: min_col, 
            g: max_col,
            b: min_col + col_range*(arg_val - 2.0*pi/3.0)/(pi/3.0)};
    } else if arg_val < 0.0 && arg_val > -pi/3.0 {
        return Color {
            r: max_col, 
            g: min_col,
            b: min_col - col_range*arg_val/(pi/3.0)};
    } else if arg_val <= -pi/3.0 && arg_val > -2.0*pi/3.0 {
        return Color {
            r: max_col + (col_range*(arg_val + pi/3.0)/(pi/3.0)),
            g: min_col, 
            b: max_col};
    } else if arg_val <= -2.0*pi/3.0 && arg_val >= -pi {
        return Color {
            r: min_col,
            g: min_col - (col_range*(arg_val + 2.0*pi/3.0)/(pi/3.0)),
            b: max_col};
    }
    else {
        return Color {r: min_col, g: max_col, b: max_col};
    }
}