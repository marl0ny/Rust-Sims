
#[derive(Copy, Clone)]
pub struct Complex<T> {
    pub real: T,
    pub imag: T,
}

impl <T: std::ops::Neg<Output=T>> Complex<T> {
    pub fn conj(self) -> Complex<T> {
        return Complex{real: self.real, imag: -self.imag};
    }
}

/* Into conversion from the size 64 (2 x f32) to size 128 (2 x f64) 
bit complex struct. This follows closely to the example given in 
the Rust documentation:

https://doc.rust-lang.org/rust-by-example/conversion/from_into.html
*/
impl std::convert::Into<Complex<f64>> for Complex<f32> {
    fn into(self) -> Complex<f64> {
        return Complex{real: self.real as f64, imag: self.imag as f64};
    }
}

/* Into conversion from the size 128 (2 x f64) to size 64 bit (2 x f32) 
complex struct. This follows closely to the example given in 
the Rust documentation:

https://doc.rust-lang.org/rust-by-example/conversion/from_into.html
*/
impl std::convert::Into<Complex<f32>> for Complex<f64> {
    fn into(self) -> Complex<f32> {
        return Complex{real: self.real as f32, imag: self.imag as f32};
    }
}

impl Complex<f32> {
    pub fn arg(self) -> f64 {
        if self.real == 0.0 {
            if self.imag >= 0.0 {
                return std::f64::consts::PI/2.0;
            } else {
                return -std::f64::consts::PI/2.0;
            }
        } else {
            let val: f64 = f64::atan((self.imag/self.real).into());
            if self.real < 0.0 {
                if self.imag >= 0.0 {
                    return std::f64::consts::PI + val;
                } else {
                    return -std::f64::consts::PI + val;
                }
            }
            return val;
        }
    }
} 

impl <T: std::ops::Neg<Output=T>
        + std::ops::Add<Output=T> 
        + std::ops::Mul<Output=T> 
        + std::ops::Div<Output=T> + Copy> Complex<T> {

    pub fn length_squared(self) -> T {
        return self.real*self.real + self.imag*self.imag;
    }

    pub fn inv(self) -> Complex<T> {
        return Complex {real: self.real/self.length_squared(),
                        imag: -self.imag/self.length_squared()};
    }

}

/* Operator add overloading for the Complex struct.
This and the other operator overloading functions are based 
on the code examples from the Rust documentation:

https://doc.rust-lang.org/std/ops/
*/
impl <T: std::ops::Add<Output=T>> std::ops::Add for Complex<T> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        return Self {real: self.real + other.real, 
                    imag: self.imag + other.imag};
    }
}

impl <T: std::ops::Sub<Output=T>> std::ops::Sub for Complex<T> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        return Self {real: self.real - other.real, 
                    imag: self.imag - other.imag};
    }
}

impl <T: std::ops::Mul<Output=T> 
        + std::ops::Add<Output=T> 
        + std::ops::Sub<Output=T>
        + Copy> std::ops::Mul for Complex<T> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        return Self {real: self.real*other.real - self.imag*other.imag,
                    imag: self.real*other.imag + self.imag*other.real};
    }
}

impl <T: std::ops::Neg<Output=T>
        + std::ops::Add<Output=T> 
        + std::ops::Sub<Output=T> 
        + std::ops::Mul<Output=T> 
        + std::ops::Div<Output=T> + Copy>std::ops::Div for Complex<T> {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        return self*other.inv();
    }
}

/* This also follows the example from
https://doc.rust-lang.org/std/ops/trait.Mul.html
*/
impl <T: std::ops::Neg<Output=T>
        + std::ops::Add<Output=T> 
        + std::ops::Sub<Output=T> 
        + std::ops::Mul<Output=T> 
        + std::ops::Div<Output=T> + Copy>std::ops::Mul<T> for Complex<T> {
    type Output = Self;
    fn mul(self, other: T) -> Self {
        return Complex {
            real: self.real*other,
            imag: self.imag*other
        };
    }
}

/* This also follows the example from
https://doc.rust-lang.org/std/ops/trait.Div.html 
*/
impl <T: std::ops::Neg<Output=T>
        + std::ops::Add<Output=T> 
        + std::ops::Sub<Output=T> 
        + std::ops::Mul<Output=T> 
        + std::ops::Div<Output=T> + Copy>std::ops::Div<T> for Complex<T> {
    type Output = Self;
    fn div(self, other: T) -> Self {
        return Complex {
            real: self.real/other,
            imag: self.imag/other
        };
    }
}

// Perhaps the next following two functions can
// be refactored in to one?
impl std::ops::Mul<Complex<f64>> for f64 {
    type Output = Complex<f64>;
    fn mul(self, other: Complex<f64>) -> Self::Output {
        return Complex {
            real: other.real*self,
            imag: other.imag*self
        };
    }
}

impl std::ops::Mul<Complex<f32>> for f64 {
    type Output = Complex<f32>;
    fn mul(self, other: Complex<f32>) -> Self::Output {
        return Complex {
            real: ((other.real as f64)*self) as f32,
            imag: ((other.imag as f64)*self) as f32
        };
    }
}

impl Complex<f64> {
    pub fn exp(self) -> Complex<f64> {
        return Complex {
            real: f64::exp(self.real)*f64::cos(self.imag),
            imag: f64::exp(self.real)*f64::sin(self.imag),
        };
    }
}

impl Complex<f32> {
    pub fn exp(self) -> Complex<f32> {
        return Complex {
            real: (f64::exp(self.real as f64)*f64::cos(self.imag as f64))
                 as f32,
            imag: (f64::exp(self.real as f64)*f64::sin(self.imag as f64))
                 as f32,
        };
    }
}
